#!/usr/bin/env python3
import pysam
import argparse
import subprocess
import sys
import os


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract BrdU modification data from a BAM file in a specific region."
    )
    parser.add_argument("bam")
    parser.add_argument("-c", "--chrom", required=True)
    parser.add_argument("-s", "--start", required=True, type=int)
    parser.add_argument("-e", "--end", required=True, type=int)
    parser.add_argument("-r", "--read_id", default=None)
    parser.add_argument("-o", "--output", default=None)
    return parser.parse_args()


def get_project_paths():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    sorted_bam_dir = os.path.abspath(os.path.join(script_dir, "../../data/sorted_bam"))
    index_dir = os.path.abspath(os.path.join(script_dir, "../../data/index_sorted_bam_bai"))
    return sorted_bam_dir, index_dir


def check_and_sort_bam(bam_path, sorted_bam_dir):
    bam = pysam.AlignmentFile(bam_path, "rb", check_sq=False)
    header = bam.header.to_dict()
    bam.close()

    sort_order = header.get("HD", {}).get("SO", "unknown")

    os.makedirs(sorted_bam_dir, exist_ok=True)

    original_base = os.path.splitext(os.path.basename(bam_path))[0]

    if original_base.endswith(".sorted"):
        original_base = original_base[:-7]

    standardized_filename = f"{original_base}.sorted.indexed.bam"
    standardized_path = os.path.join(sorted_bam_dir, standardized_filename)

    if os.path.exists(standardized_path):
        print(f"[INFO] Reusing standardized BAM: {standardized_path}", file=sys.stderr)
        return standardized_path

    if sort_order == "coordinate":
        print(f"[INFO] Copying sorted BAM → {standardized_path}", file=sys.stderr)
        subprocess.run(["cp", bam_path, standardized_path], check=True)
        return standardized_path

    print(f"[INFO] Sorting BAM → {standardized_path}", file=sys.stderr)
    subprocess.run(["samtools", "sort", "-o", standardized_path, bam_path], check=True)
    return standardized_path


def check_and_index_bam(bam_path, index_dir):
    bam_basename = os.path.basename(bam_path)
    index_path = os.path.join(index_dir, f"{bam_basename}.bai")

    os.makedirs(index_dir, exist_ok=True)

    rebuild = False

    if not os.path.exists(index_path):
        print("[INFO] No index found. Creating index...", file=sys.stderr)
        rebuild = True
    else:
        if os.path.getmtime(index_path) < os.path.getmtime(bam_path):
            print("[INFO] Index older than BAM. Rebuilding index...", file=sys.stderr)
            rebuild = True
        else:
            print(f"[INFO] Index found: {index_path}", file=sys.stderr)
            return index_path

    subprocess.run(["samtools", "index", "-o", index_path, bam_path], check=True)
    print(f"[INFO] Indexing complete: {index_path}", file=sys.stderr)
    return index_path


def extract_brdu(bam_path, chrom, start, end, read_id_filter=None):
    results = []

    _, index_dir = get_project_paths()
    index_path = os.path.join(index_dir, f"{os.path.basename(bam_path)}.bai")

    bam = pysam.AlignmentFile(bam_path, "rb", index_filename=index_path)

    for read in bam.fetch(chrom, start, end):
        if read.is_unmapped:
            continue

        if read_id_filter and read_id_filter not in read.query_name:
            continue

        if not read.modified_bases:
            continue

        for (_, _, mod_code), mod_list in read.modified_bases.items():
            if mod_code != ord('b') and mod_code != 'b':
                continue

            aligned_pairs = dict(read.get_aligned_pairs(matches_only=True))

            for query_pos, raw_prob in mod_list:
                ref_pos = aligned_pairs.get(query_pos)
                if ref_pos is None or not (start <= ref_pos < end):
                    continue

                base = read.query_sequence[query_pos]
                if base != "T":
                    continue

                results.append((
                    chrom,
                    ref_pos,
                    ref_pos + 1,
                    read.query_name,
                    base,
                    raw_prob / 255.0
                ))

    bam.close()
    results.sort(key=lambda x: (x[3], x[1]))
    return results


def main():
    args = parse_args()

    if not os.path.exists(args.bam):
        print(f"[ERROR] BAM not found: {args.bam}", file=sys.stderr)
        sys.exit(1)

    sorted_dir, index_dir = get_project_paths()

    bam_path = check_and_sort_bam(args.bam, sorted_dir)
    check_and_index_bam(bam_path, index_dir)

    results = extract_brdu(bam_path, args.chrom, args.start, args.end, args.read_id)

    out = open(args.output, "w") if args.output else sys.stdout
    for r in results:
        out.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\t{r[5]:.7f}\n")

    if args.output:
        out.close()
        print(f"[INFO] Results written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
