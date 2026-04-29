#!/usr/bin/env python3
import pysam
import argparse
import subprocess
import sys
import os
from difflib import SequenceMatcher


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("bam")
    parser.add_argument("-c", "--chrom", required=True)
    parser.add_argument("-s", "--start", required=True, type=int)
    parser.add_argument("-e", "--end", required=True, type=int)
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

    standardized_path = os.path.join(sorted_bam_dir, f"{original_base}.sorted.indexed.bam")

    if os.path.exists(standardized_path):
        print(f"[INFO] Reusing standardized BAM: {standardized_path}", file=sys.stderr)
        return standardized_path

    if sort_order == "coordinate":
        subprocess.run(["cp", bam_path, standardized_path], check=True)
        return standardized_path

    subprocess.run(["samtools", "sort", "-o", standardized_path, bam_path], check=True)
    return standardized_path


def check_and_index_bam(bam_path, index_dir):
    index_path = os.path.join(index_dir, f"{os.path.basename(bam_path)}.bai")
    os.makedirs(index_dir, exist_ok=True)

    rebuild = (
        not os.path.exists(index_path)
        or os.path.getmtime(index_path) < os.path.getmtime(bam_path)
    )

    if rebuild:
        print("[INFO] Creating/Rebuilding index...", file=sys.stderr)
        subprocess.run(["samtools", "index", "-o", index_path, bam_path], check=True)
    else:
        print(f"[INFO] Index found: {index_path}", file=sys.stderr)

    return index_path


def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()


def extract_rfb_reads(bam_path, chrom, start, end):
    motif = "TTTACCAAGAAAGATGTAAG"
    threshold = (len(motif) - 2) / len(motif)

    _, index_dir = get_project_paths()
    index_path = os.path.join(index_dir, f"{os.path.basename(bam_path)}.bai")

    bam = pysam.AlignmentFile(bam_path, "rb", index_filename=index_path)

    results = []

    for read in bam.fetch(chrom, start, end):
        if read.is_unmapped:
            continue

        seq = read.query_sequence
        if not seq:
            continue

        for i in range(len(seq) - len(motif)):
            if similar(seq[i:i+len(motif)], motif) >= threshold:
                results.append((chrom, read.reference_start, read.reference_end, read.query_name, "T", 1.0))
                break

    bam.close()
    return results


def main():
    args = parse_args()

    sorted_dir, index_dir = get_project_paths()

    bam_path = check_and_sort_bam(args.bam, sorted_dir)
    check_and_index_bam(bam_path, index_dir)

    results = extract_rfb_reads(bam_path, args.chrom, args.start, args.end)

    print(f"[INFO] {len(results)} reads containing RFB motif found.", file=sys.stderr)

    if len(results) == 0:
        print("[INFO] No RFB found — expected for non S-phase datasets.", file=sys.stderr)
        return

    out = open(args.output, "w") if args.output else sys.stdout
    for r in results:
        out.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\t{r[5]:.7f}\n")

    if args.output:
        out.close()


if __name__ == "__main__":
    main()
