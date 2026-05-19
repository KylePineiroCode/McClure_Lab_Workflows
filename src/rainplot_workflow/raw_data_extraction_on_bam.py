#!/usr/bin/env python3

import pysam
import argparse
import subprocess
import sys
import os


def parse_args():
    """
    Parse command-line arguments for extracting BrdU modification data.

    The script takes a BAM file, a genomic region, and an optional read ID filter.
    It outputs per-base BrdU modification calls in a BED-like format.

    Returns
    -------
    argparse.Namespace
        Parsed arguments containing the BAM path, chromosome, start coordinate,
        end coordinate, optional read ID, and optional output file path.
    """
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
    """
    Build paths to workflow-managed BAM output directories.

    The sorted BAMs and BAM indexes are stored in fixed workflow directories so
    that downstream scripts can reuse them instead of repeatedly sorting and
    indexing the same BAM file.

    Returns
    -------
    tuple[str, str]
        A tuple containing:
        - sorted_bam_dir: directory for standardized sorted BAM files
        - index_dir: directory for BAM index files
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    sorted_bam_dir = os.path.abspath(os.path.join(script_dir, "../../data/sorted_bam"))
    index_dir = os.path.abspath(os.path.join(script_dir, "../../data/index_sorted_bam_bai"))
    return sorted_bam_dir, index_dir


def check_and_sort_bam(bam_path, sorted_bam_dir):
    """
    Ensure the BAM file is coordinate-sorted and saved with a standardized name.

    Region-based BAM fetching requires a coordinate-sorted and indexed BAM. This
    function checks the BAM header to see whether the input is already sorted.
    If it is already coordinate-sorted, the BAM is copied into the workflow's
    sorted BAM directory. If it is not sorted, samtools sort is used to create a
    sorted copy.

    The output name is standardized as:
        <original_basename>.sorted.indexed.bam

    This makes later indexing and reuse predictable.

    Parameters
    ----------
    bam_path : str
        Path to the input BAM file.

    sorted_bam_dir : str
        Directory where the sorted/standardized BAM should be stored.

    Returns
    -------
    str
        Path to the standardized coordinate-sorted BAM file.
    """
    bam = pysam.AlignmentFile(bam_path, "rb", check_sq=False)
    header = bam.header.to_dict()
    bam.close()

    sort_order = header.get("HD", {}).get("SO", "unknown")

    os.makedirs(sorted_bam_dir, exist_ok=True)

    original_base = os.path.splitext(os.path.basename(bam_path))[0]

    # Avoid creating names like sample.sorted.sorted.indexed.bam when the input
    # BAM filename already ends with ".sorted".
    if original_base.endswith(".sorted"):
        original_base = original_base[:-7]

    standardized_filename = f"{original_base}.sorted.indexed.bam"
    standardized_path = os.path.join(sorted_bam_dir, standardized_filename)

    # Reuse an existing standardized BAM to save time on repeated workflow runs.
    if os.path.exists(standardized_path):
        print(f"[INFO] Reusing standardized BAM: {standardized_path}", file=sys.stderr)
        return standardized_path

    # If the BAM header says it is already coordinate-sorted, copying is faster
    # than running samtools sort again.
    if sort_order == "coordinate":
        print(f"[INFO] Copying sorted BAM → {standardized_path}", file=sys.stderr)
        subprocess.run(["cp", bam_path, standardized_path], check=True)
        return standardized_path

    print(f"[INFO] Sorting BAM → {standardized_path}", file=sys.stderr)
    subprocess.run(["samtools", "sort", "-o", standardized_path, bam_path], check=True)
    return standardized_path


def check_and_index_bam(bam_path, index_dir):
    """
    Ensure the sorted BAM has a current BAM index.

    pysam needs a BAM index to fetch reads from a specific genomic region. This
    function creates the index if it does not exist, or rebuilds it if the BAM
    file is newer than the index.

    The index is written to the workflow index directory rather than next to the
    BAM file so all generated indexes are kept in one place.

    Parameters
    ----------
    bam_path : str
        Path to the coordinate-sorted BAM file.

    index_dir : str
        Directory where the BAM index should be stored.

    Returns
    -------
    str
        Path to the BAM index file.
    """
    bam_basename = os.path.basename(bam_path)
    index_path = os.path.join(index_dir, f"{bam_basename}.bai")

    os.makedirs(index_dir, exist_ok=True)

    rebuild = False

    if not os.path.exists(index_path):
        print("[INFO] No index found. Creating index...", file=sys.stderr)
        rebuild = True
    else:
        # Rebuild the index if the BAM has been modified more recently than the
        # index. This prevents using a stale index with a newer BAM.
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
    """
    Extract BrdU modification calls from reads overlapping a genomic region.

    This function reads a sorted/indexed BAM file and searches for modified base
    calls with modification code "b", which represents BrdU in this workflow.
    It maps each modified query base back to its reference coordinate and keeps
    only BrdU calls that fall inside the requested region.

    Only thymine bases are kept because BrdU is interpreted as a thymidine analog
    in this pipeline.

    Parameters
    ----------
    bam_path : str
        Path to the sorted and indexed BAM file.

    chrom : str
        Chromosome or contig name to fetch from the BAM.

    start : int
        Start coordinate of the region, using 0-based BED-style coordinates.

    end : int
        End coordinate of the region, using a half-open interval.

    read_id_filter : str, optional
        Optional substring used to restrict extraction to a specific read ID.

    Returns
    -------
    list[tuple]
        A sorted list of BrdU calls. Each tuple contains:
        chrom, start, end, read_id, base, normalized_probability.
    """
    results = []

    _, index_dir = get_project_paths()
    index_path = os.path.join(index_dir, f"{os.path.basename(bam_path)}.bai")

    bam = pysam.AlignmentFile(bam_path, "rb", index_filename=index_path)

    for read in bam.fetch(chrom, start, end):
        if read.is_unmapped:
            continue

        # Allow partial matching so users can provide either the full read ID or
        # a unique substring from the read name.
        if read_id_filter and read_id_filter not in read.query_name:
            continue

        if not read.modified_bases:
            continue

        for (_, _, mod_code), mod_list in read.modified_bases.items():
            # DNAscent/modBAM BrdU calls are represented with modification code
            # "b". Depending on pysam/version behavior, the code may appear as
            # the character "b" or its ASCII integer value.
            if mod_code != ord("b") and mod_code != "b":
                continue

            # modified_bases reports positions in read/query coordinates.
            # To plot them on the genome, convert query positions to reference
            # positions using the read alignment.
            aligned_pairs = dict(read.get_aligned_pairs(matches_only=True))

            for query_pos, raw_prob in mod_list:
                ref_pos = aligned_pairs.get(query_pos)

                if ref_pos is None or not (start <= ref_pos < end):
                    continue

                base = read.query_sequence[query_pos]

                # BrdU replaces thymidine, so this keeps the output focused on
                # biologically relevant T positions instead of all modified calls.
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

    # Sort by read ID and genomic position so each read's BrdU calls are grouped
    # together in a stable order for downstream plotting.
    results.sort(key=lambda x: (x[3], x[1]))
    return results


def main():
    """
    Run the BrdU extraction workflow.

    This function validates the input BAM, prepares a sorted and indexed BAM,
    extracts BrdU calls for the requested region, and writes the results either
    to a user-provided output file or to standard output.

    Output columns are:
        chrom, start, end, read_id, base, BrdU_probability
    """
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
        # Write BED-like rows:
        # chrom, 0-based start, 1-base end, read ID, base, normalized BrdU score.
        out.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\t{r[5]:.7f}\n")

    if args.output:
        out.close()
        print(f"[INFO] Results written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()