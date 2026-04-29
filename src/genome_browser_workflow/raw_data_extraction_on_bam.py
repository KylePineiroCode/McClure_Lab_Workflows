#!/usr/bin/env python3
import pysam
import argparse
import subprocess
import sys
import os
import tempfile
import shutil
from collections import defaultdict


def parse_args():
    """
    Parse command line arguments for genome browser BrdU extraction.

    Required:
        bam: Path to input BAM file

    Optional:
        -o / --output: Output prefix for generated bedgraph files

    Returns:
        argparse.Namespace with parsed arguments
    """
    parser = argparse.ArgumentParser(
        description="Extract genome-wide BrdU bedgraph data from a BAM file using modkit."
    )
    parser.add_argument("bam")
    parser.add_argument(
        "-o", "--output", default=None,
        help="Optional output prefix for generated bedgraph files."
    )
    return parser.parse_args()


def get_project_paths():
    """
    Resolve key project directories relative to this script.

    Returns:
        sorted_bam_dir: Directory where sorted BAM files are stored
        index_dir: Directory where archived BAM index (.bai) files are stored
        bedgraph_dir: Directory where output bedgraph files are saved
        bam_dir: Directory where original BAM files are stored
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))

    sorted_bam_dir = os.path.abspath(
        os.path.join(script_dir, "../../data/sorted_bam")
    )
    index_dir = os.path.abspath(
        os.path.join(script_dir, "../../data/index_sorted_bam_bai")
    )
    bedgraph_dir = os.path.abspath(
        os.path.join(script_dir, "../../data/bedgraph")
    )
    bam_dir = os.path.abspath(
        os.path.join(script_dir, "../../data/bam")
    )

    return sorted_bam_dir, index_dir, bedgraph_dir, bam_dir


def resolve_bam_path(bam_input, bam_dir):
    """
    Resolve the BAM path provided by the user.

    If the BAM path exists as provided, use it directly.
    Otherwise, look for the BAM file in data/bam.

    Args:
        bam_input: User-provided BAM path or filename
        bam_dir: Default BAM directory in the workflow

    Returns:
        Resolved BAM path

    Raises:
        SystemExit if the BAM cannot be found
    """
    if os.path.exists(bam_input):
        return os.path.abspath(bam_input)

    candidate = os.path.join(bam_dir, bam_input)
    if os.path.exists(candidate):
        return os.path.abspath(candidate)

    print(f"[ERROR] BAM not found: {bam_input}", file=sys.stderr)
    sys.exit(1)


def check_and_sort_bam(bam_path, sorted_bam_dir):
    """
    Check if a BAM file is coordinate-sorted. If not, sort it using samtools.

    If the BAM is already sorted, it is copied into the standardized location.
    All BAM files are renamed to a consistent format:
        <original>.sorted.indexed.bam

    Args:
        bam_path: Path to input BAM file
        sorted_bam_dir: Directory to store sorted BAM files

    Returns:
        Path to sorted BAM file
    """
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
    """
    Ensure a BAM file has an up-to-date index.

    The archived index is stored in:
        data/index_sorted_bam_bai/

    The archived index filename is:
        <bam_filename>.bai

    An adjacent BAM index is also created next to the BAM so tools like modkit
    can find it automatically.

    If the archived index is missing or older than the BAM file, it is rebuilt.

    Args:
        bam_path: Path to BAM file
        index_dir: Directory where archived index files are stored

    Returns:
        Path to archived index file
    """
    bam_basename = os.path.basename(bam_path)
    archived_index_path = os.path.join(index_dir, f"{bam_basename}.bai")
    adjacent_index_path = f"{bam_path}.bai"

    os.makedirs(index_dir, exist_ok=True)

    rebuild = False

    if not os.path.exists(archived_index_path):
        print("[INFO] No index found. Creating index...", file=sys.stderr)
        rebuild = True
    else:
        if os.path.getmtime(archived_index_path) < os.path.getmtime(bam_path):
            print("[INFO] Index older than BAM. Rebuilding index...", file=sys.stderr)
            rebuild = True
        else:
            print(f"[INFO] Index found: {archived_index_path}", file=sys.stderr)

    if rebuild:
        subprocess.run(["samtools", "index", bam_path], check=True)
        print(f"[INFO] Adjacent BAM index created: {adjacent_index_path}", file=sys.stderr)

        shutil.copy2(adjacent_index_path, archived_index_path)
        print(f"[INFO] Archived index created: {archived_index_path}", file=sys.stderr)
    else:
        if not os.path.exists(adjacent_index_path):
            shutil.copy2(archived_index_path, adjacent_index_path)
            print(f"[INFO] Restored adjacent BAM index: {adjacent_index_path}", file=sys.stderr)

    return archived_index_path


def extract_brdu_with_modkit(bam_path):
    """
    Run modkit extract full on the BAM file to obtain per-base modification data.

    Uses tempfile.mkstemp to generate a safe, unique output path, then removes
    the placeholder so modkit can write to that path (modkit refuses to overwrite
    existing files).

    Args:
        bam_path: Path to sorted and indexed BAM file

    Returns:
        Path to the modkit output TSV file (caller is responsible for deletion)
    """
    fd, temp_output = tempfile.mkstemp(suffix=".tsv", prefix="modkit_extract_")
    os.close(fd)
    os.remove(temp_output)  # modkit refuses to overwrite existing files

    cmd = ["modkit", "extract", "full", bam_path, temp_output]

    print(f"[INFO] Running modkit extraction on {bam_path}", file=sys.stderr)

    # stderr=None streams modkit's progress and errors live to the terminal
    result = subprocess.run(cmd, stderr=None)

    if result.returncode != 0:
        print(f"[ERROR] modkit exited with code {result.returncode}.", file=sys.stderr)
        if os.path.exists(temp_output):
            os.remove(temp_output)
        sys.exit(1)

    if not os.path.exists(temp_output) or os.path.getsize(temp_output) == 0:
        print("[ERROR] modkit output file was not created or is empty.", file=sys.stderr)
        if os.path.exists(temp_output):
            os.remove(temp_output)
        sys.exit(1)

    return temp_output


def aggregate_bedgraph(temp_output):
    """
    Aggregate modkit output into bedgraph format for positive and negative strands.

    Streams the modkit TSV line by line to avoid loading 50M+ lines into memory
    at once. Cleans up the temp file after processing.

    For each genomic position:
        coverage = number of reads covering T
        BrdU count = number of reads with mod_code == 'b'

    The output bedgraph will contain:
        chrom, start, end, brdu_fraction, coverage

    Args:
        temp_output: Path to modkit TSV file (will be deleted after processing)

    Returns:
        positive_data: dict for '+' strand
        negative_data: dict for '-' strand
    """
    positive_data = defaultdict(lambda: [0, 0])
    negative_data = defaultdict(lambda: [0, 0])

    with open(temp_output, "r") as infile:
        for line in infile:
            line = line.strip()
            if not line:
                continue

            parts = line.split("\t")

            if len(parts) < 21:
                continue

            if parts[0] == "read_id":
                continue

            ref_position = parts[2]
            chrom = parts[3]
            ref_strand = parts[5]
            mod_code = parts[13]
            canonical_base = parts[17]

            if canonical_base != "T":
                continue

            try:
                ref_position = int(ref_position)
            except ValueError:
                continue

            key = (chrom, ref_position, ref_position + 1)

            if ref_strand == "+":
                target = positive_data
            elif ref_strand == "-":
                target = negative_data
            else:
                continue

            target[key][1] += 1

            if mod_code == "b":
                target[key][0] += 1

    os.remove(temp_output)

    return positive_data, negative_data


def write_bedgraph(data, output_file):
    """
    Write aggregated data to a bedgraph file.

    Each line contains:
        chrom, start, end, brdu_fraction, coverage

    Args:
        data: Dictionary of aggregated values
        output_file: Path to output bedgraph file
    """
    with open(output_file, "w") as out:
        for (chrom, start, end), (brdu_count, coverage) in sorted(
            data.items(), key=lambda x: (x[0][0], x[0][1])
        ):
            if coverage == 0:
                continue

            brdu_fraction = brdu_count / coverage
            out.write(f"{chrom}\t{start}\t{end}\t{brdu_fraction:.8g}\t{coverage}\n")


def main():
    """
    Main workflow:
    1. Resolve the BAM path
    2. Ensure the BAM is sorted
    3. Ensure the BAM is indexed
    4. Extract modkit data
    5. Aggregate the output into bedgraph format
    6. Write positive and negative strand bedgraph files
    """
    args = parse_args()

    sorted_dir, index_dir, bedgraph_dir, bam_dir = get_project_paths()
    os.makedirs(bedgraph_dir, exist_ok=True)

    input_bam = resolve_bam_path(args.bam, bam_dir)

    bam_path = check_and_sort_bam(input_bam, sorted_dir)
    check_and_index_bam(bam_path, index_dir)

    temp_output = extract_brdu_with_modkit(bam_path)
    positive_data, negative_data = aggregate_bedgraph(temp_output)

    bam_base = os.path.splitext(os.path.basename(args.bam))[0]

    if args.output:
        output_prefix = args.output
    else:
        output_prefix = bam_base

    positive_output = os.path.join(
        bedgraph_dir, f"{output_prefix}_positive.bedgraph"
    )
    negative_output = os.path.join(
        bedgraph_dir, f"{output_prefix}_negative.bedgraph"
    )

    write_bedgraph(positive_data, positive_output)
    write_bedgraph(negative_data, negative_output)

    print(f"[INFO] Positive strand bedgraph written to {positive_output}", file=sys.stderr)
    print(f"[INFO] Negative strand bedgraph written to {negative_output}", file=sys.stderr)


if __name__ == "__main__":
    main()
