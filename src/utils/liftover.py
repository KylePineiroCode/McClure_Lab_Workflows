#!/usr/bin/env python3
"""
liftover.py

Run UCSC LiftOver on a BED file to convert genomic coordinates from one
assembly/strain to another.

This utility is intended to be used as a reusable preprocessing step for
downstream workflows such as rain plots, genome browser plots, heatmaps, etc.

Behavior:
- If --mapped and --unmapped are provided, those exact output paths are used.
- If they are not provided, output files are automatically created in the same
  directory as the input BED file using the input BED name + chain file name.

Example automatic naming:
    input BED:   chr1_0_to_50000bases.bed
    chain file:  sacCer1ToSacCer3.over.chain.gz

    mapped:      chr1_0_to_50000bases_sacCer1_to_sacCer3.bed
    unmapped:    chr1_0_to_50000bases_sacCer1_to_sacCer3_unmapped.bed

Example:
    python src/utils/liftover.py \
        data/bed/chr1_0_to_50000bases.bed \
        --chain data/liftover_chains/sacCer1ToSacCer3.over.chain.gz
"""

import argparse
import os
import re
import shutil
import subprocess
import sys


def parse_args():
    """
    Parse command line arguments for LiftOver.
    """
    parser = argparse.ArgumentParser(
        description="Convert genomic coordinates between assemblies using UCSC LiftOver."
    )

    parser.add_argument(
        "input_bed",
        help="Input BED file containing coordinates to convert"
    )

    parser.add_argument(
        "--chain",
        required=True,
        help="Path to UCSC chain file for source -> target conversion"
    )

    parser.add_argument(
        "--mapped",
        required=False,
        default=None,
        help="Optional output BED file for successfully mapped coordinates"
    )

    parser.add_argument(
        "--unmapped",
        required=False,
        default=None,
        help="Optional output BED file for coordinates that failed to map"
    )

    parser.add_argument(
        "--liftover_bin",
        default="liftOver",
        help="Path to UCSC liftOver executable (default: liftOver from PATH)"
    )

    parser.add_argument(
        "--minmatch",
        type=float,
        default=None,
        help="Optional minimum fraction of bases that must remap (passed to liftOver as -minMatch)"
    )

    parser.add_argument(
        "--multiple",
        action="store_true",
        help="Allow multiple output regions for a single input region (-multiple)"
    )

    return parser.parse_args()


def validate_file_exists(path, label):
    """
    Validate that a required input file exists.
    """
    if not os.path.exists(path):
        print(f"[ERROR] {label} not found: {path}", file=sys.stderr)
        sys.exit(1)


def validate_liftover_executable(liftover_bin):
    """
    Validate that the UCSC liftOver executable is available.
    """
    if os.path.sep in liftover_bin or liftover_bin.startswith("."):
        if not os.path.isfile(liftover_bin):
            print(f"[ERROR] liftOver executable not found: {liftover_bin}", file=sys.stderr)
            sys.exit(1)
        if not os.access(liftover_bin, os.X_OK):
            print(f"[ERROR] liftOver is not executable: {liftover_bin}", file=sys.stderr)
            sys.exit(1)
        return

    resolved = shutil.which(liftover_bin)
    if resolved is None:
        print(
            "[ERROR] UCSC liftOver executable was not found in PATH. "
            "Load it as a module, install it, or pass it explicitly with --liftover_bin.",
            file=sys.stderr
        )
        sys.exit(1)


def ensure_parent_dir(path):
    """
    Create the parent directory for an output file if it does not exist.
    """
    parent = os.path.dirname(os.path.abspath(path))
    if parent:
        os.makedirs(parent, exist_ok=True)


def count_nonempty_lines(path):
    """
    Count non-empty, non-comment lines in a file.
    """
    if not os.path.exists(path):
        return 0

    count = 0
    with open(path, "r") as handle:
        for line in handle:
            stripped = line.strip()
            if stripped and not stripped.startswith("#"):
                count += 1
    return count


def get_chain_tag(chain_path):
    """
    Convert a UCSC chain filename into a clean suffix for BED naming.

    Example:
        sacCer1ToSacCer3.over.chain.gz -> sacCer1_to_sacCer3
        sacCer3ToSacCer2.over.chain    -> sacCer3_to_sacCer2
    """
    chain_name = os.path.basename(chain_path)

    if chain_name.endswith(".over.chain.gz"):
        chain_name = chain_name[:-14]
    elif chain_name.endswith(".over.chain"):
        chain_name = chain_name[:-11]
    elif chain_name.endswith(".chain.gz"):
        chain_name = chain_name[:-9]
    elif chain_name.endswith(".chain"):
        chain_name = chain_name[:-6]

    # Convert "...To..." into "..._to_..."
    chain_name = re.sub(r"To", "_to_", chain_name, count=1)

    return chain_name


def build_default_output_paths(input_bed, chain_path):
    """
    Build default mapped and unmapped output paths based on the input BED name
    and the chain filename.

    Example:
        input BED:  data/bed/chr1_0_to_50000bases.bed
        chain:      sacCer1ToSacCer3.over.chain.gz

        mapped:     data/bed/chr1_0_to_50000bases_sacCer1_to_sacCer3.bed
        unmapped:   data/bed/chr1_0_to_50000bases_sacCer1_to_sacCer3_unmapped.bed
    """
    input_dir = os.path.dirname(os.path.abspath(input_bed))
    input_base = os.path.splitext(os.path.basename(input_bed))[0]
    chain_tag = get_chain_tag(chain_path)

    mapped_out = os.path.join(input_dir, f"{input_base}_{chain_tag}.bed")
    unmapped_out = os.path.join(input_dir, f"{input_base}_{chain_tag}_unmapped.bed")

    return mapped_out, unmapped_out


def run_liftover(input_bed, chain_file, mapped_out, unmapped_out, liftover_bin, minmatch=None, multiple=False):
    """
    Run the UCSC LiftOver command.

    Command format:
        liftOver [options] input.bed chain.over.chain mapped.bed unmapped.bed
    """
    cmd = [liftover_bin]

    if minmatch is not None:
        cmd.append(f"-minMatch={minmatch}")

    if multiple:
        cmd.append("-multiple")

    cmd.extend([input_bed, chain_file, mapped_out, unmapped_out])

    print("[INFO] Running UCSC LiftOver...", file=sys.stderr)
    print(f"[INFO] Command: {' '.join(cmd)}", file=sys.stderr)

    try:
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True
        )
        if result.stdout.strip():
            print(result.stdout.strip(), file=sys.stderr)
        if result.stderr.strip():
            print(result.stderr.strip(), file=sys.stderr)

    except FileNotFoundError:
        print(
            "[ERROR] liftOver executable could not be run. "
            "Check your PATH or pass --liftover_bin explicitly.",
            file=sys.stderr
        )
        sys.exit(1)

    except subprocess.CalledProcessError as e:
        print("[ERROR] UCSC LiftOver failed.", file=sys.stderr)
        if e.stdout:
            print(e.stdout, file=sys.stderr)
        if e.stderr:
            print(e.stderr, file=sys.stderr)
        sys.exit(1)


def main():
    args = parse_args()

    validate_file_exists(args.input_bed, "Input BED file")
    validate_file_exists(args.chain, "Chain file")
    validate_liftover_executable(args.liftover_bin)

    # Auto-generate output names if not provided
    if args.mapped is None or args.unmapped is None:
        auto_mapped, auto_unmapped = build_default_output_paths(args.input_bed, args.chain)

        mapped_out = args.mapped if args.mapped is not None else auto_mapped
        unmapped_out = args.unmapped if args.unmapped is not None else auto_unmapped
    else:
        mapped_out = args.mapped
        unmapped_out = args.unmapped

    ensure_parent_dir(mapped_out)
    ensure_parent_dir(unmapped_out)

    input_count = count_nonempty_lines(args.input_bed)
    print(f"[INFO] Input intervals detected: {input_count}", file=sys.stderr)
    print(f"[INFO] Chain file being used: {os.path.basename(args.chain)}", file=sys.stderr)
    print(f"[INFO] Mapped output will be written to: {mapped_out}", file=sys.stderr)
    print(f"[INFO] Unmapped output will be written to: {unmapped_out}", file=sys.stderr)

    run_liftover(
        input_bed=args.input_bed,
        chain_file=args.chain,
        mapped_out=mapped_out,
        unmapped_out=unmapped_out,
        liftover_bin=args.liftover_bin,
        minmatch=args.minmatch,
        multiple=args.multiple
    )

    mapped_count = count_nonempty_lines(mapped_out)
    unmapped_count = count_nonempty_lines(unmapped_out)

    print(f"[INFO] Successfully mapped intervals: {mapped_count}", file=sys.stderr)
    print(f"[INFO] Unmapped intervals: {unmapped_count}", file=sys.stderr)

    if mapped_count == 0:
        print(
            "[WARN] No intervals were successfully lifted over. "
            "This can happen if the source and target assemblies are too different, "
            "the chain file is incorrect for this conversion, or the chromosome names "
            "in the BED file do not match the chain file naming convention.",
            file=sys.stderr
        )
    else:
        print(f"[INFO] LiftOver complete. Mapped output written to: {mapped_out}", file=sys.stderr)

    print(f"[INFO] Unmapped output written to: {unmapped_out}", file=sys.stderr)


if __name__ == "__main__":
    main()
