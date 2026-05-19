#!/usr/bin/env python3
"""
Lift over a BrdU BED while preserving the original non-coordinate columns.

Input format is expected to be at least:
    chrom  start  end  read_id  mod_base  mod_prob

The UCSC liftOver binary treats BED6 specially, which can corrupt the BrdU
columns when they do not match BED6 semantics. To avoid that, this utility:

1. Writes a temporary BED4 using a synthetic row ID
2. Runs liftOver on that BED4
3. Joins lifted coordinates back to the original BrdU columns

Mapped output keeps lifted chrom/start/end plus the original remaining columns.
Unmapped output keeps the original row content for rows that fail to map.
"""

import argparse
import csv
import os
import shutil
import subprocess
import sys
import tempfile


def parse_args():
    """
    Parse command-line arguments for BrdU BED liftOver.

    The script needs an input BrdU BED file, a UCSC chain file, mapped and
    unmapped output paths, and optionally a custom path to the liftOver binary.

    Returns
    -------
    argparse.Namespace
        Parsed arguments containing input BED, chain file, mapped output,
        unmapped output, and liftOver executable path.
    """
    parser = argparse.ArgumentParser(
        description="Lift over a BrdU BED while preserving original BrdU columns."
    )
    parser.add_argument("input_bed", help="Input BrdU BED file")
    parser.add_argument("--chain", required=True, help="Path to UCSC chain file")
    parser.add_argument("--mapped", required=True, help="Mapped output BED path")
    parser.add_argument("--unmapped", required=True, help="Unmapped output BED path")
    parser.add_argument(
        "--liftover_bin",
        default="liftOver",
        help="Path to UCSC liftOver executable"
    )
    return parser.parse_args()


def validate_path_exists(path, label):
    """
    Confirm that a required input file exists.

    Parameters
    ----------
    path : str
        File path to check.

    label : str
        Human-readable label used in the error message.

    Returns
    -------
    None
        Exits the script if the path does not exist.
    """
    if not os.path.exists(path):
        print(f"[ERROR] {label} not found: {path}", file=sys.stderr)
        sys.exit(1)


def validate_liftover_executable(path):
    """
    Validate and resolve the UCSC liftOver executable.

    The executable can be provided either as:
    1. A direct path, such as /path/to/liftOver
    2. A command available in PATH, such as liftOver

    Direct paths are checked for file existence and execute permissions.
    Commands without path separators are resolved using shutil.which().

    Parameters
    ----------
    path : str
        Path or command name for the liftOver executable.

    Returns
    -------
    str
        Resolved path or command to use for running liftOver.
    """
    # If the value looks like a path, validate it directly instead of searching
    # PATH. This prevents accidentally using a different liftOver binary.
    if os.path.sep in path or path.startswith("."):
        if not os.path.isfile(path) or not os.access(path, os.X_OK):
            print(f"[ERROR] liftOver executable not usable: {path}", file=sys.stderr)
            sys.exit(1)
        return path

    # Otherwise, treat it as a command name and look for it in PATH.
    resolved = shutil.which(path)
    if resolved is None:
        print(f"[ERROR] liftOver executable not found in PATH: {path}", file=sys.stderr)
        sys.exit(1)

    return resolved


def read_input_rows(path):
    """
    Read the original BrdU BED rows.

    The expected input is a tab-delimited BED-like file with at least six
    columns:

        chrom, start, end, read_id, mod_base, mod_prob

    Extra columns are preserved because the reconstruction step keeps all
    original columns after the first three coordinates.

    Parameters
    ----------
    path : str
        Path to the input BrdU BED file.

    Returns
    -------
    list[list[str]]
        List of original BED rows, where each row is a list of string columns.
    """
    rows = []

    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")

        for idx, row in enumerate(reader):
            if not row:
                continue

            if len(row) < 6:
                print(
                    f"[ERROR] Expected at least 6 BED columns, found {len(row)} on line {idx + 1}",
                    file=sys.stderr,
                )
                sys.exit(1)

            rows.append(row)

    return rows


def write_temp_bed4(rows, path):
    """
    Write a temporary BED4 file for UCSC liftOver.

    UCSC liftOver can treat BED6 columns as special fields such as name, score,
    and strand. BrdU BED files use columns 4-6 differently, so passing them
    directly to liftOver can corrupt the BrdU data.

    To avoid that, this function writes only:

        chrom, start, end, synthetic_row_id

    The synthetic row ID is later used to reconnect lifted coordinates to the
    original BrdU columns.

    Parameters
    ----------
    rows : list[list[str]]
        Original BrdU BED rows.

    path : str
        Path where the temporary BED4 file should be written.
    """
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")

        for idx, row in enumerate(rows):
            writer.writerow([row[0], row[1], row[2], f"row_{idx}"])


def run_liftover(liftover_bin, temp_input, chain, temp_mapped, temp_unmapped):
    """
    Run UCSC liftOver on the temporary BED4 input.

    Parameters
    ----------
    liftover_bin : str
        Path or command name for the liftOver executable.

    temp_input : str
        Temporary BED4 input file.

    chain : str
        UCSC chain file used for coordinate conversion.

    temp_mapped : str
        Temporary file where successfully lifted intervals are written.

    temp_unmapped : str
        Temporary file where failed intervals are written.

    Returns
    -------
    None
        Exits the script if liftOver fails.
    """
    cmd = [liftover_bin, temp_input, chain, temp_mapped, temp_unmapped]

    print("[INFO] Running UCSC LiftOver for BrdU BED coordinate remapping...", file=sys.stderr)
    print(f"[INFO] Command: {' '.join(cmd)}", file=sys.stderr)

    try:
        # capture_output=True lets us show liftOver's stdout/stderr if it fails,
        # which makes chain-file or coordinate errors easier to debug.
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as exc:
        print("[ERROR] UCSC LiftOver failed.", file=sys.stderr)

        if exc.stdout:
            print(exc.stdout, file=sys.stderr)
        if exc.stderr:
            print(exc.stderr, file=sys.stderr)

        sys.exit(1)

    if result.stdout.strip():
        print(result.stdout.strip(), file=sys.stderr)
    if result.stderr.strip():
        print(result.stderr.strip(), file=sys.stderr)


def read_liftover_rows(path):
    """
    Read rows produced by UCSC liftOver.

    Parameters
    ----------
    path : str
        Path to either the temporary mapped or unmapped liftOver output.

    Returns
    -------
    list[list[str]]
        Rows read from the file. Returns an empty list if the file does not
        exist.
    """
    rows = []

    if not os.path.exists(path):
        return rows

    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")

        for row in reader:
            if row:
                rows.append(row)

    return rows


def write_reconstructed_outputs(original_rows, mapped_rows, unmapped_rows, mapped_out, unmapped_out):
    """
    Rebuild mapped and unmapped BrdU BED outputs after liftOver.

    The temporary BED4 sent to liftOver contains synthetic row IDs, such as
    row_0, row_1, row_2. Those IDs are used here to look up the original BrdU
    row and restore its non-coordinate columns.

    The mapped output receives:
        lifted_chrom, lifted_start, lifted_end, original_read_id, original_base, original_prob, ...

    The unmapped output receives the original full row because those intervals
    failed to lift and therefore do not have new coordinates.

    Parameters
    ----------
    original_rows : list[list[str]]
        Original BrdU BED rows.

    mapped_rows : list[list[str]]
        liftOver mapped rows from the temporary BED4 file.

    unmapped_rows : list[list[str]]
        liftOver unmapped rows from the temporary BED4 file.

    mapped_out : str
        Final mapped output BED path.

    unmapped_out : str
        Final unmapped output BED path.

    Returns
    -------
    None
    """
    # Map each synthetic row ID back to the original row so we can restore the
    # BrdU-specific columns after liftOver changes only chrom/start/end.
    row_lookup = {f"row_{idx}": row for idx, row in enumerate(original_rows)}

    mapped_count = 0

    with open(mapped_out, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")

        for row in mapped_rows:
            if len(row) < 4:
                continue

            original = row_lookup.get(row[3])
            if original is None:
                continue

            # Use lifted coordinates from liftOver, then append the original
            # BrdU columns starting at column 4.
            writer.writerow([row[0], row[1], row[2], *original[3:]])
            mapped_count += 1

    unmapped_count = 0

    with open(unmapped_out, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")

        for row in unmapped_rows:
            if len(row) < 4:
                continue

            original = row_lookup.get(row[3])
            if original is None:
                continue

            # Unmapped intervals keep their original coordinates because they
            # did not successfully convert to the target coordinate system.
            writer.writerow(original)
            unmapped_count += 1

    print(f"[INFO] Successfully mapped intervals: {mapped_count}", file=sys.stderr)
    print(f"[INFO] Unmapped intervals: {unmapped_count}", file=sys.stderr)
    print(f"[INFO] LiftOver complete. Mapped output written to: {mapped_out}", file=sys.stderr)
    print(f"[INFO] Unmapped output written to: {unmapped_out}", file=sys.stderr)


def main():
    """
    Run the full BrdU BED liftOver workflow.

    This function validates inputs, reads the original BrdU BED rows, creates a
    temporary BED4 file, runs UCSC liftOver, then reconstructs final mapped and
    unmapped BED files while preserving the original BrdU data columns.
    """
    args = parse_args()

    validate_path_exists(args.input_bed, "Input BED file")
    validate_path_exists(args.chain, "Chain file")
    liftover_bin = validate_liftover_executable(args.liftover_bin)

    original_rows = read_input_rows(args.input_bed)

    print(f"[INFO] Input intervals detected: {len(original_rows)}", file=sys.stderr)
    print(f"[INFO] Chain file being used: {os.path.basename(args.chain)}", file=sys.stderr)
    print(f"[INFO] Mapped output will be written to: {args.mapped}", file=sys.stderr)
    print(f"[INFO] Unmapped output will be written to: {args.unmapped}", file=sys.stderr)

    # TemporaryDirectory automatically removes intermediate BED files when this
    # block finishes, keeping the workflow directory clean.
    with tempfile.TemporaryDirectory(prefix="liftover_brdu_") as temp_dir:
        temp_input = os.path.join(temp_dir, "input.bed")
        temp_mapped = os.path.join(temp_dir, "mapped.bed")
        temp_unmapped = os.path.join(temp_dir, "unmapped.bed")

        write_temp_bed4(original_rows, temp_input)
        run_liftover(liftover_bin, temp_input, args.chain, temp_mapped, temp_unmapped)

        mapped_rows = read_liftover_rows(temp_mapped)
        unmapped_rows = read_liftover_rows(temp_unmapped)

        write_reconstructed_outputs(
            original_rows=original_rows,
            mapped_rows=mapped_rows,
            unmapped_rows=unmapped_rows,
            mapped_out=args.mapped,
            unmapped_out=args.unmapped,
        )


if __name__ == "__main__":
    main()