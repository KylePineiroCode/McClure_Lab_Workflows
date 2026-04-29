#!/usr/bin/env python3
"""
Convert BED chromosome names from GenBank accession format to UCSC yeast format.

Example:
    CM007964.1 -> chrI
    CM007965.1 -> chrII

This is only intended as a preprocessing step before running UCSC liftOver.
It preserves all BED columns and only changes column 1 (chromosome name).
"""

import argparse
import sys


GENBANK_TO_UCSC = {
    "CM007964.1": "chrI",
    "CM007965.1": "chrII",
    "CM007966.1": "chrIII",
    "CM007967.1": "chrIV",
    "CM007968.1": "chrV",
    "CM007969.1": "chrVI",
    "CM007970.1": "chrVII",
    "CM007971.1": "chrVIII",
    "CM007972.1": "chrIX",
    "CM007973.1": "chrX",
    "CM007974.1": "chrXI",
    "CM007975.1": "chrXII",
    "CM007976.1": "chrXIII",
    "CM007977.1": "chrXIV",
    "CM007978.1": "chrXV",
    "CM007979.1": "chrXVI",
    "CM007980.1": "chr2micron",
    "CM007981.1": "chrM",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert BED chromosome names from GenBank accession format to UCSC yeast format."
    )
    parser.add_argument("input_bed", help="Input BED file with GenBank chromosome names")
    parser.add_argument("output_bed", help="Output BED file with UCSC chromosome names")
    return parser.parse_args()


def main():
    args = parse_args()

    converted = 0
    unchanged = 0
    unknown_chroms = set()

    with open(args.input_bed, "r") as infile, open(args.output_bed, "w") as outfile:
        for line in infile:
            stripped = line.rstrip("\n")

            if not stripped:
                outfile.write("\n")
                continue

            parts = stripped.split("\t")
            chrom = parts[0]

            if chrom in GENBANK_TO_UCSC:
                parts[0] = GENBANK_TO_UCSC[chrom]
                converted += 1
            else:
                unchanged += 1
                unknown_chroms.add(chrom)

            outfile.write("\t".join(parts) + "\n")

    print(f"[INFO] BED chromosome conversion complete.", file=sys.stderr)
    print(f"[INFO] Converted rows: {converted}", file=sys.stderr)
    print(f"[INFO] Unchanged rows: {unchanged}", file=sys.stderr)

    if unknown_chroms:
        unknown_list = ", ".join(sorted(unknown_chroms))
        print(
            f"[WARN] Some chromosome names were not found in the GenBank->UCSC mapping "
            f"and were left unchanged: {unknown_list}",
            file=sys.stderr
        )


if __name__ == "__main__":
    main()
