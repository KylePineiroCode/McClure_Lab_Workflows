#!/usr/bin/env python3
import argparse
import os
import sys


# -----------------------------
# Chromosome ID mappings
# -----------------------------
genbank_to_chr = {
    "CM007964.1": "1",  "CM007965.1": "2",  "CM007966.1": "3",  "CM007967.1": "4",
    "CM007968.1": "5",  "CM007969.1": "6",  "CM007970.1": "7",  "CM007971.1": "8",
    "CM007972.1": "9",  "CM007973.1": "10", "CM007974.1": "11", "CM007975.1": "12",
    "CM007976.1": "13", "CM007977.1": "14", "CM007978.1": "15", "CM007979.1": "16",
    "CM007980.1": "p2-micron", "CM007981.1": "MT"
}

ncbi_refseq_to_chr = {
    "NC_001133.9": "1",   "NC_001134.8": "2",   "NC_001135.5": "3",   "NC_001136.10": "4",
    "NC_001137.3": "5",   "NC_001138.5": "6",   "NC_001139.9": "7",   "NC_001140.6": "8",
    "NC_001141.2": "9",   "NC_001142.9": "10",  "NC_001143.9": "11",  "NC_001144.5": "12",
    "NC_001145.3": "13",  "NC_001146.8": "14",  "NC_001147.6": "15",  "NC_001148.4": "16",
    "NC_001224.1": "MT"
}


def normalize_chrom_id(chrom_id: str):
    """
    Convert chromosome IDs from GenBank or RefSeq to a common internal label.
    Example:
      CM007964.1 -> 1
      NC_001133.9 -> 1
    If unknown, return as-is.
    """
    if chrom_id in genbank_to_chr:
        return genbank_to_chr[chrom_id]
    if chrom_id in ncbi_refseq_to_chr:
        return ncbi_refseq_to_chr[chrom_id]
    return chrom_id


def parse_args():
    parser = argparse.ArgumentParser(
        description="Parse a GFF/GFF3 file and extract genomic features in a region to a BED-like file."
    )
    parser.add_argument("--gff", required=True, help="Path to input GFF/GFF3 file")
    parser.add_argument("--chrom", required=True, help="Chromosome/seqid to extract")
    parser.add_argument("--start", required=True, type=int, help="Region start (0-based)")
    parser.add_argument("--end", required=True, type=int, help="Region end (0-based, exclusive)")
    parser.add_argument("--output", required=True, help="Output BED file")
    return parser.parse_args()


def parse_attributes(attr_str):
    attrs = {}
    if not attr_str:
        return attrs

    for field in attr_str.strip().split(";"):
        field = field.strip()
        if not field:
            continue

        if "=" in field:
            key, value = field.split("=", 1)
            attrs[key.strip()] = value.strip()
        elif " " in field:
            key, value = field.split(" ", 1)
            attrs[key.strip()] = value.strip().strip('"')

    return attrs


def choose_feature_name(feature_type, attrs):
    for key in ("gene", "Name", "locus_tag", "ID", "product", "Parent"):
        if key in attrs and attrs[key]:
            return attrs[key]
    return feature_type


def normalize_feature_type(feature_type, attrs):
    feature_type = feature_type.strip()

    if feature_type == "gene":
        return "gene"
    if feature_type in {"mRNA", "transcript", "CDS", "exon"}:
        return feature_type
    if feature_type in {"tRNA", "rRNA", "ncRNA", "snoRNA", "snRNA"}:
        return feature_type
    if feature_type in {"centromere", "telomere", "repeat_region", "origin_of_replication"}:
        return feature_type

    # NCBI sometimes uses long_terminal_repeat instead of repeat_region
    if feature_type == "long_terminal_repeat":
        return "repeat_region"

    name = attrs.get("Name", "").lower()
    note = attrs.get("Note", "").lower()
    product = attrs.get("product", "").lower()
    gbkey = attrs.get("gbkey", "").lower()

    if "ars" in name or "autonomously replicating sequence" in note or gbkey == "rep_origin":
        return "origin_of_replication"

    return feature_type


def keep_feature(feature_type):
    wanted = {
        "gene",
        "mRNA",
        "transcript",
        "CDS",
        "exon",
        "tRNA",
        "rRNA",
        "ncRNA",
        "snoRNA",
        "snRNA",
        "centromere",
        "telomere",
        "repeat_region",
        "origin_of_replication",
    }
    return feature_type in wanted


def overlaps(feature_start, feature_end, region_start, region_end):
    return feature_end > region_start and feature_start < region_end


def main():
    args = parse_args()

    if not os.path.exists(args.gff):
        print(f"[ERROR] GFF file not found: {args.gff}", file=sys.stderr)
        sys.exit(1)

    if args.start >= args.end:
        print("[ERROR] start must be less than end.", file=sys.stderr)
        sys.exit(1)

    requested_chrom = normalize_chrom_id(args.chrom)
    rows_written = 0

    with open(args.gff, "r") as gff_in, open(args.output, "w") as out:
        for line in gff_in:
            if not line.strip() or line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue

            seqid, source, feature_type, start_1based, end_1based, score, strand, phase, attributes = parts

            seqid_normalized = normalize_chrom_id(seqid)
            if seqid_normalized != requested_chrom:
                continue

            try:
                gff_start = int(start_1based) - 1
                gff_end = int(end_1based)
            except ValueError:
                continue

            if not overlaps(gff_start, gff_end, args.start, args.end):
                continue

            attrs = parse_attributes(attributes)
            normalized_type = normalize_feature_type(feature_type, attrs)

            if not keep_feature(normalized_type):
                continue

            clipped_start = max(gff_start, args.start)
            clipped_end = min(gff_end, args.end)

            if clipped_start >= clipped_end:
                continue

            feature_name = choose_feature_name(normalized_type, attrs)
            strand = strand if strand in {"+", "-", "."} else "."

            # Write output using the USER'S chromosome ID so it matches the rest of the workflow
            out.write(
                f"{args.chrom}\t{clipped_start}\t{clipped_end}\t"
                f"{normalized_type}\t{feature_name}\t{strand}\n"
            )
            rows_written += 1

    print(f"[INFO] Wrote {rows_written} feature rows to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
