#!/usr/bin/env python3
"""
Convert a REDIportal hg38 table into the BED6 format expected by PANTRY's
query_editing_level.py.

PANTRY expects a BED-like file with at least 6 columns:
    chrom, start, end, name, score, strand
where strand (column 6) determines which mismatch counts as an edit:
    '+' strand -> G mismatches (A->I reads as A->G)
    '-' strand -> C mismatches (T->I on reverse strand reads as T->C)

REDIportal v3.0 download:
    https://rediportal.cloud.ba.infn.it/ (search / download page)
    Select hg38, download the full table (TSV/CSV).

REDIportal v3.0 TABLE1_hg38 columns:
    Accession, Region, Position, Ref, Ed, Strand, db, type, dbsnp, repeat, ...
The script auto-detects the relevant columns by header name and falls back
to positional defaults if detection fails.

Chromosome naming:
    REDIportal v3 ships chr-prefixed names (chr1, chr2, ...). By default this
    script KEEPS the chr prefix so the output matches chr-prefixed references
    (HPLRv2 GTF, GRCh38 FASTA). If your reference uses bare
    numbers (1, 2, ...), pass --strip-chr.

Usage:
    python3 rediportal_preprocess.py \\
        --input TABLE1_hg38_v3.txt.gz \\
        --output rediportal_hg38.bed

    # For a no-chr reference (e.g. NCBI GRCh38 no-alt analysis set):
    python3 rediportal_preprocess.py \\
        --input TABLE1_hg38_v3.txt.gz \\
        --output rediportal_hg38_nochr.bed \\
        --strip-chr

Notes:
    - Filters to A-to-I events (Ref=A, Ed=G on + strand; or
      Ref=T, Ed=C on - strand, which is the reverse-strand view).
      REDIportal v3 is A-to-I only, so effectively all rows are kept.
    - start is 0-based (Position - 1), end is 1-based (Position).
    - Input may be gzipped (.gz) or plain text.
"""

import argparse
import gzip
import sys
from pathlib import Path


def open_text(path: str):
    """Open a file that may be gzipped (by extension) for text reading."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def detect_columns(header_fields):
    """Map REDIportal column names to indices.

    Returns a dict with keys: chrom, position, strand, ref, base.
    Raises ValueError if required columns cannot be identified.
    """
    # Normalize headers for matching
    norm = [h.strip().lower() for h in header_fields]

    mapping = {}

    # Chromosome: "region" (REDIportal v3), "chrom", "chr", "chromosome"
    for key in ("region", "chrom", "chr", "chromosome", "seqname"):
        if key in norm:
            mapping["chrom"] = norm.index(key)
            break

    # Position: "position" (REDIportal v3), "pos", "coordinate", "genomic_position"
    for key in ("position", "pos", "coordinate", "genomic_position", "start"):
        if key in norm:
            mapping["position"] = norm.index(key)
            break

    # Strand: "strand"
    for key in ("strand",):
        if key in norm:
            mapping["strand"] = norm.index(key)
            break

    # Reference base: "ref" (REDIportal v3), "reference", "refbase", "ref_base"
    for key in ("ref", "reference", "refbase", "ref_base", "genome"):
        if key in norm:
            mapping["ref"] = norm.index(key)
            break

    # Edited base: "ed" (REDIportal v3), "base", "edited", "edit_base", "alt", "alternative"
    for key in ("ed", "base", "edited", "edit_base", "alt", "alternative", "ed_base"):
        if key in norm:
            mapping["base"] = norm.index(key)
            break

    missing = [k for k in ("chrom", "position", "strand", "ref", "base") if k not in mapping]
    if missing:
        raise ValueError(
            f"Could not identify required columns: {missing}.\n"
            f"Detected mapping: {mapping}\n"
            f"Available headers: {header_fields}\n"
            f"If your REDIportal download uses different column names, "
            f"edit the detect_columns() function in this script."
        )

    return mapping


def is_atoi(ref_base, edit_base, strand):
    """Return True if the event is an A-to-I edit.

    A-to-I on the + strand: ref=A, edit=G.
    On the - strand it appears as ref=T, edit=C (the reverse complement).
    REDIportal typically reports the strand and the reference/edited bases
    in genomic orientation, so we accept both representations.
    """
    ref_base = ref_base.upper()
    edit_base = edit_base.upper()
    strand = strand.strip()
    if strand == "+":
        return ref_base == "A" and edit_base == "G"
    elif strand == "-":
        return ref_base == "T" and edit_base == "C"
    # If strand is unknown/ambiguous, accept either A->G or T->C
    return (ref_base == "A" and edit_base == "G") or (ref_base == "T" and edit_base == "C")


def normalize_chrom(chrom, strip_chr=False):
    """Optionally remove 'chr' prefix to match a no-chr reference.

    By default the chr prefix is KEPT (REDIportal native + chr-prefixed refs).
    Pass strip_chr=True for references that use bare chromosome numbers.
    """
    if strip_chr and chrom.lower().startswith("chr"):
        return chrom[3:]
    return chrom


def main():
    parser = argparse.ArgumentParser(
        description="Convert REDIportal hg38 table to PANTRY BED6 edit-sites format"
    )
    parser.add_argument("--input", required=True, help="REDIportal hg38 table (TSV, optionally .gz)")
    parser.add_argument("--output", required=True, help="Output BED6 file (no .gz; bgzip separately if desired)")
    parser.add_argument("--name-col", default=None,
                        help="Optional column name to use for the BED 'name' field (4th col). "
                             "Defaults to 'rediportal_<chrom>_<position>'.")
    parser.add_argument("--strip-chr", action="store_true",
                        help="Remove 'chr' prefix from chromosome names (use this ONLY if your "
                             "reference uses bare numbers, e.g. NCBI GRCh38 no-alt analysis set). "
                             "Default: keep chr prefix to match chr-prefixed references.")
    parser.add_argument("--verbose", action="store_true", help="Print progress")
    args = parser.parse_args()

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    n_total = 0
    n_kept = 0
    n_skipped_non_ag = 0

    with open_text(args.input) as infile, open(output_path, "w") as outfile:
        header = None
        colmap = None
        for line in infile:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")

            if header is None:
                # First non-comment, non-empty line is the header
                header = fields
                try:
                    colmap = detect_columns(header)
                except ValueError as e:
                    # If header detection fails, try positional defaults for
                    # the common REDIportal v3.0 layout:
                    # Accession, Region, Position, Ref, Ed, Strand, ...
                    if args.verbose:
                        print(f"Header detection failed ({e}); trying positional defaults", file=sys.stderr)
                    if len(fields) >= 6:
                        colmap = {"chrom": 1, "position": 2, "ref": 3, "base": 4, "strand": 5}
                    else:
                        raise
                if args.verbose:
                    print(f"Column mapping: {colmap}", file=sys.stderr)
                    print(f"Strip chr prefix: {args.strip_chr}", file=sys.stderr)
                continue

            n_total += 1
            try:
                chrom = fields[colmap["chrom"]]
                position = int(fields[colmap["position"]])
                strand = fields[colmap["strand"]].strip()
                ref_base = fields[colmap["ref"]]
                edit_base = fields[colmap["base"]]
            except (IndexError, ValueError) as e:
                if args.verbose and n_total <= 5:
                    print(f"Skipping malformed line {n_total}: {e}", file=sys.stderr)
                continue

            # Filter to A-to-I events
            if not is_atoi(ref_base, edit_base, strand):
                n_skipped_non_ag += 1
                continue

            chrom = normalize_chrom(chrom, strip_chr=args.strip_chr)
            start = position - 1  # BED is 0-based
            end = position         # 1-based end
            name = f"rediportal_{chrom}_{position}"
            score = "0"

            outfile.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")
            n_kept += 1

            if args.verbose and n_kept > 0 and n_kept % 1000000 == 0:
                print(f"Kept {n_kept} A-to-I sites ({n_total} processed)", file=sys.stderr)

    print(
        f"Done. Processed {n_total} sites. "
        f"Kept {n_kept} A-to-I sites. "
        f"Skipped {n_skipped_non_ag} non-A-to-I sites. "
        f"Output: {output_path}",
        file=sys.stderr,
    )

    if n_kept == 0:
        print(
            "WARNING: No sites retained. Check that the input is an hg38 REDIportal "
            "table and that column detection worked (use --verbose).",
            file=sys.stderr,
        )
        sys.exit(1)


if __name__ == "__main__":
    main()
