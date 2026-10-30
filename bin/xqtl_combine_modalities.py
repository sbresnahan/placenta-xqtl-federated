#!/usr/bin/env python3
"""xqtl_combine_modalities.py -- PANTRY-style combined cross-modality BED (atlas 30_combine_modalities.py).

  xqtl_combine_modalities.py --out EUR_combined expression=EUR_expression.bed.gz splicing=EUR_splicing.bed.gz ...

phenotype_id becomes {modality}__{original_id}; sample columns are intersected across modalities;
every phenotype of a gene (original ID up to the first '__') across modalities forms one group.
Outputs OUT.bed.gz, OUT.phenotype_groups.txt (no header) and OUT.phenotype_modalities.tsv.
"""
import argparse
import sys

import pandas as pd


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("beds", nargs="+", help="modality=path")
    a = ap.parse_args(argv)
    frames, samples = [], None
    for spec in a.beds:
        mod, path = spec.split("=", 1)
        df = pd.read_csv(path, sep="\t", dtype={"#chr": str})
        df.columns = ["#chr", "start", "end", "phenotype_id"] + list(df.columns[4:])
        if df.phenotype_id.duplicated().any():
            sys.exit(f"xqtl_combine_modalities.py: duplicate phenotype IDs in {path}")
        df.insert(4, "modality", mod)
        frames.append(df)
        cols = list(df.columns[5:])
        samples = cols if samples is None else [s for s in samples if s in set(cols)]
    if not samples or len(samples) < 10:
        sys.exit("xqtl_combine_modalities.py: fewer than 10 samples shared by all modalities")
    rows = []
    for df in frames:
        d = df[["#chr", "start", "end", "phenotype_id", "modality"] + samples].copy()
        d["original_id"] = d["phenotype_id"]
        d["gene_id"] = d["original_id"].str.split("__", n=1).str[0]
        d["phenotype_id"] = d["modality"] + "__" + d["original_id"]
        rows.append(d)
    comb = pd.concat(rows, ignore_index=True)
    # grouped mapping needs one cis window per group: use the gene's first (expression) coordinates
    first = comb.groupby("gene_id")[["#chr", "start", "end"]].transform("first")
    comb[["#chr", "start", "end"]] = first
    comb = comb.sort_values(["#chr", "start", "end", "phenotype_id"])
    comb[["#chr", "start", "end", "phenotype_id"] + samples].to_csv(a.out + ".bed.gz", sep="\t", index=False)
    comb[["phenotype_id", "gene_id"]].to_csv(a.out + ".phenotype_groups.txt", sep="\t", index=False, header=False)
    comb[["phenotype_id", "modality", "gene_id", "original_id"]].to_csv(
        a.out + ".phenotype_modalities.tsv", sep="\t", index=False)
    print(f"{a.out}: {len(comb)} phenotypes from {len(frames)} modalities, {comb.gene_id.nunique()} genes, "
          f"{len(samples)} samples")


if __name__ == "__main__":
    main()
