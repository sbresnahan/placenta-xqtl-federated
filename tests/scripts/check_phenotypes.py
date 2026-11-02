#!/usr/bin/env python3
"""check_phenotypes.py -- checks on a phenotyping run (tests/scripts/test_phenotyping.sh).

  check_phenotypes.py OUTDIR --gtf GTF --metadata metadata.tsv --restricted RESTRICTED.tar.gz --modalities a,b,...

Prints one PASS/FAIL line per check and exits 1 if any check fails.
"""
import argparse
import gzip
import io
import sys
import tarfile

import numpy as np
import pandas as pd

GROUPED = {"isoform_expression", "isoform_usage", "splicing", "intron_retention", "alt_tss", "alt_polya"}
fails = 0


def report(ok, msg):
    global fails
    print(("PASS  " if ok else "FAIL  ") + msg)
    fails += 0 if ok else 1


def gene_tss(gtf):
    rows = []
    with gzip.open(gtf, "rt") if gtf.endswith(".gz") else open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            gid = f[8].split('gene_id "', 1)[1].split('"', 1)[0]
            tss = int(f[3]) - 1 if f[6] == "+" else int(f[4]) - 1   # 0-based BED start of the TSS base
            rows.append((gid, f[0], tss))
    return pd.DataFrame(rows, columns=["gene_id", "chrom", "tss"]).set_index("gene_id")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir")
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--restricted", required=True)
    ap.add_argument("--modalities", required=True)
    a = ap.parse_args()
    mods = a.modalities.split(",")
    meta = pd.read_csv(a.metadata, sep="\t", dtype=str)
    rna_ids, geno_ids = list(meta["rnaseq_id"]), sorted(meta["genotype_id"])
    tss = gene_tss(a.gtf)
    ph, nm = f"{a.outdir}/phenotypes", f"{a.outdir}/normalized"

    for m in mods:
        try:
            bed = pd.read_csv(f"{ph}/{m}.bed.gz", sep="\t")
        except FileNotFoundError:
            report(False, f"{m}: phenotypes/{m}.bed.gz missing")
            continue
        samples = list(bed.columns[4:])
        ok = (len(bed) > 0 and sorted(samples) == sorted(rna_ids) and (bed["end"] - bed["start"] == 1).all()
              and bed["phenotype_id"].is_unique)
        report(ok, f"{m}: {len(bed)} features x {len(samples)} samples, 1-bp TSS windows, unique IDs")
        gid = bed["phenotype_id"].str.replace(r"__.*$", "", regex=True)
        known = gid.isin(tss.index)
        exp = tss.reindex(gid[known])["tss"].to_numpy()
        report(known.mean() > 0.99 and (bed.loc[known, "start"].to_numpy() == exp).all(),
               f"{m}: positions are the strand-aware gene TSS from the annotation ({known.sum()} genes)")
        try:
            norm = pd.read_csv(f"{nm}/ALL_{m}.bed.gz", sep="\t")
        except FileNotFoundError:
            report(False, f"{m}: normalized/ALL_{m}.bed.gz missing")
            continue
        vals = norm.iloc[:, 4:].to_numpy(dtype=float)
        report(sorted(norm.columns[4:]) == geno_ids and len(norm) > 0 and np.isfinite(vals).all()
               and set(norm["phenotype_id"]) <= set(bed["phenotype_id"]),
               f"{m}: normalized {len(norm)} features, genotype IDs as columns, all values finite")
        if m in GROUPED:
            g = pd.read_csv(f"{nm}/ALL_{m}.phenotype_groups.txt", sep="\t", header=None, names=["pid", "grp"])
            report(set(g["pid"]) == set(norm["phenotype_id"])
                   and (g["pid"].str.replace(r"__.*$", "", regex=True) == g["grp"]).all(),
                   f"{m}: phenotype groups cover the normalized features and map each to its gene")

    if "splicing" in mods:
        with tarfile.open(a.restricted) as t:
            fz = pd.read_csv(io.BytesIO(t.extractfile("splicing/frozen_clusters.tsv.gz").read()),
                             sep="\t", compression="gzip")
        bed = pd.read_csv(f"{ph}/splicing.bed.gz", sep="\t")
        clu = bed["phenotype_id"].str.extract(r"_(n?clu_\d+_[+\-.])$")[0]
        novel = clu.str.startswith("n")
        frozen_ok = set(clu[~novel]) <= set(fz["cluster_id"])
        report(novel.sum() > 0 and (~novel).sum() > 0 and frozen_ok,
               f"splicing: {(~novel).sum()} junctions in frozen clusters (IDs from the restricted file), "
               f"{novel.sum()} in novel clusters")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
