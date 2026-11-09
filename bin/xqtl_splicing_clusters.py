#!/usr/bin/env python3
"""xqtl_splicing_clusters.py -- frozen + site-local leafCutter clusters for the splicing modality.

Splicing at every site is quantified against the frozen cluster set from the restricted file
(leafCutter clusters pooled across the MDACC discovery cohorts), so the same phenotype IDs exist at every site and can be
meta-analysed. Site-local leafCutter clusters that share no intron with the frozen set are kept
as novel clusters (cluster token nclu_N_strand) and are flagged in the export.

  frozen-refined  frozen_clusters.tsv.gz -> leafCutter "refined" cluster file (input to
                  leafcutter_cluster_regtools_py3.py -c) + line-number -> cluster_id map
  merge           frozen-quantified numerators + site-local numerators -> one numerator table
                  (leafCutter perind_numers format) for xqtl_assemble_bed.py splicing

Intron coordinates are leafCutter intron coordinates (the start/end in leafCutter intron IDs
chrom:start:end:clu_N_strand), as stored in splicing/frozen_clusters.tsv.gz.
"""
import argparse
import gzip
import re
import sys

import pandas as pd

TOKEN = re.compile(r"^(?P<chrom>.+):(?P<start>\d+):(?P<end>\d+):(?P<clu>n?clu_\d+_[+\-?])$")


def frozen_refined(a):
    fz = pd.read_csv(a.frozen, sep="\t", dtype={"chrom": str, "cluster_id": str})
    need = {"chrom", "start", "end", "strand", "cluster_id"}
    if not need <= set(fz.columns):
        sys.exit(f"frozen clusters table lacks columns {sorted(need - set(fz.columns))}")
    bad = ~fz["cluster_id"].str.match(r"^clu_\d+_[+\-]$")
    if bad.any():
        sys.exit(f"frozen cluster_id values must look like clu_<N>_<strand>; e.g. {fz.loc[bad, 'cluster_id'].iloc[0]}")
    if a.chroms:
        fz = fz[fz["chrom"].isin(a.chroms.split(","))]
    fz = fz.drop_duplicates(["chrom", "start", "end", "strand", "cluster_id"])
    with open(a.out, "w") as out, open(a.map, "w") as mp:
        mp.write("line\tcluster_id\n")
        n = 0
        for (cid, chrom, strand), g in fz.groupby(["cluster_id", "chrom", "strand"], sort=False):
            n += 1
            g = g.sort_values(["start", "end"])
            out.write(f"{chrom}:{strand} " + " ".join(f"{s}:{e}:0" for s, e in zip(g["start"], g["end"])) + "\n")
            mp.write(f"{n}\t{cid}\n")
    print(f"frozen clusters: {n} written to {a.out}", file=sys.stderr)


def read_numers(path):
    with gzip.open(path, "rt") as fh:
        samples = fh.readline().split()
    try:
        df = pd.read_csv(path, sep=" ", skiprows=1, header=None, index_col=0)
    except pd.errors.EmptyDataError:  # no clusters passed leafCutter's filters
        return pd.DataFrame(columns=samples, index=pd.Index([], dtype=object), dtype="int64")
    df.columns = samples
    df.index.name = None
    return df


def merge(a):
    fz = read_numers(a.frozen_numers)
    cmap = pd.read_csv(a.frozen_map, sep="\t", dtype={"cluster_id": str}).set_index("line")["cluster_id"]
    p = fz.index.to_series().astype(object).str.extract(TOKEN)
    line = p["clu"].str.extract(r"^clu_(\d+)_")[0].astype(int)
    fz.index = p["chrom"] + ":" + p["start"] + ":" + p["end"] + ":" + line.map(cmap).values
    fz_introns = set(zip(p["chrom"], p["start"].astype(int), p["end"].astype(int), p["clu"].str[-1]))

    loc = read_numers(a.local_numers)
    q = loc.index.to_series().astype(object).str.extract(TOKEN)
    q["key"] = list(zip(q["chrom"], q["start"].astype(int), q["end"].astype(int), q["clu"].str[-1]))
    q["in_frozen"] = q["key"].isin(fz_introns)
    matched = q.groupby("clu")["in_frozen"].transform("any")
    nov = loc[~matched.values].copy()
    qn = q[~matched.values]
    nov.index = qn["chrom"] + ":" + qn["start"] + ":" + qn["end"] + ":n" + qn["clu"]
    if list(nov.columns) != list(fz.columns):
        nov = nov[fz.columns]
    out = pd.concat([fz, nov])
    with gzip.open(a.out, "wt") as fh:
        fh.write(" ".join(out.columns) + "\n")
        out.to_csv(fh, sep=" ", header=False)
    fz_tot = fz.groupby(line.values).sum().sum(axis=1)
    summ = pd.DataFrame([
        {"set": "frozen", "clusters": int(p["clu"].nunique()), "introns": len(fz),
         "clusters_with_reads": int((fz_tot > 0).sum())},
        {"set": "local_matched_to_frozen", "clusters": int(q.loc[matched.values, "clu"].nunique()),
         "introns": int(matched.sum()), "clusters_with_reads": int(q.loc[matched.values, "clu"].nunique())},
        {"set": "novel", "clusters": int(qn["clu"].nunique()), "introns": len(nov),
         "clusters_with_reads": int(qn["clu"].nunique())},
    ])
    summ.to_csv(a.summary, sep="\t", index=False)
    print(summ.to_string(index=False), file=sys.stderr)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("frozen-refined")
    f.add_argument("--frozen", required=True)
    f.add_argument("--chroms", default=None, help="comma-separated chromosomes to keep (default all)")
    f.add_argument("--out", required=True)
    f.add_argument("--map", required=True)
    m = sub.add_parser("merge")
    m.add_argument("--frozen-numers", required=True)
    m.add_argument("--frozen-map", required=True)
    m.add_argument("--local-numers", required=True)
    m.add_argument("--out", required=True)
    m.add_argument("--summary", required=True)
    a = ap.parse_args(argv)
    {"frozen-refined": frozen_refined, "merge": merge}[a.cmd](a)


if __name__ == "__main__":
    main()
