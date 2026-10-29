#!/usr/bin/env python3
"""Build the two tables of the restricted file from atlas outputs (MDACC only).

Called by scripts/make_restricted_bundle.sh.

frozen_clusters.tsv.gz  chrom, start, end, strand, cluster_id, atlas_clusters
    The atlas pools leafCutter clusters across cohorts within each ancestry (step 17:
    clusters sharing a junction are merged). Here the per-ancestry pooled clusters are
    merged once more across ancestries the same way (connected components over shared
    junctions), giving one frozen cluster set with IDs clu_1, clu_2, ... ordered by
    position. atlas_clusters lists the per-ancestry clusters it came from (EAS:clu_12;...).

mdacc_signals.tsv.gz  phenotype_id, variant_id, modality, ancestry, signal_rank,
                      atlas_phenotype_id, scan, feature_key
    Lead signals with q <= QVALUE_MAX from {ANC}_{MOD}_cisqtl.parquet plus conditionally
    independent signals from {ANC}_{MOD}_cisqtl_independent.parquet, for the nine atlas
    modalities. Modality names are converted to the federated names. Splicing phenotype
    IDs (gene__chrom:start:end:clu_N_strand) are rewritten to the frozen cluster IDs so
    they match what every site produces. Intron-retention IDs come from MAJIQ in the
    atlas and from IRFinder-S at the sites; they are kept as atlas IDs here and matched
    at the site by intron coordinates (feature_key = chrom:start:end, when parseable).
"""
from __future__ import annotations

import argparse
import gzip
import os
import re
import sys

import pandas as pd

# atlas modality label -> federated modality name
MODALITIES = {
    "expression": "expression",
    "isoform_expression": "isoform_expression",
    "isoforms": "isoform_usage",
    "splicing": "splicing",
    "intron_retention": "intron_retention",
    "alt_TSS": "alt_tss",
    "alt_polyA": "alt_polya",
    "RNA_editing": "rna_editing",
    "stability": "stability",
}
JUNC = re.compile(r"^(?P<chrom>[^:]+):(?P<start>\d+):(?P<end>\d+):clu_(?P<clu>\d+)_(?P<strand>[+\-.NA]+)$")
COORD = re.compile(r"(chr[0-9XYM]+)[_:](\d+)[-_:](\d+)")


def chrom_key(c: str):
    c = c[3:] if c.startswith("chr") else c
    return (0, int(c), "") if c.isdigit() else (1, 0, c)


# ---------------------------------------------------------------------------- clusters
def read_junctions(path: str):
    """Yield (chrom, start, end, strand, clu_n) from a leafCutter perind counts file."""
    with gzip.open(path, "rt") as fh:
        fh.readline()  # header: sample names
        for line in fh:
            rid = line.split(" ", 1)[0].split("\t", 1)[0]
            m = JUNC.match(rid)
            if not m:
                raise ValueError(f"{path}: unexpected junction id {rid!r}")
            yield m["chrom"], int(m["start"]), int(m["end"]), m["strand"], int(m["clu"])


def frozen_clusters(counts: dict[str, str]) -> pd.DataFrame:
    parent: dict = {}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    junc_nodes: dict[tuple, list] = {}
    for anc, path in counts.items():
        for chrom, s, e, strand, clu in read_junctions(path):
            node = (anc, clu)
            parent.setdefault(node, node)
            junc_nodes.setdefault((chrom, s, e, strand), []).append(node)
    for nodes in junc_nodes.values():
        for n in nodes[1:]:
            union(nodes[0], n)
    comp_juncs: dict = {}
    for j, nodes in junc_nodes.items():
        comp_juncs.setdefault(find(nodes[0]), []).append(j)
    comps = sorted(comp_juncs.values(), key=lambda js: (chrom_key(min(js)[0]), min(j[1] for j in js)))
    rows = []
    for k, js in enumerate(comps, start=1):
        for chrom, s, e, strand in sorted(js, key=lambda j: (j[1], j[2])):
            src = sorted({f"{a}:clu_{c}" for a, c in junc_nodes[(chrom, s, e, strand)]})
            rows.append((chrom, s, e, strand, f"clu_{k}", ";".join(src)))
    df = pd.DataFrame(rows, columns=["chrom", "start", "end", "strand", "cluster_id", "atlas_clusters"])
    return df.sort_values(["chrom", "start", "end"], key=lambda s: s.map(chrom_key) if s.name == "chrom" else s,
                          kind="stable").reset_index(drop=True)


# ---------------------------------------------------------------------------- signals
def _ids(df: pd.DataFrame) -> pd.DataFrame:
    if "phenotype_id" not in df.columns:
        df = df.reset_index().rename(columns={"index": "phenotype_id"})
    return df


def read_signals(results: str, anc: str, label: str, qmax: float) -> pd.DataFrame | None:
    lead_p = os.path.join(results, f"{anc}_{label}_cisqtl.parquet")
    ind_p = os.path.join(results, f"{anc}_{label}_cisqtl_independent.parquet")
    parts = []
    if os.path.exists(ind_p):
        ind = _ids(pd.read_parquet(ind_p))
        if len(ind):
            rank = ind["rank"] if "rank" in ind.columns else 1
            parts.append(pd.DataFrame({"atlas_phenotype_id": ind["phenotype_id"], "variant_id": ind["variant_id"],
                                       "signal_rank": rank}))
    if os.path.exists(lead_p):
        lead = _ids(pd.read_parquet(lead_p))
        if "qval" in lead.columns:
            lead = lead[lead["qval"] <= qmax]
            parts.append(pd.DataFrame({"atlas_phenotype_id": lead["phenotype_id"],
                                       "variant_id": lead["variant_id"], "signal_rank": 1}))
    if not parts:
        return None
    df = pd.concat(parts, ignore_index=True)
    df["signal_rank"] = df["signal_rank"].astype(int)
    return (df.sort_values("signal_rank").drop_duplicates(["atlas_phenotype_id", "variant_id"])
              .reset_index(drop=True))


def canonical_splicing(pid: str, jmap: dict) -> str:
    gene, _, junc = pid.rpartition("__")
    m = JUNC.match(junc)
    if not m:
        raise ValueError(f"unexpected splicing phenotype id {pid!r}")
    key = (m["chrom"], int(m["start"]), int(m["end"]), m["strand"])
    if key not in jmap:
        raise ValueError(f"splicing signal {pid!r}: junction not in the frozen cluster set")
    new = f"{m['chrom']}:{m['start']}:{m['end']}:{jmap[key]}_{m['strand']}"
    return f"{gene}__{new}" if gene else new


def feature_key(pid: str) -> str:
    m = COORD.search(pid)
    return f"{m[1]}:{m[2]}:{m[3]}" if m else ""


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ancestries", required=True, help="space-separated, e.g. 'EAS EUR'")
    ap.add_argument("--splicing-counts", required=True, help="path with {ANC} placeholder")
    ap.add_argument("--qtl-results", required=True)
    ap.add_argument("--qvalue-max", type=float, default=0.05)
    ap.add_argument("--include-ungrouped", action="store_true")
    ap.add_argument("--out-clusters", required=True)
    ap.add_argument("--out-signals", required=True)
    a = ap.parse_args(argv)
    ancs = a.ancestries.split()

    counts = {anc: a.splicing_counts.replace("{ANC}", anc) for anc in ancs}
    for anc, p in counts.items():
        if not os.path.exists(p):
            print(f"missing splicing counts for {anc}: {p}", file=sys.stderr)
            return 1
    fc = frozen_clusters(counts)
    jmap = {(r.chrom, r.start, r.end, r.strand): r.cluster_id for r in fc.itertuples()}
    print(f"frozen clusters: {fc['cluster_id'].nunique()} clusters, {len(fc)} junctions "
          f"(from {', '.join(ancs)})")

    out = []
    for anc in ancs:
        for label, mod in MODALITIES.items():
            for scan, lab in [("grouped", label)] + ([("ungrouped", f"{label}_ungrouped")] if a.include_ungrouped else []):
                df = read_signals(a.qtl_results, anc, lab, a.qvalue_max)
                if df is None:
                    print(f"  {anc} {lab}: no results found", file=sys.stderr)
                    continue
                df["modality"], df["ancestry"], df["scan"] = mod, anc, scan
                if mod == "splicing":
                    df["phenotype_id"] = [canonical_splicing(p, jmap) for p in df["atlas_phenotype_id"]]
                else:
                    df["phenotype_id"] = df["atlas_phenotype_id"]
                df["feature_key"] = [feature_key(p) if mod == "intron_retention" else "" for p in df["atlas_phenotype_id"]]
                print(f"  {anc} {lab}: {len(df)} signals")
                out.append(df)
    if not out:
        print("no signals found", file=sys.stderr)
        return 1
    sig = pd.concat(out, ignore_index=True)
    sig = (sig.sort_values(["modality", "ancestry", "phenotype_id", "signal_rank", "scan"])
              .drop_duplicates(["phenotype_id", "variant_id", "modality", "ancestry"]))
    cols = ["phenotype_id", "variant_id", "modality", "ancestry", "signal_rank", "atlas_phenotype_id", "scan", "feature_key"]
    os.makedirs(os.path.dirname(os.path.abspath(a.out_signals)), exist_ok=True)
    # gzip with mtime 0 so the archive is byte-reproducible
    for df, path in ((fc, a.out_clusters), (sig[cols], a.out_signals)):
        with open(path, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0, filename="") as gz:
            gz.write(df.to_csv(sep="\t", index=False).encode())
    print(f"MDACC signal list: {len(sig)} signal(s), {sig['phenotype_id'].nunique()} features, "
          f"{sig.groupby(['modality', 'ancestry']).ngroups} modality x ancestry sets")
    return 0


if __name__ == "__main__":
    sys.exit(main())
