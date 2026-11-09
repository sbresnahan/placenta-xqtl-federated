#!/usr/bin/env python3
"""xqtl_signals.py -- signal test sets for G x C and G x E tier 2, per stratum x modality.

  xqtl_signals.py --anc EUR --modality expression --cis EUR_expression.cis.parquet \
      --bed EUR_expression.bed.gz --pvar EUR.pvar --signals mdacc_signals.tsv.gz \
      [--fdr 0.05] --out EUR_expression

--cis is the stratum's lead table: the expression permutation table for expression, the
ungrouped (feature-level) permutation table for the other modalities. Its top variant is the site lead of each phenotype with
qval <= --fdr.

Writes
  OUT.gx_pairs.tsv       G x C pairs: site leads + MDACC pairs with signal_rank 1 (any discovery
                         ancestry, kept in mdacc_ancestry). source = site | mdacc | both (same
                         phenotype_id and variant_id). untested_reason = phenotype_absent (not in
                         the stratum BED) | variant_absent (not in the stratum genotypes) | empty.
  OUT.gx_phenotypes.tsv  tier-2 phenotypes: site leads + every phenotype in the MDACC list (any
                         rank), with the site lead and MDACC lead variants; phenotype_absent rows
                         are kept as untested.

Intron retention: the MDACC list uses MAJIQ segments (GENE__IR_GENE_chr_s_e_strand_n), the
site IRFinder-S introns (GENE__IR_chr_s_e_strand). An MDACC segment maps to the shortest site
intron of the same gene, chromosome and strand that contains it (site start <= segment start and
site end >= segment end - 1); the original ID is kept in mdacc_phenotype_id.
"""
import argparse
import gzip
import re

import numpy as np
import pandas as pd

IR_COORD = re.compile(r"(chr[^_]+)_(\d+)_(\d+)_([+-])")
PAIR_COLS = ["ancestry", "modality", "phenotype_id", "variant_id", "source", "site_lead", "site_qval",
             "mdacc_lead", "mdacc_ancestry", "mdacc_phenotype_id", "untested_reason"]
PHEN_COLS = ["ancestry", "modality", "phenotype_id", "source", "site_lead_variant", "site_qval",
             "mdacc_lead_variants", "mdacc_min_rank", "mdacc_phenotype_id", "untested_reason"]


def bed_ids(path):
    ids = []
    with gzip.open(path, "rt") if path.endswith(".gz") else open(path) as fh:
        next(fh)
        for line in fh:
            ids.append(line.split("\t", 4)[3])
    return ids


def pvar_ids(path):
    ids = set()
    with open(path) as fh:
        hdr = None
        for line in fh:
            if line.startswith("##"):
                continue
            t = line.rstrip("\n").split("\t")
            if hdr is None and line.startswith("#"):
                hdr = [h.lstrip("#") for h in t]
                continue
            ids.add(t[hdr.index("ID")] if hdr else t[1])
    return ids


def ir_key(pid):
    m = IR_COORD.search(pid)
    if not m:
        return None
    return pid.split("__", 1)[0], m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)


def map_ir(mdacc_ids, site_ids):
    """MDACC IR segment -> shortest containing site intron (same gene, chrom, strand), or None."""
    site = [(s, ir_key(s)) for s in site_ids]
    site = [(s, k) for s, k in site if k]
    out = {}
    for m in mdacc_ids:
        k = ir_key(m)
        hit = None
        if k:
            g, c, s0, e0, st = k
            cand = [(e - s, sid) for sid, (gg, cc, s, e, ss) in site
                    if gg == g and cc == c and ss == st and s <= s0 and e >= e0 - 1]
            if cand:
                hit = min(cand)[1]
        out[m] = hit
    return out


def build(anc, mod, cis, site_ids, variants, signals, fdr=0.05):
    """Return (pairs, phenotypes) data frames. cis: permutation lead table (may be empty)."""
    site_set = set(site_ids)
    sig = signals[signals["modality"] == mod].copy()
    sig["mdacc_phenotype_id"] = sig["phenotype_id"]
    if mod == "intron_retention" and len(sig):
        mp = map_ir(sig["phenotype_id"].unique(), site_ids)
        sig["phenotype_id"] = [mp.get(p) or p for p in sig["mdacc_phenotype_id"]]
        sig["mapped"] = [mp.get(p) is not None for p in sig["mdacc_phenotype_id"]]
    else:
        sig["mapped"] = sig["phenotype_id"].isin(site_set)
    sig["signal_rank"] = pd.to_numeric(sig["signal_rank"], errors="coerce")

    leads = cis[pd.to_numeric(cis["qval"], errors="coerce") <= fdr][["phenotype_id", "variant_id", "qval"]] \
        if len(cis) else pd.DataFrame(columns=["phenotype_id", "variant_id", "qval"])
    leads = leads.drop_duplicates("phenotype_id")

    # ---- G x C pairs
    pairs = {}
    for r in leads.itertuples():
        pairs[(r.phenotype_id, r.variant_id)] = {"site_lead": True, "site_qval": float(r.qval),
                                                 "mdacc_lead": False, "mdacc_ancestry": "", "mdacc_phenotype_id": ""}
    for (pid, vid), sub in sig[sig["signal_rank"] == 1].groupby(["phenotype_id", "variant_id"], sort=False):
        rec = pairs.setdefault((pid, vid), {"site_lead": False, "site_qval": np.nan, "mdacc_lead": False,
                                            "mdacc_ancestry": "", "mdacc_phenotype_id": ""})
        rec["mdacc_lead"] = True
        rec["mdacc_ancestry"] = ",".join(sorted(set(sub["ancestry"].astype(str))))
        rec["mdacc_phenotype_id"] = ",".join(sorted(set(sub["mdacc_phenotype_id"].astype(str))))
    rows = []
    for (pid, vid), r in pairs.items():
        src = "both" if r["site_lead"] and r["mdacc_lead"] else ("site" if r["site_lead"] else "mdacc")
        reason = "" if pid in site_set else "phenotype_absent"
        if not reason and vid not in variants:
            reason = "variant_absent"
        rows.append({"ancestry": anc, "modality": mod, "phenotype_id": pid, "variant_id": vid, "source": src,
                     **r, "untested_reason": reason})
    pairs_df = pd.DataFrame(rows, columns=PAIR_COLS).sort_values(["phenotype_id", "variant_id"]).reset_index(drop=True)

    # ---- tier-2 phenotypes
    ph = {}
    for r in leads.itertuples():
        ph[r.phenotype_id] = {"site_lead_variant": r.variant_id, "site_qval": float(r.qval),
                              "mdacc_lead_variants": "", "mdacc_min_rank": np.nan, "mdacc_phenotype_id": ""}
    for pid, sub in sig.groupby("phenotype_id", sort=False):
        rec = ph.setdefault(pid, {"site_lead_variant": "", "site_qval": np.nan, "mdacc_lead_variants": "",
                                  "mdacc_min_rank": np.nan, "mdacc_phenotype_id": ""})
        rec["mdacc_lead_variants"] = ",".join(sorted(set(sub.loc[sub["signal_rank"] == 1, "variant_id"])))
        rec["mdacc_min_rank"] = float(sub["signal_rank"].min())
        rec["mdacc_phenotype_id"] = ",".join(sorted(set(sub["mdacc_phenotype_id"].astype(str))))
    prow = []
    for pid, r in ph.items():
        site = bool(r["site_lead_variant"])
        mdacc = bool(r["mdacc_phenotype_id"])
        src = "both" if site and mdacc else ("site" if site else "mdacc")
        prow.append({"ancestry": anc, "modality": mod, "phenotype_id": pid, "source": src, **r,
                     "untested_reason": "" if pid in site_set else "phenotype_absent"})
    phen_df = pd.DataFrame(prow, columns=PHEN_COLS).sort_values("phenotype_id").reset_index(drop=True)
    return pairs_df, phen_df


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--anc", required=True)
    p.add_argument("--modality", required=True)
    p.add_argument("--cis", required=True)
    p.add_argument("--bed", required=True)
    p.add_argument("--pvar", required=True)
    p.add_argument("--signals", required=True)
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--out", required=True)
    a = p.parse_args(argv)
    cis = pd.read_parquet(a.cis)
    sig = pd.read_csv(a.signals, sep="\t", dtype={"phenotype_id": str, "variant_id": str, "modality": str,
                                                    "ancestry": str})
    pairs, phen = build(a.anc, a.modality, cis, bed_ids(a.bed), pvar_ids(a.pvar), sig, a.fdr)
    pairs.to_csv(a.out + ".gx_pairs.tsv", sep="\t", index=False)
    phen.to_csv(a.out + ".gx_phenotypes.tsv", sep="\t", index=False)
    print(f"{a.out}: {len(pairs)} G x C pairs ({(pairs.untested_reason == '').sum()} testable), "
          f"{len(phen)} tier-2 phenotypes ({(phen.untested_reason == '').sum()} present)")


if __name__ == "__main__":
    main()
