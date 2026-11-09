#!/usr/bin/env python3
"""Central (MDACC) meta-analysis of site exports (--mode central). See docs/CENTRAL_ANALYSIS.md.

A *contribution* is one study x ancestry stratum (an export folder <ANC>/ with status 'mapped' in the
study's I_strata table). Contributions are numbered in contributions.tsv; `contrib_mask` columns are
bit masks over that index (bit i set = contribution i contributed).

Subcommands
  ingest         MANIFEST checksums, versions, file schemas, contributions, allele harmonization map
                 (per chromosome) with palindromic drops and frequency-outlier flags
  meta-nominal   fixed- and random-effects IVW with Q / I^2 of file A, per modality x chromosome,
                 within each ancestry and across all contributions (scope ALL)
  features       ACAT of per-contribution pval_beta (file B permutation tables) per arm and scope,
                 Storey q (BH fallback), meta lead variant per feature
  finemap-prep   per significant feature: within-ancestry meta z-scores and N-weighted LD (file C)
                 on a shared variant set -> .npz inputs for xqtl_finemap.py (SuShiE summary mode)
  coloc-prep     xQTL x GWAS (catalog by path) inputs for xqtl_coloc.R (SuSiE-coloc, colocBoost)
  gxc            G x C pair x cell-type IVW (ACAT / single p otherwise), BH, central primary cell type
  gxe-tier1      tier-1 GxE synthesis across contributions (atlas 54_tier1_meta), scope labels
  gxe-tier2      tier-2 variant-level IVW after rescaling to the pooled exposure SD, top variant per
                 phenotype, BH within modality x exposure
  twas-inputs    OTTERS / mvOTTERS stage-I inputs (annotation + eQTL summary statistics per ancestry)
                 and FUSION-format weight tables from SuShiE weights (written to .wgt.RDat by
                 xqtl_fusion_wgt.R)
  concat         concatenate shard tables (fine-mapping, colocalization)
  summary        aggregation tables and PIP enrichment (Fisher, PIP-weighted) for annotation BEDs
"""
from __future__ import annotations

import argparse
import glob
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd
from scipy import stats

COMP = str.maketrans("ACGT", "TGCA")
NOVEL_RE = re.compile(r"nclu_")


def log(msg):
    print(msg, file=sys.stderr, flush=True)


# ---------------------------------------------------------------------------------------------
# statistics
# ---------------------------------------------------------------------------------------------
def acat(p):
    """Cauchy combination of p-values (equal weights); NaNs ignored (atlas 54_tier1_meta)."""
    p = np.asarray(p, dtype=float)
    p = p[np.isfinite(p)]
    if len(p) == 0:
        return np.nan
    if len(p) == 1:
        return float(p[0])
    p = np.clip(p, 1e-300, 1 - 1e-15)
    small = p < 1e-15
    terms = np.where(small, 1.0 / (p * np.pi), np.tan((0.5 - p) * np.pi))
    t = terms.mean()
    if t > 1e15:
        return float(1.0 / (t * np.pi))
    return float(stats.cauchy.sf(t))


def bh(p):
    p = np.asarray(p, dtype=float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return out
    x = p[ok]
    o = np.argsort(x)
    q = x[o] * len(x) / np.arange(1, len(x) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    tmp = np.empty(len(x))
    tmp[o] = np.clip(q, 0, 1)
    out[ok] = tmp
    return out


def storey(p):
    """Storey q-values via bin/compute_qvalues.R (qvalue, default lambda); BH if that fails."""
    p = np.asarray(p, dtype=float)
    ok = np.isfinite(p)
    out = np.full(len(p), np.nan)
    if ok.sum() == 0:
        return out, "none"
    bridge = shutil.which("compute_qvalues.R") or os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                               "compute_qvalues.R")
    with tempfile.TemporaryDirectory() as td:
        fi, fo = os.path.join(td, "in.tsv"), os.path.join(td, "out.tsv")
        pd.DataFrame({"pval_beta": p[ok]}).to_csv(fi, sep="\t", index=False)
        try:
            subprocess.run(["Rscript", bridge, fi, fo], check=True, capture_output=True, timeout=600)
            q = pd.read_csv(fo, sep="\t")["qval"].to_numpy(float)
            if len(q) == ok.sum() and np.isfinite(q).all():
                out[ok] = q
                return out, "storey"
        except Exception as e:  # noqa: BLE001  (fallback documented in the output column)
            log(f"Storey q failed ({e}); using BH")
    out[ok] = bh(p[ok])
    return out, "bh"


def ivw(gid, beta, se, ngroups):
    """Grouped fixed- and random-effects (DerSimonian-Laird) IVW. gid: group index per row."""
    beta = np.asarray(beta, float)
    se = np.asarray(se, float)
    w = 1.0 / se ** 2
    k = np.bincount(gid, minlength=ngroups).astype(float)
    sw = np.bincount(gid, w, ngroups)
    swb = np.bincount(gid, w * beta, ngroups)
    swbb = np.bincount(gid, w * beta ** 2, ngroups)
    sww = np.bincount(gid, w ** 2, ngroups)
    with np.errstate(divide="ignore", invalid="ignore"):
        b = swb / sw
        s = np.sqrt(1.0 / sw)
        q = np.clip(swbb - swb ** 2 / sw, 0, None)
        dfq = k - 1
        i2 = np.where(q > 0, np.clip((q - dfq) / q, 0, None), 0.0)
        p_het = stats.chi2.sf(q, np.maximum(dfq, 1))
        c = sw - sww / sw
        tau2 = np.where((k >= 2) & (c > 0), np.clip((q - dfq) / c, 0, None), 0.0)
        multi = k >= 2
        i2[~multi] = np.nan
        p_het[~multi] = np.nan
        q[~multi] = np.nan
        wr = 1.0 / (se ** 2 + tau2[gid])
        swr = np.bincount(gid, wr, ngroups)
        br = np.bincount(gid, wr * beta, ngroups) / swr
        sr = np.sqrt(1.0 / swr)
    z = b / s
    return {"k": k.astype(int), "beta_fe": b, "se_fe": s, "z_fe": z, "p_fe": 2 * stats.norm.sf(np.abs(z)),
            "beta_re": br, "se_re": sr, "p_re": 2 * stats.norm.sf(np.abs(br / sr)),
            "tau2": tau2, "q_het": q, "p_het": p_het, "i2": i2}


def group_index(df, keys):
    """(gid per row, first-row frame of each group in gid order)."""
    gid = df.groupby(keys, sort=False, dropna=False).ngroup().to_numpy()
    first = df.drop_duplicates(keys).reset_index(drop=True)
    return gid, first


# ---------------------------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------------------------
def read_any(stem, columns=None, filters=None):
    """Read <stem>.parquet or <stem>.tsv.gz; None if neither exists."""
    if os.path.exists(stem + ".parquet"):
        return pd.read_parquet(stem + ".parquet", columns=columns, filters=filters)
    if os.path.exists(stem + ".tsv.gz"):
        d = pd.read_csv(stem + ".tsv.gz", sep="\t", usecols=columns)
        if filters:
            for col, _, val in filters:
                d = d[d[col] == val]
        return d
    return None


def table_columns(stem):
    if os.path.exists(stem + ".parquet"):
        import pyarrow.parquet as pq
        return pq.read_schema(stem + ".parquet").names
    if os.path.exists(stem + ".tsv.gz"):
        return list(pd.read_csv(stem + ".tsv.gz", sep="\t", nrows=0).columns)
    return None


def write(df, path):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    if path.endswith(".parquet"):
        df.to_parquet(path, index=False)
    else:
        df.to_csv(path, sep="\t", index=False, na_rep="NA")


def cfg_get(cfg, dotted, default=None):
    d = cfg
    for k in dotted.split("."):
        if not isinstance(d, dict) or k not in d:
            return default
        d = d[k]
    return default if d is None else d


def load_cfg(path):
    with open(path) as fh:
        return json.load(fh)


def contributions(path):
    c = pd.read_csv(path, sep="\t", dtype={"study": str, "cohort": str, "ancestry": str})
    c["prefix"] = [os.path.join(r.export_dir, r.ancestry, f"{r.cohort}_{r.ancestry}") for r in c.itertuples()]
    return c


def harm_map(harm_dir, chrom=None):
    files = [os.path.join(harm_dir, f"harm_{chrom}.parquet")] if chrom else \
        sorted(glob.glob(os.path.join(harm_dir, "harm_*.parquet")))
    files = [f for f in files if os.path.exists(f)]
    if not files:
        return pd.DataFrame(columns=["contrib", "variant_id", "canonical_id", "sign", "freq_outlier"])
    return pd.concat([pd.read_parquet(f, columns=["contrib", "variant_id", "canonical_id", "sign", "freq_outlier"])
                      for f in files], ignore_index=True)


def apply_harm(df, hm, contrib, beta_cols=(), eaf_cols=()):
    """Join a contribution's rows to the harmonization map; drop excluded variants, flip signs."""
    h = hm[hm["contrib"] == contrib][["variant_id", "canonical_id", "sign", "freq_outlier"]]
    d = df.merge(h, on="variant_id", how="inner")
    d = d[d["sign"] != 0].copy()
    neg = d["sign"] < 0
    for c in beta_cols:
        if c in d.columns:
            d.loc[neg, c] = -d.loc[neg, c]
    for c in eaf_cols:
        if c in d.columns:
            d.loc[neg, c] = 1 - d.loc[neg, c]
    d["variant_id"] = d.pop("canonical_id")
    return d


# ---------------------------------------------------------------------------------------------
# ingest
# ---------------------------------------------------------------------------------------------
REQUIRED = {
    "A": ["phenotype_id", "variant_id", "tss_distance", "n", "beta", "se", "pval"],
    "B_permutation": ["phenotype_id", "variant_id", "pval_beta", "slope", "slope_se"],
    "D": ["variant_id", "chrom", "pos", "ref", "alt", "eaf", "n"],
    "E": ["modality", "phenotype_id", "chr", "start", "end", "novel"],
    "G": ["modality", "phenotype_id", "variant_id", "cell_type", "status", "beta_GxC", "se_GxC", "p_GxC"],
    "H_tier1": ["modality", "exposure", "phenotype_id", "variant_id", "b_gi", "b_gi_se", "pval_beta"],
    "H_tier2": ["phenotype_id", "variant_id", "n", "b_gi", "b_gi_se", "pval_nominal"],
}


class Report:
    def __init__(self):
        self.rows = []

    def add(self, check, status, detail=""):
        self.rows.append({"check": check, "status": status, "detail": detail})
        if status == "FAIL":
            log(f"FAIL {check}: {detail}")

    @property
    def failed(self):
        return any(r["status"] == "FAIL" for r in self.rows)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def harmonize(d, pal_maf_max, freq_diff_max):
    """d: contrib, ancestry, variant_id, chrom, pos, ref, alt, eaf, n. Returns the map rows."""
    d = d.copy()
    ref, alt = d["ref"].astype(str), d["alt"].astype(str)
    snp = ref.str.fullmatch(r"[ACGT]") & alt.str.fullmatch(r"[ACGT]")
    cref = np.where(snp, ref.str.translate(COMP), ref)
    calt = np.where(snp, alt.str.translate(COMP), alt)
    a = np.where(ref < alt, ref + "/" + alt, alt + "/" + ref)
    b = np.where(cref < calt, pd.Series(cref) + "/" + pd.Series(calt), pd.Series(calt) + "/" + pd.Series(cref))
    d["pal"] = snp & (pd.Series(cref, index=d.index) == alt)
    d["ckey"] = d["chrom"].astype(str) + ":" + d["pos"].astype(str) + ":" + np.minimum(a, b)
    d["orient"] = ref + ">" + alt
    d["order"] = np.arange(len(d))
    cnt = d.groupby(["ckey", "orient"], sort=False).agg(c=("order", "size"), o=("order", "min")).reset_index()
    cnt = cnt.sort_values(["ckey", "c", "o"], ascending=[True, False, True]).drop_duplicates("ckey")
    d = d.merge(cnt[["ckey", "orient"]].rename(columns={"orient": "canon"}), on="ckey", how="left")
    can = d["canon"].str.split(">", expand=True)
    c_ref, c_alt = can[0], can[1]
    ref, alt = d["ref"].astype(str), d["alt"].astype(str)
    cref = ref.str.translate(COMP)
    calt = alt.str.translate(COMP)
    keep = (ref == c_ref) & (alt == c_alt)
    swap = (ref == c_alt) & (alt == c_ref)
    flip = ~keep & ~swap & (cref == c_ref) & (calt == c_alt)
    flip_swap = ~keep & ~swap & ~flip & (cref == c_alt) & (calt == c_ref)
    d["action"] = np.select([keep, swap, flip, flip_swap], ["keep", "swap", "flip", "flip_swap"], "mismatch")
    d["sign"] = np.select([keep | flip, swap | flip_swap], [1, -1], 0)
    eaf = d["eaf"].astype(float)
    maf = np.minimum(eaf, 1 - eaf)
    # Palindromic (A/T, C/G) sites: exports carry GRCh38 forward-strand IDs, so a site with the same
    # orientation in every contribution is unambiguous and kept. Where orientations differ, a swap cannot
    # be told from a strand flip; the site is dropped from every contribution when its MAF exceeds the
    # limit in any of them (frequencies cannot resolve it), else read literally (swap).
    disagree = set(d.loc[d["pal"] & ~keep, "ckey"])
    high = set(d.loc[d["pal"] & (maf > pal_maf_max), "ckey"])
    pal_drop = d["pal"] & d["ckey"].isin(disagree & high)
    d.loc[pal_drop, "action"] = "drop_palindromic"
    d.loc[pal_drop, "sign"] = 0
    d["canonical_id"] = d["chrom"].astype(str) + ":" + d["pos"].astype(str) + ":" + c_ref + ":" + c_alt
    d["eaf_h"] = np.where(d["sign"] < 0, 1 - eaf, eaf)
    # frequency outliers: |eaf - N-weighted mean eaf of the other contributions of the ancestry|
    d["freq_outlier"] = False
    d["eaf_others"] = np.nan
    ok = d["sign"] != 0
    sub = d[ok]
    g = sub.groupby(["ancestry", "canonical_id"])
    ntot = g["n"].transform("sum").astype(float)
    sne = (sub["n"] * sub["eaf_h"]).groupby([sub["ancestry"], sub["canonical_id"]]).transform("sum")
    k = g["n"].transform("size")
    with np.errstate(divide="ignore", invalid="ignore"):
        other = (sne - sub["n"] * sub["eaf_h"]) / (ntot - sub["n"])
    other = other.where(k >= 2)
    d.loc[ok, "eaf_others"] = other.to_numpy()
    d.loc[ok, "freq_outlier"] = (np.abs(sub["eaf_h"] - other) > freq_diff_max).fillna(False).to_numpy()
    return d[["contrib", "ancestry", "variant_id", "canonical_id", "action", "sign", "eaf_h", "eaf_others",
              "freq_outlier", "n"]]


def cmd_ingest(a):
    cfg = load_cfg(a.config)
    rep = Report()
    studies = pd.read_csv(a.studies, sep="\t", dtype=str)
    if not {"study_id", "export_dir"} <= set(studies.columns):
        sys.exit("central.studies must have columns study_id, export_dir")
    if studies["study_id"].duplicated().any():
        sys.exit("duplicate study_id in central.studies")
    verify = bool(cfg_get(cfg, "central.verify_checksums", True))
    rows, versions = [], []
    for i, s in enumerate(studies.itertuples()):
        ed = s.export_dir
        if a.staged_prefix:
            # Nextflow stages the i-th study as <prefix><i+1>/<dir>, zero-padded by the stageAs pattern
            cands = [os.path.join(a.staged_prefix + str(i + 1).zfill(w), os.path.basename(s.export_dir.rstrip("/")))
                     for w in (1, 2, 3)]
            ed = next((c for c in cands if os.path.isdir(c)), cands[0])
        man_p = os.path.join(ed, "MANIFEST.tsv")
        if not os.path.exists(man_p):
            rep.add(f"manifest.{s.study_id}", "FAIL", f"no MANIFEST.tsv in {s.export_dir}")
            continue
        man = pd.read_csv(man_p, sep="\t", dtype=str)
        base = os.path.dirname(os.path.abspath(ed))
        bad, missing = [], []
        for r in man.itertuples():
            fp = os.path.join(base, r.path)
            if not os.path.exists(fp):
                missing.append(r.path)
            elif verify and sha256(fp) != r.sha256:
                bad.append(r.path)
        if missing or bad:
            rep.add(f"manifest.{s.study_id}", "FAIL", f"missing {missing[:3]} / checksum mismatch {bad[:3]} "
                    f"({len(missing)} missing, {len(bad)} mismatched of {len(man)})")
        else:
            rep.add(f"manifest.{s.study_id}", "PASS", f"{len(man)} files" + (", sha256 verified" if verify else ""))
        for col in ("pipeline_version", "restricted_version", "image_sha256"):
            for v in man[col].dropna().unique():
                versions.append({"study": s.study_id, "item": col, "value": v})
        strata = glob.glob(os.path.join(ed, "*_cohort_I_strata.tsv"))
        if not strata:
            rep.add(f"strata.{s.study_id}", "FAIL", "no *_cohort_I_strata.tsv")
            continue
        st = pd.read_csv(strata[0], sep="\t", dtype=str)
        for r in st.itertuples():
            if r.status != "mapped":
                rep.add(f"strata.{s.study_id}.{r.ancestry}", "INFO", f"not used: status {r.status}")
                continue
            rows.append({"study": s.study_id, "cohort": r.cohort, "ancestry": r.ancestry,
                         "n_individuals": int(r.n_individuals), "n_rnaseq": int(r.n_rnaseq),
                         "gxc": getattr(r, "gxc", ""), "gxe": getattr(r, "gxe", ""),
                         "export_dir": ed, "source_dir": s.export_dir})
    ver = pd.DataFrame(versions, columns=["study", "item", "value"])
    for item in ("pipeline_version", "restricted_version"):
        vals = ver[ver["item"] == item].groupby("study")["value"].agg(lambda x: ";".join(sorted(set(x))))
        if item == "pipeline_version":
            majors = {str(v).split(".")[0] for v in vals}
            ok = len(majors) <= 1 and all(";" not in v for v in vals)
        else:
            ok = vals.nunique() <= 1 and all(";" not in v for v in vals)
        rep.add(f"versions.{item}", "PASS" if ok else "FAIL",
                ", ".join(f"{k}={v}" for k, v in vals.items()) +
                ("" if ok else " (exports must share the export format: same major pipeline version and "
                 "restricted-file version)"))
    img = ver[ver["item"] == "image_sha256"]["value"].unique()
    if len(img) > 1:
        rep.add("versions.image_sha256", "WARN", f"{len(img)} different images: " + ", ".join(sorted(img)))
    con = pd.DataFrame(rows)
    if con.empty:
        rep.add("contributions", "FAIL", "no mapped study x ancestry strata")
    else:
        con.insert(0, "contrib", np.arange(len(con)))
        rep.add("contributions", "PASS", f"{len(con)} study x ancestry contributions: " +
                ", ".join(f"{r.study}/{r.ancestry} (n={r.n_individuals})" for r in con.itertuples()))
    os.makedirs(a.out_dir, exist_ok=True)

    # schemas + units
    units = set()
    if not con.empty:
        con["prefix"] = [os.path.join(r.export_dir, r.ancestry, f"{r.cohort}_{r.ancestry}") for r in con.itertuples()]
        for r in con.itertuples():
            pre, d = r.prefix, os.path.dirname(r.prefix)
            checks = {"D": [pre + "_D_variant_qc"], "E": [pre + "_E_features"]}
            checks["A"] = [p[:-len(".parquet")] if p.endswith(".parquet") else p[:-len(".tsv.gz")]
                           for p in sorted(glob.glob(os.path.join(d, "A_nominal", "*_A_nominal_*")))]
            checks["B_permutation"] = [p.rsplit(".", 2 if p.endswith(".gz") else 1)[0]
                                       for p in sorted(glob.glob(os.path.join(d, "B_mapping", "*_B_permutation_*")))]
            if os.path.exists(pre + "_G_gxc.parquet") or os.path.exists(pre + "_G_gxc.tsv.gz"):
                checks["G"] = [pre + "_G_gxc"]
            if os.path.exists(pre + "_H_gxe_tier1.parquet") or os.path.exists(pre + "_H_gxe_tier1.tsv.gz"):
                checks["H_tier1"] = [pre + "_H_gxe_tier1"]
            checks["H_tier2"] = [p.rsplit(".", 2 if p.endswith(".gz") else 1)[0]
                                 for p in sorted(glob.glob(os.path.join(d, "H_gxe_tier2", "*_H_gxe_tier2_*")))]
            for cls, stems in checks.items():
                if cls in ("A", "B_permutation", "D", "E") and not stems:
                    rep.add(f"schema.{r.study}.{r.ancestry}.{cls}", "FAIL", "no files")
                for stem in stems:
                    cols = table_columns(stem)
                    miss = [c for c in REQUIRED[cls] if cols is None or c not in cols]
                    if miss:
                        rep.add(f"schema.{r.study}.{r.ancestry}.{cls}", "FAIL",
                                f"{os.path.basename(stem)}: missing columns {miss}")
            for stem in checks["A"]:
                m = re.search(r"_A_nominal_(.+)_(chr[0-9XY]+)$", stem)
                if m:
                    units.add((m.group(1), m.group(2)))
        if not any(r["check"].startswith("schema.") for r in rep.rows):
            rep.add("schema", "PASS", "required columns present in A, B, D, E (and G, H when exported)")
    pd.DataFrame(sorted(units), columns=["modality", "chrom"]).to_csv(
        os.path.join(a.out_dir, "units.tsv"), sep="\t", index=False)

    # harmonization per chromosome
    hdir = os.path.join(a.out_dir, "harmonization")
    os.makedirs(hdir, exist_ok=True)
    summ = []
    if not con.empty and not rep.failed:
        pal = float(cfg_get(cfg, "central.palindromic_maf_max", 0.4))
        fdm = float(cfg_get(cfg, "central.freq_diff_max", 0.2))
        chroms = sorted({c for _, c in units}, key=lambda c: (len(c), c))
        for chrom in chroms:
            parts = []
            for r in con.itertuples():
                d = read_any(r.prefix + "_D_variant_qc", columns=REQUIRED["D"], filters=[("chrom", "=", chrom)])
                if d is None or d.empty:
                    continue
                d.insert(0, "contrib", r.contrib)
                d.insert(1, "ancestry", r.ancestry)
                parts.append(d)
            if not parts:
                continue
            h = harmonize(pd.concat(parts, ignore_index=True), pal, fdm)
            h.to_parquet(os.path.join(hdir, f"harm_{chrom}.parquet"), index=False)
            s = h.groupby(["contrib", "action"]).size().rename("n").reset_index()
            fo = h[h["freq_outlier"]].groupby("contrib").size().rename("n").reset_index()
            fo["action"] = "freq_outlier_flag"
            s = pd.concat([s, fo], ignore_index=True)
            s.insert(0, "chrom", chrom)
            summ.append(s)
        hs = pd.concat(summ, ignore_index=True) if summ else pd.DataFrame(columns=["chrom", "contrib", "action", "n"])
        hs = hs.groupby(["contrib", "action"], as_index=False)["n"].sum()
        hs = hs.merge(con[["contrib", "study", "ancestry"]], on="contrib")[["contrib", "study", "ancestry", "action", "n"]]
        hs.to_csv(os.path.join(a.out_dir, "harmonization_summary.tsv"), sep="\t", index=False)
        tot = hs.groupby("action")["n"].sum().to_dict()
        rep.add("harmonization", "PASS", "; ".join(f"{k} {v}" for k, v in sorted(tot.items())))
    out_con = con.drop(columns=["prefix"], errors="ignore")
    out_con.to_csv(os.path.join(a.out_dir, "contributions.tsv"), sep="\t", index=False)
    ver.to_csv(os.path.join(a.out_dir, "versions.tsv"), sep="\t", index=False)
    pd.DataFrame(rep.rows).to_csv(os.path.join(a.out_dir, "ingest_report.tsv"), sep="\t", index=False)
    if rep.failed:
        sys.exit("ERROR: central ingest failed; see ingest_report.tsv")


# ---------------------------------------------------------------------------------------------
# nominal meta-analysis (file A)
# ---------------------------------------------------------------------------------------------
def meta_frame(df, keys, extra_first=()):
    """IVW per group of `keys`. df has beta, se, n, eaf, contrib, freq_outlier."""
    gid, first = group_index(df, list(keys))
    ng = len(first)
    res = ivw(gid, df["beta"].to_numpy(), df["se"].to_numpy(), ng)
    out = first[list(keys) + list(extra_first)].copy()
    n = np.bincount(gid, df["n"].to_numpy(float), ng)
    out["n"] = n.astype(int)
    if "eaf" in df.columns:
        out["eaf"] = np.bincount(gid, (df["n"] * df["eaf"]).to_numpy(float), ng) / n
    for k, v in res.items():
        out[k] = v
    if "freq_outlier" in df.columns:
        out["n_freq_outlier"] = np.bincount(gid, df["freq_outlier"].to_numpy(float), ng).astype(int)
    mask = np.zeros(ng, dtype=np.int64)
    np.bitwise_or.at(mask, gid, (1 << df["contrib"].to_numpy(np.int64)))
    out["contrib_mask"] = mask
    return out


def scoped_meta(df, keys, extra_first=()):
    parts = []
    for anc, sub in df.groupby("ancestry", sort=True):
        m = meta_frame(sub, keys, extra_first)
        m.insert(0, "scope", anc)
        parts.append(m)
    m = meta_frame(df, keys, extra_first)
    m.insert(0, "scope", "ALL")
    parts.append(m)
    return pd.concat(parts, ignore_index=True)


def cmd_meta_nominal(a):
    con = contributions(a.contributions)
    hm = harm_map(a.harmonization, a.chrom)
    parts = []
    for r in con.itertuples():
        stem = os.path.join(os.path.dirname(r.prefix), "A_nominal", f"{r.cohort}_{r.ancestry}_A_nominal_{a.modality}_{a.chrom}")
        d = read_any(stem, columns=["phenotype_id", "variant_id", "tss_distance", "eaf", "n", "beta", "se"])
        if d is None or d.empty:
            continue
        d = d[~d["phenotype_id"].astype(str).str.contains(NOVEL_RE)]
        d = apply_harm(d, hm, r.contrib, beta_cols=("beta",), eaf_cols=("eaf",))
        d = d[np.isfinite(d["se"]) & (d["se"] > 0) & np.isfinite(d["beta"])]
        d["contrib"] = r.contrib
        d["ancestry"] = r.ancestry
        parts.append(d)
    if not parts:
        write(pd.DataFrame(), a.out)
        return
    df = pd.concat(parts, ignore_index=True)
    m = scoped_meta(df, ["phenotype_id", "variant_id"], extra_first=["tss_distance"])
    m.insert(1, "modality", a.modality)
    m.insert(2, "chrom", a.chrom)
    write(m, a.out)
    log(f"{a.modality} {a.chrom}: {len(df)} contribution rows -> {len(m)} meta rows")


# ---------------------------------------------------------------------------------------------
# features (file B permutation, ACAT + Storey)
# ---------------------------------------------------------------------------------------------
def cmd_features(a):
    cfg = load_cfg(a.config)
    fdr = float(cfg_get(cfg, "central.fdr", 0.05))
    con = contributions(a.contributions)
    hm = harm_map(a.harmonization)
    parts = []
    for r in con.itertuples():
        for p in sorted(glob.glob(os.path.join(os.path.dirname(r.prefix), "B_mapping", "*_B_permutation_*"))):
            arm = re.sub(r"\.(parquet|tsv\.gz)$", "", p).split("_B_permutation_", 1)[1]
            d = read_any(re.sub(r"\.(parquet|tsv\.gz)$", "", p),
                         columns=["phenotype_id", "variant_id", "slope", "pval_beta", "qval"])
            d = d[~d["phenotype_id"].astype(str).str.contains(NOVEL_RE)].copy()
            # pval_beta is feature-level: keep every feature; the lead ID is harmonized where possible
            h = hm[hm["contrib"] == r.contrib].set_index("variant_id")["canonical_id"]
            d["variant_id"] = d["variant_id"].map(h).fillna(d["variant_id"])
            d["arm"], d["contrib"], d["ancestry"], d["study"] = arm, r.contrib, r.ancestry, r.study
            parts.append(d)
    b = pd.concat(parts, ignore_index=True)
    arms = set(b["arm"])
    b["modality"] = b["arm"].str.replace(r"_ungrouped$", "", regex=True)
    b["grouped"] = [(x == "combined") or (f"{x}_ungrouped" in arms) for x in b["arm"]]
    b["sig_contrib"] = b["qval"] <= fdr
    rows = []
    for scope, sub in [(anc, s) for anc, s in b.groupby("ancestry")] + [("ALL", b)]:
        g = sub.groupby(["arm", "modality", "grouped", "phenotype_id"], sort=True)
        f = g.agg(k=("pval_beta", "size"), n_sig_contrib=("sig_contrib", "sum"),
                  n_studies=("study", "nunique"), min_pval_beta=("pval_beta", "min")).reset_index()
        f["p_acat"] = g["pval_beta"].apply(lambda x: acat(x.to_numpy())).to_numpy()
        f["lead_concordant"] = g["variant_id"].nunique().to_numpy() == 1
        f["contrib_mask"] = g["contrib"].apply(lambda x: int(np.bitwise_or.reduce(1 << x.to_numpy(np.int64)))).to_numpy()
        f.insert(0, "scope", scope)
        rows.append(f)
    feats = pd.concat(rows, ignore_index=True)
    feats["q_acat"], feats["qvalue_method"] = np.nan, ""
    for (scope, arm), idx in feats.groupby(["scope", "arm"]).groups.items():
        q, meth = storey(feats.loc[idx, "p_acat"].to_numpy())
        feats.loc[idx, "q_acat"] = q
        feats.loc[idx, "qvalue_method"] = meth
    feats["significant"] = feats["q_acat"] <= fdr
    # meta lead variant (file A meta) for feature-level arms
    leads = []
    for p in sorted(glob.glob(os.path.join(a.meta_dir, "*.parquet"))):
        m = pd.read_parquet(p, columns=["scope", "modality", "phenotype_id", "variant_id", "beta_fe", "se_fe",
                                        "p_fe", "p_re", "i2", "k"])
        if m.empty:
            continue
        m = m.sort_values("p_fe").drop_duplicates(["scope", "modality", "phenotype_id"])
        leads.append(m)
    if leads:
        ld = pd.concat(leads).sort_values("p_fe").drop_duplicates(["scope", "modality", "phenotype_id"])
        ld = ld.rename(columns={"variant_id": "lead_variant", "beta_fe": "lead_beta_fe", "se_fe": "lead_se_fe",
                                "p_fe": "lead_p_fe", "p_re": "lead_p_re", "i2": "lead_i2", "k": "lead_k"})
        fl = feats[~feats["grouped"]].merge(ld, on=["scope", "modality", "phenotype_id"], how="left")
        feats = pd.concat([fl, feats[feats["grouped"]]], ignore_index=True)
    feats = feats.sort_values(["arm", "scope", "p_acat"]).reset_index(drop=True)
    write(feats, a.out)
    s = feats.groupby(["arm", "scope"]).agg(tested=("phenotype_id", "size"), significant=("significant", "sum"))
    log(s.to_string())


# ---------------------------------------------------------------------------------------------
# fine-mapping inputs (SuShiE summary mode)
# ---------------------------------------------------------------------------------------------
def sanitize(x):
    return re.sub(r"[^A-Za-z0-9._+-]", "_", str(x))


def load_ld(prefix, positions, hm, contrib):
    """Pick the C window covering most of `positions`; return (canonical ids, R float64) or None."""
    d = os.path.dirname(prefix)
    wins = glob.glob(os.path.join(d, "C_ld", "*_C_ld.windows.tsv"))
    if not wins:
        return None
    w = pd.read_csv(wins[0], sep="\t")
    pos = np.asarray(positions)
    cov = [((pos >= r.start) & (pos <= r.end)).sum() for r in w.itertuples()]
    if not cov or max(cov) == 0:
        return None
    r = w.iloc[int(np.argmax(cov))]
    npz = os.path.join(d, "C_ld", os.path.basename(str(r["file"])))
    if not os.path.exists(npz):
        cand = glob.glob(os.path.join(d, "C_ld", f"*_C_ld_{r['chrom']}_{r['start']}_{r['end']}.npz"))
        if not cand:
            return None
        npz = cand[0]
    var = pd.read_csv(npz[:-4] + ".variants.tsv.gz", sep="\t", usecols=["variant_id"])
    R = np.load(npz)["R"].astype(np.float64)
    h = hm[hm["contrib"] == contrib].set_index("variant_id")
    j = var["variant_id"].map(h["canonical_id"])
    s = var["variant_id"].map(h["sign"]).fillna(0).to_numpy()
    keep = j.notna().to_numpy() & (s != 0)
    R = R[np.ix_(keep, keep)] * np.outer(s[keep], s[keep])
    return list(j[keep]), R


def cmd_finemap_prep(a):
    cfg = load_cfg(a.config)
    min_var = int(cfg_get(cfg, "central.finemap.min_variants", 50))
    ridge = float(cfg_get(cfg, "central.finemap.ld_ridge", 0.001))
    max_loci = cfg_get(cfg, "central.finemap.max_loci", None)
    con = contributions(a.contributions)
    hm = harm_map(a.harmonization)
    feats = pd.read_csv(a.features, sep="\t")
    loci = feats[(feats["scope"] == "ALL") & (~feats["grouped"]) & (feats["significant"])]
    loci = loci.sort_values("p_acat")
    if max_loci:
        loci = loci.head(int(max_loci))
    meta = {}
    for p in sorted(glob.glob(os.path.join(a.meta_dir, "*.parquet"))):
        m = pd.read_parquet(p)
        if len(m):
            for (mod, pid), sub in m[m["phenotype_id"].isin(set(loci["phenotype_id"]))].groupby(["modality", "phenotype_id"]):
                meta[(mod, pid)] = sub
    os.makedirs(a.out_dir, exist_ok=True)
    rows = []
    for f in loci.itertuples():
        lid = sanitize(f"{f.arm}__{f.phenotype_id}")
        rec = {"locus_id": lid, "arm": f.arm, "modality": f.modality, "phenotype_id": f.phenotype_id}
        m = meta.get((f.modality, f.phenotype_id))
        if m is None:
            rows.append({**rec, "status": "no_meta_rows"})
            continue
        anc_sets, anc_data = {}, {}
        for anc in sorted(set(m["scope"]) - {"ALL"}):
            ma = m[m["scope"] == anc].set_index("variant_id")
            pos = ma.index.to_series().str.split(":").str[1].astype(int).to_numpy()
            avail = set(ma.index)
            lds = []
            for r in con[con["ancestry"] == anc].itertuples():
                ld = load_ld(r.prefix, pos, hm, r.contrib)
                if ld is None:
                    avail = set()
                    break
                lds.append((r.n_individuals, ld))
                avail &= set(ld[0])
            if len(avail) >= min_var:
                anc_sets[anc] = avail
                anc_data[anc] = (ma, lds)
        common, kept = None, []
        for anc in sorted(anc_sets, key=lambda x: -len(anc_sets[x])):
            nxt = anc_sets[anc] if common is None else common & anc_sets[anc]
            if len(nxt) >= min_var:
                common, kept = nxt, kept + [anc]
        if not kept:
            rows.append({**rec, "status": "too_few_variants_with_ld"})
            continue
        kept = sorted(kept)
        vids = sorted(common, key=lambda v: int(v.split(":")[1]))
        zs, lds, ns, betas, ses, eafs = [], [], [], [], [], []
        for anc in kept:
            ma, ldl = anc_data[anc]
            sub = ma.loc[vids]
            zz, nn = sub["z_fe"].to_numpy(), float(np.median(sub["n"]))
            # PVE-adjusted z (susieR convention): marginal correlations are linear in the joint effects,
            # raw z are not when effects are large
            zs.append(zz * np.sqrt((nn - 1) / (zz ** 2 + nn - 2)))
            betas.append(sub["beta_fe"].to_numpy())
            ses.append(sub["se_fe"].to_numpy())
            eafs.append(sub["eaf"].to_numpy())
            ns.append(float(np.median(sub["n"])))
            acc, wsum = np.zeros((len(vids), len(vids))), 0.0
            for n, (ids, R) in ldl:
                ix = pd.Index(ids).get_indexer(vids)
                acc += n * R[np.ix_(ix, ix)]
                wsum += n
            R = acc / wsum
            np.fill_diagonal(R, 1.0)
            # float16 storage and averaging can leave small negative eigenvalues: shrink towards I
            lam = max(ridge, ridge - 1.1 * float(np.linalg.eigvalsh(R).min()))
            R = (R + lam * np.eye(len(vids))) / (1 + lam)
            lds.append(R.astype(np.float32))
        np.savez_compressed(os.path.join(a.out_dir, f"{lid}.npz"), variants=np.array(vids), ancestries=np.array(kept),
                            zs=np.vstack(zs), lds=np.stack(lds), ns=np.array(ns), betas=np.vstack(betas),
                            ses=np.vstack(ses), eafs=np.vstack(eafs))
        rows.append({**rec, "status": "prepared", "n_variants": len(vids), "ancestries": ";".join(kept),
                     "chrom": vids[0].split(":")[0], "start": int(vids[0].split(":")[1]),
                     "end": int(vids[-1].split(":")[1])})
    pd.DataFrame(rows, columns=["locus_id", "arm", "modality", "phenotype_id", "status", "n_variants", "ancestries",
                                "chrom", "start", "end"]).to_csv(os.path.join(a.out_dir, "loci.tsv"), sep="\t", index=False)


# ---------------------------------------------------------------------------------------------
# colocalization inputs
# ---------------------------------------------------------------------------------------------
def resolve(path, staged_dir):
    if staged_dir:
        cand = os.path.join(staged_dir, os.path.basename(path))
        if os.path.exists(cand):
            return cand
    return path


def read_gwas(path, chroms):
    d = pd.read_csv(path, sep="\t")
    cols = {c.lower(): c for c in d.columns}
    d = d.rename(columns={v: k for k, v in cols.items()})
    if "variant_id" not in d.columns:
        need = {"chrom", "pos", "ref", "alt"}
        if not need <= set(d.columns):
            raise ValueError(f"{path}: needs variant_id or chrom, pos, ref, alt")
        ch = d["chrom"].astype(str)
        ch = np.where(ch.str.startswith("chr"), ch, "chr" + ch)
        d["variant_id"] = pd.Series(ch, index=d.index) + ":" + d["pos"].astype(str) + ":" + d["ref"] + ":" + d["alt"]
    d["chrom"] = d["variant_id"].str.split(":").str[0]
    d = d[d["chrom"].isin(chroms)]
    if "beta" not in d.columns or "se" not in d.columns:
        raise ValueError(f"{path}: needs beta and se columns (GRCh38, effect = alt allele)")
    return d


def cmd_coloc_prep(a):
    cfg = load_cfg(a.config)
    min_var = int(cfg_get(cfg, "central.coloc.min_variants", 50))
    cat = pd.read_csv(a.catalog, sep="\t", dtype=str)
    loci = pd.read_csv(a.loci, sep="\t")
    loci = loci[loci["status"] == "prepared"]
    os.makedirs(a.out_dir, exist_ok=True)
    tasks = []
    chroms = set(loci["chrom"])
    for t in cat.itertuples():
        g = read_gwas(resolve(t.path, a.staged_dir), chroms)
        g = g.drop_duplicates("variant_id").set_index("variant_id")
        sp = g.index.to_series().str.split(":", expand=True)
        swapped = (sp[0] + ":" + sp[1] + ":" + sp[3] + ":" + sp[2]).to_numpy()
        gs = g.copy()
        gs.index = swapped
        gs["beta"] = -gs["beta"]
        if "eaf" in gs.columns:
            gs["eaf"] = 1 - gs["eaf"]
        n_default = float(t.n) if getattr(t, "n", None) not in (None, "", np.nan) and pd.notna(t.n) else np.nan
        for L in loci.itertuples():
            ancs = str(L.ancestries).split(";")
            if t.ancestry not in ancs:
                continue
            z = np.load(os.path.join(a.loci_dir, f"{L.locus_id}.npz"))
            ai = ancs.index(t.ancestry)
            vids = list(z["variants"])
            x = pd.DataFrame({"variant_id": vids, "beta_x": z["betas"][ai], "se_x": z["ses"][ai], "eaf_x": z["eafs"][ai]})
            gg = pd.concat([g.reindex(vids).dropna(subset=["beta"]), gs.reindex(vids).dropna(subset=["beta"])])
            gg = gg[~gg.index.duplicated()]
            x = x.merge(gg[["beta", "se"] + (["n"] if "n" in gg.columns else [])].rename(
                columns={"beta": "beta_g", "se": "se_g", "n": "n_g"}), left_on="variant_id", right_index=True)
            if len(x) < min_var:
                tasks.append({"task_id": f"{L.locus_id}__{t.ancestry}__{sanitize(t.trait_id)}", "locus_id": L.locus_id,
                              "status": "too_few_shared_variants", "nsnps": len(x)})
                continue
            if "n_g" not in x.columns:
                x["n_g"] = n_default
            x["pos"] = x["variant_id"].str.split(":").str[1].astype(int)
            idx = pd.Index(vids).get_indexer(x["variant_id"])
            R = z["lds"][ai][np.ix_(idx, idx)].astype("<f8")
            tid = f"{L.locus_id}__{t.ancestry}__{sanitize(t.trait_id)}"
            x.to_csv(os.path.join(a.out_dir, f"{tid}.tsv"), sep="\t", index=False)
            R.tofile(os.path.join(a.out_dir, f"{tid}.ld.bin"))
            tasks.append({"task_id": tid, "locus_id": L.locus_id, "arm": L.arm, "modality": L.modality,
                          "phenotype_id": L.phenotype_id, "ancestry": t.ancestry, "trait_id": t.trait_id,
                          "trait_type": getattr(t, "trait_type", "quantitative") or "quantitative",
                          "prop_cases": getattr(t, "prop_cases", ""), "n_xqtl": float(z["ns"][ai]),
                          "nsnps": len(x), "status": "prepared"})
    pd.DataFrame(tasks, columns=["task_id", "locus_id", "arm", "modality", "phenotype_id", "ancestry", "trait_id",
                                 "trait_type", "prop_cases", "n_xqtl", "nsnps", "status"]).to_csv(
        os.path.join(a.out_dir, "coloc_tasks.tsv"), sep="\t", index=False)


# ---------------------------------------------------------------------------------------------
# G x C
# ---------------------------------------------------------------------------------------------
def cmd_gxc(a):
    cfg = load_cfg(a.config)
    fdr = float(cfg_get(cfg, "central.fdr", 0.05))
    min_k = int(cfg_get(cfg, "central.gxc_min_contributions", 2))
    con = contributions(a.contributions)
    hm = harm_map(a.harmonization)
    parts = []
    for r in con.itertuples():
        d = read_any(r.prefix + "_G_gxc")
        if d is None or d.empty:
            continue
        d = d[(d["status"] == "ok") & d["p_GxC"].notna()].copy()
        cols = ["modality", "phenotype_id", "variant_id", "cell_type", "source", "n", "beta_GxC", "se_GxC", "p_GxC",
                "maternal_sensitivity_beta_GxC", "maternal_sensitivity_se_GxC"]
        d = d[[c for c in cols if c in d.columns]]
        d = apply_harm(d, hm, r.contrib, beta_cols=("beta_GxC", "maternal_sensitivity_beta_GxC"))
        d["contrib"], d["ancestry"] = r.contrib, r.ancestry
        parts.append(d)
    os.makedirs(a.out_dir, exist_ok=True)
    if not parts:
        write(pd.DataFrame(), os.path.join(a.out_dir, "gxc_meta.tsv.gz"))
        write(pd.DataFrame(), os.path.join(a.out_dir, "gxc_annotation.tsv.gz"))
        return
    df = pd.concat(parts, ignore_index=True)
    keys = ["modality", "phenotype_id", "variant_id", "cell_type"]
    out = []
    for scope, sub in [(anc, s) for anc, s in df.groupby("ancestry")] + [("ALL", df)]:
        sub = sub.rename(columns={"beta_GxC": "beta", "se_GxC": "se"})
        m = meta_frame(sub, keys)
        g = sub.groupby(keys, sort=False)
        m["p_acat"] = g["p_GxC"].apply(lambda x: acat(x.to_numpy())).to_numpy()
        m["sources"] = g["source"].apply(lambda x: ";".join(sorted(set(map(str, x))))).to_numpy() \
            if "source" in sub.columns else ""
        m["method"] = np.where(m["k"] >= min_k, "ivw", "single")
        m["p_meta"] = np.where(m["k"] >= min_k, m["p_fe"], m["p_acat"])
        if "maternal_sensitivity_beta_GxC" in sub.columns:
            ms = sub.dropna(subset=["maternal_sensitivity_beta_GxC", "maternal_sensitivity_se_GxC"])
            ms = ms[ms["maternal_sensitivity_se_GxC"] > 0].rename(
                columns={"beta": "_b", "se": "_s", "maternal_sensitivity_beta_GxC": "beta",
                         "maternal_sensitivity_se_GxC": "se"})
            if len(ms):
                mm = meta_frame(ms, keys)[keys + ["k", "beta_fe", "se_fe", "p_fe"]].rename(columns={
                    "k": "maternal_k", "beta_fe": "maternal_beta_fe", "se_fe": "maternal_se_fe", "p_fe": "maternal_p_fe"})
                m = m.merge(mm, on=keys, how="left")
        m.insert(0, "scope", scope)
        m["q_meta"] = np.nan
        for _, idx in m.groupby("modality").groups.items():
            m.loc[idx, "q_meta"] = bh(m.loc[idx, "p_meta"].to_numpy())
        out.append(m)
    res = pd.concat(out, ignore_index=True)
    write(res, os.path.join(a.out_dir, "gxc_meta.tsv.gz"))
    ann = []
    for key, sub in res.groupby(["scope", "modality", "phenotype_id", "variant_id"], sort=True):
        s = sub.sort_values(["q_meta", "p_meta"])
        b = s.iloc[0]
        ann.append({"scope": key[0], "modality": key[1], "phenotype_id": key[2], "variant_id": key[3],
                    "primary_cell_type": b["cell_type"] if b["q_meta"] <= fdr else "unassigned",
                    "best_cell_type": b["cell_type"], "primary_method": b["method"], "primary_k": int(b["k"]),
                    "primary_beta_fe": b["beta_fe"], "primary_p_meta": b["p_meta"], "primary_q_meta": b["q_meta"],
                    "primary_i2": b["i2"], "maternal_p_fe": b.get("maternal_p_fe", np.nan),
                    "n_celltypes": len(s)})
    write(pd.DataFrame(ann), os.path.join(a.out_dir, "gxc_annotation.tsv.gz"))


# ---------------------------------------------------------------------------------------------
# G x E
# ---------------------------------------------------------------------------------------------
def exposure_scale(con):
    """Per exposure: pooled raw SD over all contributions' strata and each contribution's site scale SD.
    Rescaling applies to continuous, untransformed exposures only (b_gi is per site-pooled SD there)."""
    rows = []
    for r in con.itertuples():
        p = r.prefix + "_H_exposures.tsv"
        if os.path.exists(p):
            e = pd.read_csv(p, sep="\t")
            e["contrib"] = r.contrib
            rows.append(e)
    if not rows:
        return {}, pd.DataFrame()
    e = pd.concat(rows, ignore_index=True)
    strata = e[(e["ancestry"] != "pooled") & (e["status"] == "available")]
    out = {}
    for ex, s in strata.groupby("exposure"):
        typ = s["type"].iloc[0]
        tr = str(s["transform"].iloc[0])
        n, m, sd = s["n"].astype(float), s["mean_raw"].astype(float), s["sd_raw"].astype(float)
        N = n.sum()
        mbar = (n * m).sum() / N
        pooled_sd = float(np.sqrt((n * (sd ** 2 + (m - mbar) ** 2)).sum() / N))
        site_sd = s.set_index("contrib")["scale_sd"].astype(float).to_dict()
        rescale = typ == "continuous" and tr in ("none", "nan", "")
        out[ex] = {"pooled_sd": pooled_sd, "site_sd": site_sd, "rescale": rescale, "type": typ, "transform": tr,
                   "n": int(N), "pooled_mean": float(mbar)}
    tab = pd.DataFrame([{"exposure": k, **{kk: vv for kk, vv in v.items() if kk != "site_sd"}} for k, v in out.items()])
    return out, tab


def rescale(d, contrib, ex, scales, cols=("b_gi", "b_gi_se")):
    s = scales.get(ex)
    if not s or not s["rescale"]:
        return d, 1.0
    f = s["pooled_sd"] / s["site_sd"][contrib]
    for c in cols:
        d[c] = d[c] * f
    return d, f


def cmd_gxe_tier1(a):
    cfg = load_cfg(a.config)
    fdr = float(cfg_get(cfg, "central.fdr", 0.05))
    con = contributions(a.contributions)
    hm = harm_map(a.harmonization)
    scales, tab = exposure_scale(con)
    os.makedirs(a.out_dir, exist_ok=True)
    write(tab, os.path.join(a.out_dir, "gxe_exposure_scales.tsv"))
    parts = []
    for r in con.itertuples():
        d = read_any(r.prefix + "_H_gxe_tier1")
        if d is None or d.empty:
            continue
        d = d[d["pval_beta"].notna()]
        d = d[["modality", "exposure", "phenotype_id", "variant_id", "b_gi", "b_gi_se", "pval_beta", "qval", "n"]].copy()
        d = apply_harm(d, hm, r.contrib, beta_cols=("b_gi",))
        for ex in d["exposure"].unique():
            sel = d["exposure"] == ex
            sub, _ = rescale(d.loc[sel].copy(), r.contrib, ex, scales)
            d.loc[sel, ["b_gi", "b_gi_se"]] = sub[["b_gi", "b_gi_se"]].to_numpy()
        d["contrib"], d["study"], d["ancestry"] = r.contrib, r.study, r.ancestry
        parts.append(d)
    if not parts:
        write(pd.DataFrame(), os.path.join(a.out_dir, "gxe_tier1_meta.tsv.gz"))
        return
    df = pd.concat(parts, ignore_index=True)
    rows = []
    for (mod, ex, pid), sub in df.groupby(["modality", "exposure", "phenotype_id"], sort=True):
        k = len(sub)
        variants = sub["variant_id"].unique()
        rec = {"modality": mod, "exposure": ex, "phenotype_id": pid, "n_contributions": k,
               "contributions": ";".join(f"{s}/{an}" for s, an in zip(sub["study"], sub["ancestry"])),
               "p_feature_acat": acat(sub["pval_beta"].to_numpy()),
               "significant_any": bool((sub["qval"] <= fdr).any()),
               "significant_all": bool(k >= 2 and (sub["qval"] <= fdr).all()),
               "lead_concordant": bool(k >= 2 and len(variants) == 1)}
        rec["scope"] = "single_study" if k == 1 else ("shared_concordant_lead" if rec["lead_concordant"]
                                                      else "shared_discordant_lead")
        if rec["lead_concordant"]:
            m = ivw(np.zeros(k, dtype=int), sub["b_gi"].to_numpy(), sub["b_gi_se"].to_numpy(), 1)
            rec.update(meta_variant_id=variants[0], beta_meta=m["beta_fe"][0], se_meta=m["se_fe"][0],
                       p_meta=m["p_fe"][0], q_het=m["q_het"][0], p_het=m["p_het"][0], i2=m["i2"][0])
        elif k == 1:
            rec.update(meta_variant_id=variants[0], beta_meta=sub["b_gi"].iloc[0], se_meta=sub["b_gi_se"].iloc[0])
        rows.append(rec)
    out = pd.DataFrame(rows)
    out["q_feature_acat"] = np.nan
    for _, idx in out.groupby(["modality", "exposure"]).groups.items():
        out.loc[idx, "q_feature_acat"] = bh(out.loc[idx, "p_feature_acat"].to_numpy())
    out = out.sort_values(["modality", "exposure", "p_feature_acat"])
    write(out, os.path.join(a.out_dir, "gxe_tier1_meta.tsv.gz"))
    log(out.groupby(["exposure", "scope"]).size().to_string())


def cmd_gxe_tier2(a):
    con = contributions(a.contributions)
    hm = harm_map(a.harmonization)
    scales, _ = exposure_scale(con)
    os.makedirs(a.out_dir, exist_ok=True)
    units = {}
    for r in con.itertuples():
        for p in glob.glob(os.path.join(os.path.dirname(r.prefix), "H_gxe_tier2", "*_H_gxe_tier2_*")):
            stem = re.sub(r"\.(parquet|tsv\.gz)$", "", p)
            units.setdefault(stem.split("_H_gxe_tier2_", 1)[1], []).append((r, stem))
    exposures = sorted(scales, key=len, reverse=True)
    tops = []
    for unit, items in sorted(units.items()):
        ex = next((e for e in exposures if unit.endswith("_" + e)), unit.rsplit("_", 1)[-1])
        mod = unit[: -(len(ex) + 1)]
        parts = []
        for r, stem in items:
            d = read_any(stem, columns=["phenotype_id", "variant_id", "n", "af", "b_gi", "b_gi_se", "site_lead_variant",
                                        "mdacc_lead_variant"])
            d = apply_harm(d, hm, r.contrib, beta_cols=("b_gi",), eaf_cols=("af",))
            d, _ = rescale(d, r.contrib, ex, scales)
            d = d[np.isfinite(d["b_gi_se"]) & (d["b_gi_se"] > 0)]
            d["contrib"], d["ancestry"] = r.contrib, r.ancestry
            parts.append(d)
        df = pd.concat(parts, ignore_index=True).rename(columns={"b_gi": "beta", "b_gi_se": "se", "af": "eaf"})
        m = meta_frame(df, ["phenotype_id", "variant_id"])
        g = df.groupby(["phenotype_id", "variant_id"], sort=False)
        m["site_lead_any"] = g["site_lead_variant"].any().to_numpy()
        m["mdacc_lead_any"] = g["mdacc_lead_variant"].any().to_numpy()
        m.insert(0, "modality", mod)
        m.insert(1, "exposure", ex)
        src = df.groupby("phenotype_id").agg(site=("site_lead_variant", "any"), mdacc=("mdacc_lead_variant", "any"))
        m["source"] = m["phenotype_id"].map(
            (src["site"].map({True: "site", False: ""}) + src["mdacc"].map({True: "+mdacc", False: ""}))
            .str.strip("+").replace("", "none"))
        m["top_variant"] = False
        top = m.sort_values("p_fe").drop_duplicates("phenotype_id").index
        m.loc[top, "top_variant"] = True
        m["q_top_bh"] = np.nan
        m.loc[top, "q_top_bh"] = bh(m.loc[top, "p_fe"].to_numpy())
        write(m, os.path.join(a.out_dir, "tier2", f"gxe_tier2_meta_{mod}_{ex}.parquet"))
        tops.append(m.loc[top])
    t = pd.concat(tops, ignore_index=True) if tops else pd.DataFrame()
    write(t, os.path.join(a.out_dir, "gxe_tier2_top.tsv.gz"))


# ---------------------------------------------------------------------------------------------
# TWAS inputs (OTTERS / mvOTTERS stage I; FUSION weights from SuShiE)
# ---------------------------------------------------------------------------------------------
def cmd_twas_inputs(a):
    cfg = load_cfg(a.config)
    which = cfg_get(cfg, "central.twas.features", "significant")
    con = contributions(a.contributions)
    feats = pd.read_csv(a.features, sep="\t")
    feats = feats[~feats["grouped"]]
    if which == "significant":
        sig = set(feats.loc[feats["significant"], "phenotype_id"])
    else:
        sig = set(feats["phenotype_id"])
    coords = []
    for r in con.itertuples():
        e = read_any(r.prefix + "_E_features", columns=["modality", "phenotype_id", "chr", "start", "end"])
        if e is not None:
            coords.append(e)
    coords = pd.concat(coords).drop_duplicates(["modality", "phenotype_id"])
    os.makedirs(a.out_dir, exist_ok=True)
    meta = pd.concat([pd.read_parquet(p) for p in sorted(glob.glob(os.path.join(a.meta_dir, "*.parquet")))
                      if os.path.getsize(p) > 0], ignore_index=True)
    meta = meta[meta["phenotype_id"].isin(sig) & (meta["scope"] != "ALL")]
    for anc, sub in meta.groupby("scope"):
        d = os.path.join(a.out_dir, "otters", anc)
        os.makedirs(d, exist_ok=True)
        v = sub["variant_id"].str.split(":", expand=True)
        sst = pd.DataFrame({"CHROM": v[0].str.replace("chr", "", regex=False), "POS": v[1].astype(int),
                            "A1": v[3], "A2": v[2], "Zscore": sub["z_fe"].to_numpy(),
                            "TargetID": sub["modality"] + ":" + sub["phenotype_id"], "N": sub["n"].astype(int)})
        sst = sst.sort_values(["CHROM", "POS"], key=lambda s: s.map(lambda x: int(x) if str(x).isdigit() else 99)
                              if s.name == "CHROM" else s)
        sst.to_csv(os.path.join(d, "eQTLSumStats.txt"), sep="\t", index=False)
        an = coords[(coords["modality"] + ":" + coords["phenotype_id"]).isin(set(sst["TargetID"]))]
        pd.DataFrame({"CHROM": an["chr"].str.replace("chr", "", regex=False), "GeneStart": an["start"],
                      "GeneEnd": an["end"], "TargetID": an["modality"] + ":" + an["phenotype_id"]}).to_csv(
            os.path.join(d, "exp_anno.txt"), sep="\t", index=False)
    # SuShiE weights -> OTTERS weight format per ancestry (input of xqtl_fusion_wgt.R)
    if a.weights and os.path.exists(a.weights) and os.path.getsize(a.weights) > 0:
        w = pd.read_csv(a.weights, sep="\t")
        loci = pd.read_csv(a.loci, sep="\t").set_index("locus_id")
        for anc, sub in w.groupby("ancestry"):
            c = "weight"
            v = sub["variant_id"].str.split(":", expand=True)
            tid = sub["locus_id"].map(loci["modality"]) + ":" + sub["locus_id"].map(loci["phenotype_id"])
            out = pd.DataFrame({"CHROM": v[0].str.replace("chr", "", regex=False), "POS": v[1].astype(int),
                                "A1": v[3], "A2": v[2], "TargetID": tid, "ES": sub[c]})
            d = os.path.join(a.out_dir, "fusion", anc)
            os.makedirs(d, exist_ok=True)
            out.to_csv(os.path.join(d, "sushie_weights.txt"), sep="\t", index=False)


# ---------------------------------------------------------------------------------------------
# summary and PIP enrichment
# ---------------------------------------------------------------------------------------------
def read_bed(path):
    b = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2], comment="#", names=["chrom", "start", "end"],
                    compression="infer")
    b["chrom"] = np.where(b["chrom"].astype(str).str.startswith("chr"), b["chrom"].astype(str),
                          "chr" + b["chrom"].astype(str))
    return b


def in_bed(chrom, pos, bed):
    hit = np.zeros(len(pos), dtype=bool)
    for c, sub in bed.groupby("chrom"):
        sel = np.where(chrom == c)[0]
        if not len(sel):
            continue
        s = sub.sort_values("start")
        st, en = s["start"].to_numpy(), s["end"].to_numpy()
        en_max = np.maximum.accumulate(en)
        p = pos[sel] - 1   # BED is 0-based half-open
        i = np.searchsorted(st, p, side="right") - 1
        ok = i >= 0
        hit[sel[ok]] = en_max[i[ok]] > p[ok]
    return hit


def enrichment(pips, annot, staged_dir, thr):
    """Variant-level enrichment of high-PIP variants (pip >= thr) per annotation class (atlas 34)."""
    v = pips.sort_values("pip", ascending=False).drop_duplicates("variant_id").copy()
    sp = v["variant_id"].str.split(":", expand=True)
    chrom, pos = sp[0].to_numpy(), sp[1].astype(int).to_numpy()
    v["high"] = v["pip"] >= thr
    classes = {}
    if "tss_distance" in v.columns:
        d = v["tss_distance"].abs()
        classes["tss_le_1kb"] = (d <= 1000).to_numpy()
        classes["tss_1kb_10kb"] = ((d > 1000) & (d <= 10000)).to_numpy()
        classes["tss_10kb_100kb"] = ((d > 10000) & (d <= 100000)).to_numpy()
        classes["tss_gt_100kb"] = (d > 100000).to_numpy()
    if annot:
        an = pd.read_csv(annot, sep="\t", dtype=str)
        for r in an.itertuples():
            classes[r.name] = in_bed(chrom, pos, read_bed(resolve(r.bed, staged_dir)))
    rows = []
    hi = v["high"].to_numpy()
    tot_pip = v["pip"].sum()
    for name, x in classes.items():
        tab = [[int((x & hi).sum()), int((~x & hi).sum())], [int((x & ~hi).sum()), int((~x & ~hi).sum())]]
        orr, p = stats.fisher_exact(tab) if min(sum(tab[0]), sum(tab[1])) > 0 else (np.nan, np.nan)
        frac = x.mean() if len(x) else np.nan
        pw = (v.loc[x, "pip"].sum() / tot_pip) / frac if frac and tot_pip else np.nan
        rows.append({"annotation": name, "n_variants": len(x), "n_in_class": int(x.sum()),
                     "n_high_pip": int(hi.sum()), "n_high_pip_in_class": int((x & hi).sum()),
                     "odds_ratio": orr, "p_fisher": p, "pip_weighted_enrichment": pw})
    out = pd.DataFrame(rows)
    if len(out):
        out["q_fisher"] = bh(out["p_fisher"].to_numpy())
    return out


def cmd_summary(a):
    cfg = load_cfg(a.config)
    thr = float(cfg_get(cfg, "central.pip_threshold", 0.9))
    rows = []

    def add(section, item, value):
        if isinstance(value, (int, np.integer)) or (isinstance(value, float) and float(value).is_integer()):
            value = str(int(value))
        rows.append({"section": section, "item": item, "value": str(value)})

    con = contributions(a.contributions)
    add("ingest", "contributions", len(con))
    add("ingest", "studies", con["study"].nunique())
    for anc, n in con.groupby("ancestry")["n_individuals"].sum().items():
        add("ingest", f"n_individuals_{anc}", int(n))
    if a.features and os.path.exists(a.features):
        f = pd.read_csv(a.features, sep="\t")
        for (arm, scope), s in f.groupby(["arm", "scope"]):
            add("features", f"{arm}|{scope}|tested", len(s))
            add("features", f"{arm}|{scope}|significant", int(s["significant"].sum()))
    pips = None
    if a.pips and os.path.exists(a.pips) and os.path.getsize(a.pips) > 0:
        pips = pd.read_csv(a.pips, sep="\t")
        if a.loci and os.path.exists(a.loci):
            # distance to the feature anchor (phenotype BED end = TSS for genes, cluster/site position otherwise)
            loci = pd.read_csv(a.loci, sep="\t")
            coords = [read_any(r.prefix + "_E_features", columns=["modality", "phenotype_id", "end"])
                      for r in con.itertuples()]
            coords = pd.concat([c for c in coords if c is not None]).drop_duplicates(["modality", "phenotype_id"])
            anchor = loci.merge(coords, on=["modality", "phenotype_id"], how="left").set_index("locus_id")["end_y"]
            vpos = pips["variant_id"].str.split(":").str[1].astype(int)
            pips["tss_distance"] = vpos - pips["locus_id"].map(anchor)
        cs = pd.read_csv(a.cs, sep="\t") if a.cs and os.path.exists(a.cs) and os.path.getsize(a.cs) > 0 else None
        add("finemap", "loci", pips["locus_id"].nunique())
        add("finemap", f"variants_pip_ge_{thr}", int((pips["pip"] >= thr).sum()))
        if cs is not None and len(cs):
            sizes = cs.groupby(["locus_id", "cs_index"]).size()
            add("finemap", "credible_sets", len(sizes))
            add("finemap", "credible_set_size_median", float(sizes.median()))
            add("finemap", "loci_with_cs", cs["locus_id"].nunique())
    if a.status and os.path.exists(a.status) and os.path.getsize(a.status) > 0:
        st = pd.read_csv(a.status, sep="\t")
        for s, n in st["status"].value_counts().items():
            add("finemap", f"status_{s}", int(n))
    pp_h4 = float(cfg_get(cfg, "central.coloc.pp_h4", 0.7))
    if a.coloc and os.path.exists(a.coloc) and os.path.getsize(a.coloc) > 0:
        c = pd.read_csv(a.coloc, sep="\t")
        add("coloc_susie", "tasks", len(c))
        for m, n in c["method"].value_counts().items():
            add("coloc_susie", f"method_{m}", int(n))
        add("coloc_susie", f"tasks_PP.H4_ge_{pp_h4}", int((c["PP.H4"] >= pp_h4).sum()))
    if a.colocboost and os.path.exists(a.colocboost) and os.path.getsize(a.colocboost) > 0:
        c = pd.read_csv(a.colocboost, sep="\t")
        add("colocboost", "tasks", c["task_id"].nunique())
        add("colocboost", "tasks_with_colocalized_set", int(c.loc[c["status"] == "ok", "task_id"].nunique()))
    if a.gxc and os.path.exists(a.gxc):
        try:
            g = pd.read_csv(a.gxc, sep="\t")
            for scope, s in g.groupby("scope"):
                add("gxc", f"{scope}|pairs", len(s))
                add("gxc", f"{scope}|primary_assigned", int((s["primary_cell_type"] != "unassigned").sum()))
        except pd.errors.EmptyDataError:
            pass
    if a.tier1 and os.path.exists(a.tier1):
        try:
            t = pd.read_csv(a.tier1, sep="\t")
            for (ex, sc), s in t.groupby(["exposure", "scope"]):
                add("gxe_tier1", f"{ex}|{sc}", len(s))
            add("gxe_tier1", "q_feature_acat_le_fdr", int((t["q_feature_acat"] <= float(cfg_get(cfg, "central.fdr", 0.05))).sum()))
        except pd.errors.EmptyDataError:
            pass
    os.makedirs(a.out_dir, exist_ok=True)
    pd.DataFrame(rows).to_csv(os.path.join(a.out_dir, "central_summary.tsv"), sep="\t", index=False)
    if pips is not None and len(pips):
        enrichment(pips, a.annotations, a.staged_dir, thr).to_csv(
            os.path.join(a.out_dir, "pip_enrichment.tsv"), sep="\t", index=False, na_rep="NA")


# ---------------------------------------------------------------------------------------------
def cmd_concat(a):
    """Concatenate shard tables (same header); empty or header-only shards are skipped."""
    parts = []
    for p in sorted(a.inputs):
        if os.path.getsize(p) > 0:
            try:
                parts.append(pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False))
            except pd.errors.EmptyDataError:
                continue
    out = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    out.to_csv(a.out, sep="\t", index=False)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("ingest")
    s.add_argument("--studies", required=True)
    s.add_argument("--config", required=True)
    s.add_argument("--staged-prefix", default="", help="studies staged as <prefix>1/<dir>, <prefix>2/<dir>, ...")
    s.add_argument("--out-dir", required=True)
    s = sub.add_parser("meta-nominal")
    for x in ("--contributions", "--harmonization", "--modality", "--chrom", "--out"):
        s.add_argument(x, required=True)
    s = sub.add_parser("features")
    for x in ("--contributions", "--harmonization", "--meta-dir", "--config", "--out"):
        s.add_argument(x, required=True)
    s = sub.add_parser("finemap-prep")
    for x in ("--contributions", "--harmonization", "--meta-dir", "--features", "--config", "--out-dir"):
        s.add_argument(x, required=True)
    s = sub.add_parser("coloc-prep")
    for x in ("--catalog", "--loci", "--loci-dir", "--config", "--out-dir"):
        s.add_argument(x, required=True)
    s.add_argument("--staged-dir", default="")
    for name in ("gxc", "gxe-tier1", "gxe-tier2"):
        s = sub.add_parser(name)
        for x in ("--contributions", "--harmonization", "--config", "--out-dir"):
            s.add_argument(x, required=True)
    s = sub.add_parser("twas-inputs")
    for x in ("--contributions", "--meta-dir", "--features", "--config", "--out-dir"):
        s.add_argument(x, required=True)
    s.add_argument("--weights")
    s.add_argument("--loci")
    s = sub.add_parser("summary")
    s.add_argument("--contributions", required=True)
    s.add_argument("--config", required=True)
    s.add_argument("--out-dir", required=True)
    for x in ("--features", "--loci", "--pips", "--cs", "--status", "--coloc", "--colocboost", "--gxc", "--tier1",
              "--annotations", "--staged-dir"):
        s.add_argument(x)
    s = sub.add_parser("concat")
    s.add_argument("--out", required=True)
    s.add_argument("inputs", nargs="+")
    a = ap.parse_args(argv)
    {"concat": cmd_concat, "ingest": cmd_ingest, "meta-nominal": cmd_meta_nominal, "features": cmd_features,
     "finemap-prep": cmd_finemap_prep, "coloc-prep": cmd_coloc_prep, "gxc": cmd_gxc, "gxe-tier1": cmd_gxe_tier1,
     "gxe-tier2": cmd_gxe_tier2, "twas-inputs": cmd_twas_inputs, "summary": cmd_summary}[a.cmd](a)


if __name__ == "__main__":
    main()
