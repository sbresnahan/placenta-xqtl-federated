#!/usr/bin/env python3
"""xqtl_gxe.py -- SNP x exposure (G x E) interaction scans at the site .

  xqtl_gxe.py exposures   --metadata M --config gxe.json --samples EUR=EUR.samples.tsv ... [--out-dir .]
  xqtl_gxe.py scan        --tier 1|2 --anc EUR --modality expression --exposure GA --chrom chr22
                          --bed EUR_expression.bed.gz --cov EUR_expression.covariates.tsv --pfile EUR
                          --exposures gxe_exposures.tsv --summary EUR.gxe_exposures.tsv
                          [--phenotypes EUR_expression.gx_phenotypes.tsv] [options] --out PREFIX
  xqtl_gxe.py merge       --tier 1|2 --anc EUR --modality expression --exposure GA --shards P.parquet ...
                          [--phenotypes ...] [--fdr 0.05] --out EUR_expression_GA
  xqtl_gxe.py sensitivity --anc EUR --modality expression --exposure GA --tier1 T1.parquet --bed --cov
                          --pfile --exposures --metadata --samples --config gxe.json --out EUR_expression_GA

exposures: per enabled exposure in cohort.yml gxe.exposures, the metadata
  column is averaged over technical replicates to genotype IDs, range-filtered on the raw scale
  (min/max), transformed (none | log | inverse_normal), then continuous exposures are z-scored across
  all samples of the mapped strata, so per-stratum b_GxE share one scale.
  Binary exposures (0/1) are not scaled. Writes gxe_exposures.tsv (per individual; never exported)
  and <ANC>.gxe_exposures.tsv per stratum: n, raw-unit mean/SD, the pooled scale used for z, cases,
  controls, prevalence, and status (available | exposure_unavailable | insufficient_n |
  insufficient_cases), plus a 'pooled' row.

scan (cis permutation / nominal, gxe_core math, mirrors tensorQTL calculate_interaction_nominal):
  Y = G + E + G x E + covariates over the +-window cis window, one chromosome per call.
  Y: the stratum's normalized phenotype, autosomes only, z-scored within the stratum; rows with SD
  < 1e-8 dropped. Covariates: the modality's final mapping covariates; a row named like the exposure
  is dropped. Genotypes: hard calls by default, --dosages for dosages; missing values are
  mean-imputed; variants need MAF >= --maf in the analysed samples and MAF >= --maf-interaction in
  both halves of the stable-sorted exposure.
  Tier 1: every phenotype; adaptive phenotype permutations (--blocks, stop when running p > --stop-p),
  per-phenotype seed = --seed + row index in the full BED (independent of sharding); Beta-
  approximated pval_beta. Tier 2: every cis variant of the tier-2 phenotypes (signal sets), with
  flags for the site-lead and MDACC-lead variants (--tier2-variants pairs keeps the flagged rows
  only); untested phenotypes / lead variants listed with a reason.
merge: tier 1 -> Storey q on pval_beta within stratum x modality x exposure (BH fallback recorded);
  tier 2 -> per-phenotype top variant and BH q within the stratum (summary only; the cross-ancestry
  and cross-study IVW runs centrally). Each merge writes a status row.
sensitivity: tier-1 hits with q <= --fdr are refitted with Y ~ G + E + G:E + sum_j G:C_j +
  covariates, C_j = enabled gxe.sensitivity covariates (new columns also as main effects; ancestry_pcs
  uses the PC rows of the covariate table as interactions only). Status OFF (none enabled) |
  BLOCKED (metadata columns missing) | DONE.
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

T1_COLS = ["phenotype_id", "num_var", "variant_id", "start_distance", "end_distance", "af", "ma_samples",
           "ma_count", "maf_lower", "maf_upper", "n", "b_g", "b_g_se", "b_e", "b_e_se", "b_gi", "b_gi_se",
           "tstat_gi", "pval_g", "pval_e", "pval_nominal", "dof", "nperm", "pval_perm", "pval_beta",
           "beta_shape1", "beta_shape2", "true_df", "pval_true_df"]
T2_COLS = ["phenotype_id", "variant_id", "tss_distance", "af", "n", "b_gi", "b_gi_se", "tstat_gi",
           "pval_nominal", "dof", "site_lead_variant", "mdacc_lead_variant"]
UNT_COLS = ["phenotype_id", "variant_id", "untested_reason"]
STATUS_COLS = ["ancestry", "modality", "exposure", "component", "status", "reason", "n", "n_covariates",
               "dof", "n_tested", "n_sig"]


def log(msg):
    print(msg, flush=True)


# ------------------------------------------------------------------------------------------ exposures
def inverse_normal(v):
    from scipy import stats
    out = v.copy()
    ok = v.notna()
    r = v[ok].rank(method="average")
    out[ok] = stats.norm.ppf((r - 0.5) / ok.sum())
    return out


def per_individual(meta, samples_by_anc, col):
    """Metadata column averaged over replicate RNA-seq samples -> genotype ID (all strata)."""
    rows = []
    for anc, smp in samples_by_anc.items():
        m = meta[meta["rnaseq_id"].isin(smp["rnaseq_id"])].copy()
        m["gid"] = m["rnaseq_id"].map(dict(zip(smp["rnaseq_id"], smp["genotype_id"])))
        v = pd.to_numeric(m[col], errors="coerce").groupby(m["gid"]).mean()
        rows.append(pd.DataFrame({"genotype_id": v.index, "ancestry": anc, "value": v.values}))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=["genotype_id", "ancestry", "value"])


def build_exposures(meta, samples_by_anc, cfg):
    exps = {k: v for k, v in (cfg.get("exposures") or {}).items() if v.get("enabled")}
    min_n, min_cases = int(cfg.get("min_n", 50)), int(cfg.get("min_cases", 20))
    ids = pd.concat([s[["genotype_id"]].assign(ancestry=a) for a, s in samples_by_anc.items()]).drop_duplicates("genotype_id")
    table = pd.DataFrame(index=pd.Index(ids["genotype_id"], name="genotype_id"))
    summ = []
    for name, spec in exps.items():
        col = spec.get("column") or name
        typ = spec.get("type", "continuous")
        tr = spec.get("transform", "none")
        lo, hi = spec.get("min"), spec.get("max")
        base = {"exposure": name, "column": col, "type": typ, "transform": tr, "min": lo, "max": hi}
        if col not in meta.columns:
            for anc in list(samples_by_anc) + ["pooled"]:
                summ.append({**base, "ancestry": anc, "n": 0, "status": "exposure_unavailable",
                             "reason": f"column {col} not in metadata"})
            continue
        ind = per_individual(meta, samples_by_anc, col).drop_duplicates("genotype_id").set_index("genotype_id")
        raw = ind["value"].copy()
        if typ == "binary":
            raw = raw.where(raw.isin([0.0, 1.0]))      # replicate disagreement -> missing
        if lo is not None:
            raw = raw.where(raw >= float(lo))
        if hi is not None:
            raw = raw.where(raw <= float(hi))
        val = raw.copy()
        if typ == "continuous":
            if tr == "log":
                val = np.log(val.where(val > 0))
            elif tr == "inverse_normal":
                val = inverse_normal(val)
            mu, sd = float(val.mean()), float(val.std(ddof=0))
            val = (val - mu) / sd if np.isfinite(sd) and sd > 0 else val * np.nan
        else:
            mu, sd = np.nan, np.nan
        table[name] = val.reindex(table.index)
        groups = [(a, ind.index[ind["ancestry"] == a]) for a in samples_by_anc] + [("pooled", ind.index)]
        for anc, idx in groups:
            r = raw.reindex(idx)
            v = val.reindex(idx)
            n = int(v.notna().sum())
            rec = {**base, "ancestry": anc, "n": n, "n_missing": int(v.isna().sum()),
                   "mean_raw": float(r[v.notna()].mean()) if n else np.nan,
                   "sd_raw": float(r[v.notna()].std(ddof=0)) if n else np.nan,
                   "scale_mean": mu, "scale_sd": sd}
            status, reason = "available", ""
            if typ == "binary":
                cases = int((r[v.notna()] == 1).sum())
                rec.update({"n_cases": cases, "n_controls": n - cases, "prevalence": cases / n if n else np.nan})
                if anc != "pooled" and min(cases, n - cases) < min_cases:
                    status, reason = "insufficient_cases", f"cases {cases}, controls {n - cases} < gxe.min_cases {min_cases}"
            if n == 0:
                status, reason = "exposure_unavailable", "no non-missing values"
            elif anc != "pooled" and n < min_n and status == "available":
                status, reason = "insufficient_n", f"n {n} < gxe.min_n {min_n}"
            rec.update({"status": status, "reason": reason})
            summ.append(rec)
    cols = ["ancestry", "exposure", "column", "type", "transform", "min", "max", "n", "n_missing", "mean_raw",
            "sd_raw", "scale_mean", "scale_sd", "n_cases", "n_controls", "prevalence", "status", "reason"]
    return table, pd.DataFrame(summ).reindex(columns=cols)


def cmd_exposures(a):
    meta = pd.read_csv(a.metadata, sep="\t", dtype={"rnaseq_id": str, "genotype_id": str})
    with open(a.config) as fh:
        cfg = json.load(fh)
    samples = {}
    for s in a.samples:
        anc, path = s.split("=", 1)
        samples[anc] = pd.read_csv(path, sep="\t", dtype=str)
    table, summ = build_exposures(meta, samples, cfg)
    table.to_csv(os.path.join(a.out_dir, "gxe_exposures.tsv"), sep="\t", float_format="%.8g")
    for anc in samples:
        summ[summ["ancestry"].isin([anc, "pooled"])].to_csv(
            os.path.join(a.out_dir, f"{anc}.gxe_exposures.tsv"), sep="\t", index=False)
    log(f"gxe_exposures.tsv: {table.shape[1]} exposure(s) x {table.shape[0]} individuals")
    for r in summ.itertuples():
        log(f"  {r.ancestry} {r.exposure}: n={r.n} {r.status} {r.reason or ''}")


# ------------------------------------------------------------------------------------------ inputs
def interaction_maf_mask(G, interaction, maf_threshold_interaction=0.05):
    """tensorqtl.core.filter_maf_interaction for one interaction term."""
    if not (0.0 <= maf_threshold_interaction <= 0.5):
        raise ValueError(f"maf_threshold_interaction must be between 0 and 0.5; got {maf_threshold_interaction}")
    if maf_threshold_interaction <= 0:
        return np.ones(G.shape[0], dtype=bool)
    interaction = np.asarray(interaction, dtype=float)
    if interaction.ndim != 1 or interaction.shape[0] != G.shape[1]:
        raise ValueError("interaction must be a 1D vector with one value per genotype sample")
    if not np.isfinite(interaction).all():
        raise ValueError("interaction contains non-finite values after sample alignment")
    if G.shape[1] < 2:
        return np.zeros(G.shape[0], dtype=bool)
    lo, hi = halves(interaction)

    def _maf(ix):
        af = G[:, ix].sum(1) / (2.0 * len(ix))
        return np.where(af > 0.5, 1.0 - af, af)
    eps = 1e-7
    return (_maf(lo) >= maf_threshold_interaction - eps) & (_maf(hi) >= maf_threshold_interaction - eps)


def halves(interaction):
    order = np.argsort(np.asarray(interaction, float), kind="mergesort")
    n = len(order)
    return order[:n // 2], order[n // 2:]


def zscore_rows(x):
    mu = x.mean(1)
    sd = x.std(1, ddof=0)
    sd = np.where(sd < 1e-8, np.nan, sd)
    return (x - mu[:, None]) / sd[:, None]


def chrom_key(c):
    return str(c).lower().replace("chr", "")


def read_bed(path):
    bed = pd.read_csv(path, sep="\t", dtype={0: str, 3: str})
    bed.columns = ["chr", "start", "end", "phenotype_id"] + [str(c) for c in bed.columns[4:]]
    bed["row"] = np.arange(len(bed))
    return bed


def load_genotypes(pfile, chrom, samples, dosages=False):
    """Variants on one chromosome x samples (float32, mean-imputed); ids, positions."""
    from tensorqtl import pgen as tq
    pvar = tq.read_pvar(pfile + ".pvar")
    psam = tq.read_psam(pfile + ".psam")
    iids = [str(i) for i in psam.index]
    sel = np.where(pvar["chrom"].astype(str).map(chrom_key) == chrom_key(chrom))[0]
    if len(sel) == 0:
        return np.empty((0, len(samples)), np.float32), np.array([], dtype=object), np.array([], dtype=np.int64), set()
    reader = tq.read_dosages_list if dosages else tq.read_list
    g = reader(pfile + ".pgen", np.asarray(sel, dtype=np.uint32)).astype(np.float32)
    col = {s: i for i, s in enumerate(iids)}
    miss = [s for s in samples if s not in col]
    if miss:
        sys.exit(f"xqtl_gxe.py: {len(miss)} analysis samples missing from {pfile}.psam")
    g = g[:, [col[s] for s in samples]]
    if dosages:
        g[(g < 0) | (g > 2)] = np.nan
    else:
        g[g == -9] = np.nan
    if np.isnan(g).any():
        mu = np.nanmean(g, axis=1)
        ix = np.where(np.isnan(g) & np.isfinite(mu)[:, None])
        g[ix] = np.take(mu, ix[0])
    ids = pvar["id"].to_numpy()[sel].astype(str)
    pos = pvar["pos"].to_numpy()[sel].astype(np.int64)
    ok = np.isfinite(g).all(1)
    return g[ok], ids[ok], pos[ok], set(ids)


def status_row(a, component, status, reason="", **kw):
    return {"ancestry": a.anc, "modality": a.modality, "exposure": a.exposure, "component": component,
            "status": status, "reason": reason, **kw}


def write_scan(a, rows, cols, status, untested=None):
    pd.DataFrame(rows, columns=cols).to_parquet(a.out + ".parquet", index=False)
    pd.DataFrame(untested or [], columns=UNT_COLS).to_csv(a.out + ".untested.tsv", sep="\t", index=False)
    pd.DataFrame([status], columns=STATUS_COLS).to_csv(a.out + ".status.tsv", sep="\t", index=False)
    log(f"{a.out}: {len(rows)} rows; {status['status']} {status.get('reason', '')}")


def cmd_scan(a):
    import torch
    import xqtl_gxe_core as core
    torch.set_num_threads(max(1, a.threads))
    comp = f"tier{a.tier}"
    tier2 = a.tier == 2
    cols = T2_COLS if tier2 else T1_COLS

    t2 = None
    if tier2:
        t2 = pd.read_csv(a.phenotypes, sep="\t", dtype=str, keep_default_na=False)
        t2 = t2[t2["untested_reason"] == ""]

    summ = pd.read_csv(a.summary, sep="\t")
    srow = summ[(summ["ancestry"] == a.anc) & (summ["exposure"] == a.exposure)]
    bed = read_bed(a.bed)
    onchr = bed[bed["chr"].map(chrom_key) == chrom_key(a.chrom)]
    auto = chrom_key(a.chrom).isdigit() and 1 <= int(chrom_key(a.chrom)) <= 22

    def t2_untested(reason):
        if t2 is None:
            return []
        ids = set(onchr["phenotype_id"])
        return [{"phenotype_id": p, "variant_id": "", "untested_reason": reason} for p in t2["phenotype_id"] if p in ids]

    if len(srow) == 0 or srow["status"].iloc[0] != "available":
        reason = "exposure_unavailable" if len(srow) == 0 else str(srow["status"].iloc[0])
        return write_scan(a, [], cols, status_row(a, comp, "SKIPPED", reason), t2_untested(reason))
    if not auto or len(onchr) == 0:
        return write_scan(a, [], cols, status_row(a, comp, "DONE", "no autosomal phenotypes on " + a.chrom,
                                                  n_tested=0), [])

    exp = pd.read_csv(a.exposures, sep="\t", index_col=0, dtype={0: str})[a.exposure].dropna()
    cov = pd.read_csv(a.cov, sep="\t", index_col=0)
    cov.columns = cov.columns.astype(str)
    if a.exposure in cov.index:
        cov = cov.drop(index=a.exposure)
    bed_samples = [c for c in bed.columns[4:-1]]
    samples = [s for s in bed_samples if s in set(cov.columns) and s in set(exp.index)]
    n = len(samples)
    if n < a.min_n:
        return write_scan(a, [], cols, status_row(a, comp, "SKIPPED", "insufficient_n", n=n), t2_untested("insufficient_n"))
    # phenotypes z-scored within the stratum, then aligned to the analysis samples
    vals = zscore_rows(onchr[bed_samples].to_numpy(dtype=np.float64))
    zdf = pd.DataFrame(vals, index=onchr["phenotype_id"].values, columns=bed_samples)
    good = np.isfinite(vals).all(1)
    zero_var = set(onchr["phenotype_id"].values[~good])
    ph = zdf[good][samples]
    rowidx = dict(zip(onchr["phenotype_id"], onchr["row"]))
    pos_s = onchr.set_index("phenotype_id")[["start", "end"]]

    covm = cov[samples].T.to_numpy(dtype=np.float64)
    if not np.isfinite(covm).all():
        sys.exit("xqtl_gxe.py: non-finite covariate values after sample alignment")
    e = exp[samples].to_numpy(dtype=np.float64)
    G, vids, vpos, present = load_genotypes(a.pfile, a.chrom, samples, a.dosages)
    order = np.argsort(vpos, kind="mergesort")
    G, vids, vpos = G[order], vids[order], vpos[order]
    af = G.sum(1) / (2.0 * n) if len(G) else np.array([])
    maf = np.where(af > 0.5, 1 - af, af)
    mono = (G == G[:, [0]]).all(1) if len(G) else np.array([], bool)
    keep = (maf >= a.maf) & ~mono
    if a.maf_interaction > 0 and len(G):
        keep &= interaction_maf_mask(G, e, a.maf_interaction)
    G, vids, vpos = G[keep], vids[keep], vpos[keep]
    lo_ix, hi_ix = halves(e)

    device = core.get_device()
    k = covm.shape[1]
    Q_t, dof_base = core.make_residualizer(covm, device) if k else (None, n - 2)
    dof = dof_base - 2
    base = {"n": n, "n_covariates": k, "dof": dof}
    if dof < 10:
        return write_scan(a, [], cols, status_row(a, comp, "SKIPPED", "insufficient_df", **base), t2_untested("insufficient_df"))
    G_t = torch.tensor(G, dtype=torch.float32, device=device)
    e_t = torch.tensor(e, dtype=torch.float32, device=device)
    g0_t, e0_t, ge0_t = core.prepare_interaction_terms(G_t, e_t, Q_t)
    af_t, ma_s_t, ma_c_t = core.allele_stats(G_t)
    af, ma_s, ma_c = af_t.cpu().numpy(), ma_s_t.cpu().numpy(), ma_c_t.cpu().numpy()

    rows, untested = [], []
    if tier2:
        lead_site = dict(zip(t2["phenotype_id"], t2["site_lead_variant"]))
        lead_md = {p: set(filter(None, v.split(","))) for p, v in zip(t2["phenotype_id"], t2["mdacc_lead_variants"])}
        targets = [p for p in t2["phenotype_id"] if p in set(onchr["phenotype_id"])]
        for p in targets:
            if p in zero_var:
                untested.append({"phenotype_id": p, "variant_id": "", "untested_reason": "phenotype_absent"})
        targets = [p for p in targets if p in ph.index]
    else:
        targets = list(ph.index)
    blocks = tuple(int(b) for b in a.blocks)
    for i, pid in enumerate(targets, 1):
        p_raw_t = torch.tensor(ph.loc[pid].to_numpy(np.float32), device=device)
        p0_t = core.prepare_phenotype(p_raw_t, Q_t)
        start, end = int(pos_s.loc[pid, "start"]), int(pos_s.loc[pid, "end"])
        lo = np.searchsorted(vpos, start - a.window, side="left")
        hi = np.searchsorted(vpos, end + a.window, side="right")
        if tier2:
            leads = {lead_site.get(pid, "")} | lead_md.get(pid, set())
            win = set(vids[lo:hi])
            for v in sorted(x for x in leads if x):
                if v not in present:
                    untested.append({"phenotype_id": pid, "variant_id": v, "untested_reason": "variant_absent"})
                elif v not in set(vids):
                    untested.append({"phenotype_id": pid, "variant_id": v, "untested_reason": "interaction_maf"})
                elif v not in win:
                    untested.append({"phenotype_id": pid, "variant_id": v, "untested_reason": "outside_cis_window"})
        if hi - lo < 1:
            if tier2:
                untested.append({"phenotype_id": pid, "variant_id": "", "untested_reason": "no_cis_variants"})
            continue
        sl = slice(lo, hi)
        Xt_w, Xinv_w, ok_w = core.build_window_design(g0_t[sl], e0_t, ge0_t[sl])
        b, b_se, tstat = core.window_fit(Xt_w, Xinv_w, p0_t, dof)
        t_gi = tstat[:, 2].cpu().numpy()
        if tier2:
            pv = core.r2_to_pval(core.t_to_r2(t_gi, dof), dof)
            df = pd.DataFrame({"phenotype_id": pid, "variant_id": vids[sl], "tss_distance": (vpos[sl] - start).astype(int),
                               "af": af[sl], "n": n, "b_gi": b[:, 2].cpu().numpy(), "b_gi_se": b_se[:, 2].cpu().numpy(),
                               "tstat_gi": t_gi, "pval_nominal": pv, "dof": int(dof)})
            df = df[np.isfinite(df["pval_nominal"])]
            df["site_lead_variant"] = df["variant_id"] == lead_site.get(pid, "")
            df["mdacc_lead_variant"] = df["variant_id"].isin(lead_md.get(pid, set()))
            if a.tier2_variants == "pairs":
                df = df[df["site_lead_variant"] | df["mdacc_lead_variant"]]
            rows.extend(df.to_dict("records"))
            continue
        r2 = core.t_to_r2(t_gi, dof)
        r2 = np.where(np.isfinite(r2), r2, -1.0)
        ix = int(np.argmax(r2))
        if r2[ix] < 0:
            continue
        r2max = float(r2[ix])
        rng = np.random.default_rng(a.seed + int(rowidx[pid]))
        r2_perm, nperm, pval_perm = core.adaptive_permutation_scan(Xt_w, Xinv_w, p_raw_t, Q_t, dof, r2max, rng,
                                                                   blocks=blocks, stop_p=a.stop_p)
        beta = core.beta_approx_pval(r2_perm, r2max, dof)
        gv = lo + ix
        gl, gu = G[gv, lo_ix], G[gv, hi_ix]
        mlo, mup = gl.mean() / 2, gu.mean() / 2
        bb, bs, tt = b[ix].cpu().numpy(), b_se[ix].cpu().numpy(), tstat[ix].cpu().numpy()
        rows.append(dict(phenotype_id=pid, num_var=int(hi - lo), variant_id=vids[gv],
                         start_distance=int(vpos[gv] - start), end_distance=int(vpos[gv] - end),
                         af=float(af[gv]), ma_samples=int(ma_s[gv]), ma_count=float(ma_c[gv]),
                         maf_lower=float(min(mlo, 1 - mlo)), maf_upper=float(min(mup, 1 - mup)), n=n,
                         b_g=float(bb[0]), b_g_se=float(bs[0]), b_e=float(bb[1]), b_e_se=float(bs[1]),
                         b_gi=float(bb[2]), b_gi_se=float(bs[2]), tstat_gi=float(tt[2]),
                         pval_g=float(core.r2_to_pval(core.t_to_r2(tt[0], dof), dof)),
                         pval_e=float(core.r2_to_pval(core.t_to_r2(tt[1], dof), dof)),
                         pval_nominal=float(core.r2_to_pval(r2max, dof)), dof=int(dof), nperm=int(nperm),
                         pval_perm=float(pval_perm), **beta))
        if i % 1000 == 0:
            log(f"  {i}/{len(targets)} phenotypes")
    write_scan(a, rows, cols, status_row(a, comp, "DONE", "", n_tested=len(targets), **base), untested)


# ------------------------------------------------------------------------------------------ merge
def cmd_merge(a):
    comp = f"tier{a.tier}"
    stats = pd.concat([pd.read_csv(p.replace(".parquet", ".status.tsv"), sep="\t") for p in a.shards], ignore_index=True) \
        if a.shards else pd.DataFrame(columns=STATUS_COLS)
    parts = [pd.read_parquet(p) for p in a.shards]
    parts = [p for p in parts if len(p)]
    df = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=T2_COLS if a.tier == 2 else T1_COLS)
    done = stats[stats["status"] == "DONE"]
    if len(done):
        status, reason = "DONE", ""
    else:
        status = "SKIPPED"
        reason = ",".join(sorted(set(stats["reason"].dropna().astype(str)))) or "no shards"
    n = done["n"].dropna().max() if len(done) and done["n"].notna().any() else (stats["n"].dropna().max() if stats["n"].notna().any() else np.nan)
    rec = {"ancestry": a.anc, "modality": a.modality, "exposure": a.exposure, "component": comp, "status": status,
           "reason": reason, "n": n, "n_covariates": done["n_covariates"].dropna().max() if len(done) else np.nan,
           "dof": done["dof"].dropna().max() if len(done) else np.nan,
           "n_tested": int(done["n_tested"].fillna(0).sum()) if len(done) else 0}
    if a.tier == 1:
        if len(df):
            import xqtl_tensorqtl as tq
            df = tq.add_qvalues(df)
        else:
            df["qval"], df["qvalue_method"] = pd.Series(dtype=float), pd.Series(dtype=str)
        rec["n_sig"] = int((pd.to_numeric(df["qval"], errors="coerce") <= a.fdr).sum())
        df.to_parquet(a.out + ".gxe_tier1.parquet", index=False)
    else:
        unt = [pd.read_csv(p.replace(".parquet", ".untested.tsv"), sep="\t", dtype=str, keep_default_na=False) for p in a.shards]
        if a.phenotypes:
            t2 = pd.read_csv(a.phenotypes, sep="\t", dtype=str, keep_default_na=False)
            ab = t2[t2["untested_reason"] != ""]
            unt.append(pd.DataFrame({"phenotype_id": ab["phenotype_id"], "variant_id": "", "untested_reason": ab["untested_reason"]}))
            # tier-2 phenotypes missing from every shard (e.g. no shard for their chromosome)
            seen = set(df["phenotype_id"]) | set().union(*[set(u["phenotype_id"]) for u in unt]) if unt else set(df["phenotype_id"])
            rest = t2[(t2["untested_reason"] == "") & ~t2["phenotype_id"].isin(seen)]
            unt.append(pd.DataFrame({"phenotype_id": rest["phenotype_id"], "variant_id": "",
                                     "untested_reason": "exposure_unavailable" if status == "SKIPPED" else "phenotype_absent"}))
        unt = pd.concat(unt, ignore_index=True).reindex(columns=UNT_COLS).drop_duplicates() if unt else pd.DataFrame(columns=UNT_COLS)
        df["top_variant"] = False
        df["q_top_bh"] = np.nan
        if len(df):
            import xqtl_gxc as gx
            top = df.sort_values(["phenotype_id", "pval_nominal"]).drop_duplicates("phenotype_id").index
            df.loc[top, "top_variant"] = True
            df.loc[top, "q_top_bh"] = gx.bh(df.loc[top, "pval_nominal"].to_numpy(float))
        rec["n_sig"] = int((df["q_top_bh"] <= a.fdr).sum())
        df.to_parquet(a.out + ".gxe_tier2.parquet", index=False)
        unt.to_csv(a.out + ".gxe_tier2_untested.tsv", sep="\t", index=False)
    pd.DataFrame([rec], columns=STATUS_COLS).to_csv(a.out + f".gxe_tier{a.tier}.status.tsv", sep="\t", index=False)
    log(f"{a.out} tier {a.tier}: {len(df)} rows; {status} {reason}; n_sig={rec['n_sig']}")


# ------------------------------------------------------------------------------------------ sensitivity
def ols_b_se(y, X, idx):
    """OLS via lstsq; (b, se, p) for column idx."""
    from scipy import stats
    beta, _, rank, _ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = max(len(y) - rank, 1)
    sigma2 = float(resid @ resid) / dof
    se = np.sqrt(max(np.linalg.pinv(X.T @ X)[idx, idx] * sigma2, 1e-300))
    t = beta[idx] / se
    return float(beta[idx]), float(se), float(2 * stats.t.sf(abs(t), dof))


SENS_COLS = ["ancestry", "modality", "exposure", "phenotype_id", "variant_id", "covariates", "n", "b_int_base",
             "se_base", "p_base", "b_int_adj", "se_adj", "p_adj", "delta_b", "rel_change", "stable"]


def cmd_sensitivity(a):
    with open(a.config) as fh:
        cfg = json.load(fh)
    sens = cfg.get("sensitivity") or {}
    exp_col = ((cfg.get("exposures") or {}).get(a.exposure) or {}).get("column") or a.exposure
    enabled = {k: v for k, v in sens.items() if isinstance(v, dict) and v.get("enabled")}
    # an exposure is never its own sensitivity covariate (GA is used for non-GA exposures only)
    enabled = {k: v for k, v in enabled.items() if k == "ancestry_pcs" or (v.get("column") or k) != exp_col}
    st = {"ancestry": a.anc, "modality": a.modality, "exposure": a.exposure, "component": "sensitivity"}
    out = pd.DataFrame(columns=SENS_COLS)

    def finish(status, reason, n_tested=0, n_sig=np.nan, n=np.nan):
        out.to_csv(a.out + ".gxe_sensitivity.tsv", sep="\t", index=False)
        pd.DataFrame([{**st, "status": status, "reason": reason, "n": n, "n_tested": n_tested, "n_sig": n_sig}],
                     columns=STATUS_COLS).to_csv(a.out + ".gxe_sensitivity.status.tsv", sep="\t", index=False)
        log(f"{a.out} sensitivity: {status} {reason}")

    if not enabled:
        return finish("OFF", "no gxe.sensitivity covariates enabled")
    meta = pd.read_csv(a.metadata, sep="\t", dtype={"rnaseq_id": str, "genotype_id": str})
    cols = {k: (v.get("column") or k) for k, v in enabled.items() if k != "ancestry_pcs"}
    missing = sorted(c for c in cols.values() if c not in meta.columns)
    if missing:
        return finish("BLOCKED", "metadata columns missing: " + ",".join(missing))
    t1 = pd.read_parquet(a.tier1)
    hits = t1[pd.to_numeric(t1.get("qval", pd.Series(dtype=float)), errors="coerce") <= a.fdr] if len(t1) else t1
    if len(hits) == 0:
        return finish("DONE", f"no tier-1 hits at q <= {a.fdr}")
    smp = pd.read_csv(a.samples, sep="\t", dtype=str)
    exp = pd.read_csv(a.exposures, sep="\t", index_col=0, dtype={0: str})[a.exposure].dropna()
    cov = pd.read_csv(a.cov, sep="\t", index_col=0)
    cov.columns = cov.columns.astype(str)
    if a.exposure in cov.index:
        cov = cov.drop(index=a.exposure)
    bed = read_bed(a.bed)
    bed_samples = list(bed.columns[4:-1])
    samples = [s for s in bed_samples if s in set(cov.columns) and s in set(exp.index)]
    vals = zscore_rows(bed[bed_samples].to_numpy(dtype=np.float64))
    ph = pd.DataFrame(vals, index=bed["phenotype_id"].values, columns=bed_samples)[samples]
    extra = {}
    for k, c in cols.items():
        v = per_individual(meta, {a.anc: smp}, c).set_index("genotype_id")["value"].reindex(samples)
        extra[k] = ((v - v.mean()) / v.std(ddof=0)).fillna(0.0).to_numpy()
    pcs = [r for r in cov.index if str(r).startswith("PC") and str(r)[2:].isdigit()]
    rows = []
    for chrom, hh in hits.groupby(hits["variant_id"].str.split(":").str[0]):
        G, vids, _, _ = load_genotypes(a.pfile, chrom, samples, dosages=False)
        vix = {v: i for i, v in enumerate(vids)}
        for h in hh.itertuples():
            if h.phenotype_id not in ph.index or h.variant_id not in vix:
                continue
            y = ph.loc[h.phenotype_id].to_numpy(float)
            g = G[vix[h.variant_id]].astype(float)
            e = exp[samples].to_numpy(float)
            Z = cov[samples].T.to_numpy(float)
            mains, gx = [], []
            for k in enabled:
                if k == "ancestry_pcs":
                    C = cov.loc[pcs, samples].T.to_numpy(float)
                else:
                    C = extra[k].reshape(-1, 1)
                    mains.append(C)
                gx.append(g.reshape(-1, 1) * C)
            X0 = np.column_stack([np.ones(len(samples)), g, e, g * e, Z])
            X1 = np.column_stack([X0] + mains + gx)
            b0, se0, p0 = ols_b_se(y, X0, 3)
            b1, se1, p1 = ols_b_se(y, X1, 3)
            rows.append({**{k: st[k] for k in ("ancestry", "modality", "exposure")}, "phenotype_id": h.phenotype_id,
                         "variant_id": h.variant_id, "covariates": ",".join(enabled), "n": len(samples),
                         "b_int_base": b0, "se_base": se0, "p_base": p0, "b_int_adj": b1, "se_adj": se1, "p_adj": p1,
                         "delta_b": b1 - b0, "rel_change": abs(b1 - b0) / abs(b0) if b0 != 0 else np.nan,
                         "stable": bool(abs(b1 - b0) <= 0.2 * abs(b0) and ((p0 < 0.05) == (p1 < 0.05)))})
    out = pd.DataFrame(rows, columns=SENS_COLS)
    finish("DONE", f"{len(out)} refits; {int(out['stable'].sum()) if len(out) else 0} stable", len(out),
           int(out["stable"].sum()) if len(out) else 0, len(samples))


# ------------------------------------------------------------------------------------------ main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("exposures")
    e.add_argument("--metadata", required=True)
    e.add_argument("--config", required=True)
    e.add_argument("--samples", nargs="+", required=True)
    e.add_argument("--out-dir", default=".")
    s = sub.add_parser("scan")
    for k in ("--anc", "--modality", "--exposure", "--chrom", "--bed", "--cov", "--pfile", "--exposures", "--summary", "--out"):
        s.add_argument(k, required=True)
    s.add_argument("--tier", type=int, choices=[1, 2], required=True)
    s.add_argument("--phenotypes")
    s.add_argument("--window", type=int, default=1_000_000)
    s.add_argument("--maf", type=float, default=0.01)
    s.add_argument("--maf-interaction", type=float, default=0.05)
    s.add_argument("--min-n", type=int, default=50)
    s.add_argument("--blocks", type=int, nargs="+", default=[100, 400, 500, 9000])
    s.add_argument("--stop-p", type=float, default=0.10)
    s.add_argument("--seed", type=int, default=12345)
    s.add_argument("--dosages", action="store_true")
    s.add_argument("--tier2-variants", choices=["all_cis_union", "pairs"], default="all_cis_union",
                   help="tier 2: every cis variant (default) or only the site-lead / MDACC-lead variants")
    s.add_argument("--threads", type=int, default=1)
    m = sub.add_parser("merge")
    for k in ("--anc", "--modality", "--exposure", "--out"):
        m.add_argument(k, required=True)
    m.add_argument("--tier", type=int, choices=[1, 2], required=True)
    m.add_argument("--shards", nargs="*", default=[])
    m.add_argument("--phenotypes")
    m.add_argument("--fdr", type=float, default=0.05)
    n = sub.add_parser("sensitivity")
    for k in ("--anc", "--modality", "--exposure", "--tier1", "--bed", "--cov", "--pfile", "--exposures",
              "--metadata", "--samples", "--config", "--out"):
        n.add_argument(k, required=True)
    n.add_argument("--fdr", type=float, default=0.05)
    a = ap.parse_args(argv)
    if a.cmd == "scan" and a.tier == 2 and not a.phenotypes:
        ap.error("scan --tier 2 needs --phenotypes")
    {"exposures": cmd_exposures, "scan": cmd_scan, "merge": cmd_merge, "sensitivity": cmd_sensitivity}[a.cmd](a)


if __name__ == "__main__":
    main()
