#!/usr/bin/env python3
"""check_central.py -- acceptance checks on a --mode central run of the two test cohorts.

  check_central.py extract CENTRAL_DIR OUT_DIR [--pairs 30 --seed 1]
      writes, for the planted pairs and a seeded random sample of meta rows, the raw file A rows of every
      contribution (ivw_inputs.tsv, read straight from the study exports, before harmonization; the test
      cohorts share one reference, so the harmonization is the identity) and the pipeline's meta rows
      (ivw_meta.tsv), for the independent IVW recalculation in check_central_ivw.R
  check_central.py check CENTRAL_DIR TRUTH.tsv
      planted signals (meta, features, fine-mapping, colocalization, G x C, GxE tiers 1-2), tier-1 scope
      labels, ingest report, TWAS input files and the summary. Exit 1 on any failure.
"""
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd

ARM = {"expression": "expression", "splicing": "splicing_ungrouped"}
COLOC_LEAD = "chr22:19730097:A:G"           # tests/fixtures/central/make_gwas.py
COLOC_GENE, H3_GENE = "ENSGFIX0000011", "ENSGFIX0000006"
SCOPES = {"single_study", "shared_concordant_lead", "shared_discordant_lead"}


def read_any(stem):
    for ext in (".parquet", ".tsv.gz"):
        if os.path.exists(stem + ext):
            return pd.read_parquet(stem + ext) if ext == ".parquet" else pd.read_csv(stem + ext, sep="\t")
    return None


def tsv(path):
    return pd.read_csv(path, sep="\t")


def extract(a):
    cen = a.central_dir
    con = tsv(os.path.join(cen, "ingest", "contributions.tsv"))
    metas = {os.path.basename(p)[:-8]: pd.read_parquet(p) for p in glob.glob(os.path.join(cen, "meta", "*.parquet"))}
    truth = tsv(a.truth) if a.truth else pd.DataFrame(columns=["phenotype_id", "variant_id", "modality"])
    rng = np.random.default_rng(a.seed)
    picks = []
    for unit, m in sorted(metas.items()):
        m = m.copy()
        m["unit"] = unit
        planted = m[m.set_index(["phenotype_id", "variant_id"]).index.isin(
            list(zip(truth["phenotype_id"], truth["variant_id"])))]
        # random rows with k >= 2 in every scope, including pairs missing from some contributions
        rand = []
        for _, s in m[m["k"] >= 2].groupby("scope"):
            rand.append(s.iloc[rng.choice(len(s), size=min(a.pairs, len(s)), replace=False)])
        picks.append(pd.concat([planted] + rand))
    meta = pd.concat(picks).drop_duplicates(["unit", "scope", "phenotype_id", "variant_id"])
    rows = []
    for r in con.itertuples():
        for unit in sorted(set(meta["unit"])):
            # export_dir is the staged copy inside the Nextflow work dir; source_dir is the study's export
            stem = os.path.join(r.source_dir, r.ancestry, "A_nominal", f"{r.cohort}_{r.ancestry}_A_nominal_{unit}")
            d = read_any(stem)
            if d is None:
                continue
            want = meta[meta["unit"] == unit][["phenotype_id", "variant_id"]].drop_duplicates()
            d = d.merge(want, on=["phenotype_id", "variant_id"])
            d = d[["phenotype_id", "variant_id", "beta", "se", "n", "eaf"]].copy()
            d.insert(0, "unit", unit)
            d["contrib"], d["ancestry"] = r.contrib, r.ancestry
            rows.append(d)
    os.makedirs(a.out_dir, exist_ok=True)
    pd.concat(rows).to_csv(os.path.join(a.out_dir, "ivw_inputs.tsv"), sep="\t", index=False)
    meta[["unit", "scope", "phenotype_id", "variant_id", "k", "beta_fe", "se_fe", "p_fe", "beta_re", "se_re", "p_re",
          "tau2", "q_het", "i2", "contrib_mask"]].to_csv(os.path.join(a.out_dir, "ivw_meta.tsv"), sep="\t", index=False)
    print(f"extracted {len(meta)} meta rows ({meta['scope'].value_counts().to_dict()})")


def check(a):
    cen, truth = a.central_dir, tsv(a.truth)
    fails = 0

    def ok(cond, msg):
        nonlocal fails
        print(("PASS  " if cond else "FAIL  ") + msg)
        fails += 0 if cond else 1

    # ---- ingest ------------------------------------------------------------------------------------
    rep = tsv(os.path.join(cen, "ingest", "ingest_report.tsv"))
    ok(not (rep["status"] == "FAIL").any(), f"ingest report has no FAIL ({rep['status'].value_counts().to_dict()})")
    con = tsv(os.path.join(cen, "ingest", "contributions.tsv"))
    got = sorted(zip(con["study"], con["ancestry"]))
    ok(len(con) == 4 and con["study"].nunique() == 2 and set(con["ancestry"]) == {"EUR", "EAS"},
       f"4 contributions, 2 studies x EUR/EAS: {got}")
    hs = tsv(os.path.join(cen, "ingest", "harmonization_summary.tsv"))
    bad = hs[hs["action"].isin(["mismatch", "drop_palindromic"])]["n"].sum()
    ok(bad == 0, f"harmonization: no mismatched or dropped variants in reference-aligned exports ({bad})")

    # ---- planted cis effects: meta and feature level --------------------------------------------
    meta = pd.concat([pd.read_parquet(p) for p in glob.glob(os.path.join(cen, "meta", "*.parquet"))])
    feats = pd.read_csv(os.path.join(cen, "features.tsv.gz"), sep="\t")
    for t in truth.itertuples():
        m = meta[(meta["phenotype_id"] == t.phenotype_id) & (meta["variant_id"] == t.variant_id)].set_index("scope")
        if "ALL" not in m.index:
            ok(False, f"{t.phenotype_id} @ {t.variant_id}: no ALL meta row")
            continue
        al = m.loc["ALL"]
        sign_ok = np.sign(al["beta_fe"]) == np.sign(t.effect_per_alt)
        ok(al["k"] == 4 and al["p_fe"] <= 1e-6 and sign_ok,
           f"meta ALL {t.phenotype_id} @ {t.variant_id}: k={al['k']} beta={al['beta_fe']:.3f} p={al['p_fe']:.2e}")
        within = [m.loc[s, "p_fe"] for s in ("EUR", "EAS") if s in m.index]
        ok(len(within) == 2 and al["p_fe"] < min(within),
           f"meta ALL p < within-ancestry p for {t.phenotype_id} ({al['p_fe']:.1e} < {min(within):.1e})")
        f = feats[(feats["scope"] == "ALL") & (feats["arm"] == ARM[t.modality]) & (feats["phenotype_id"] == t.phenotype_id)]
        ok(len(f) == 1 and bool(f["significant"].iloc[0]) and int(f["k"].iloc[0]) == 4,
           f"feature {ARM[t.modality]} {t.phenotype_id}: significant across 4 contributions "
           f"(q={f['q_acat'].iloc[0] if len(f) else float('nan'):.1e})")

    # ---- fine-mapping ------------------------------------------------------------------------------
    fm = os.path.join(cen, "finemap")
    status, cs, loci = tsv(os.path.join(fm, "status.tsv")), tsv(os.path.join(fm, "cs.tsv")), tsv(os.path.join(fm, "loci.tsv"))
    for t in truth.itertuples():
        loc = loci[(loci["arm"] == ARM[t.modality]) & (loci["phenotype_id"] == t.phenotype_id)]
        if loc.empty:
            ok(False, f"fine-mapping locus for {t.phenotype_id}")
            continue
        lid = loc["locus_id"].iloc[0]
        st = status[status["locus_id"] == lid]
        in_cs = t.variant_id in set(cs.loc[cs["locus_id"] == lid, "variant_id"])
        ok(len(st) == 1 and st["status"].iloc[0] == "ok" and loc["ancestries"].iloc[0] == "EAS;EUR" and in_cs,
           f"fine-mapping {lid}: both ancestries, planted variant in a credible set")

    # ---- colocalization ----------------------------------------------------------------------------
    co = tsv(os.path.join(cen, "coloc", "coloc.tsv"))
    pos = co[co["phenotype_id"] == COLOC_GENE]
    ok(len(pos) == 1 and bool(pos["coloc"].iloc[0]) and pos["hit_xqtl"].iloc[0] == COLOC_LEAD,
       f"SuSiE-coloc {COLOC_GENE}: PP.H4={pos['PP.H4'].iloc[0] if len(pos) else float('nan'):.3f} at {COLOC_LEAD}")
    # coloc.tsv keeps the credible-set pair with the highest PP.H4 per locus; the H3 evidence is in coloc_pairs
    neg = co[co["phenotype_id"] == H3_GENE]
    pairs = tsv(os.path.join(cen, "coloc", "coloc_pairs.tsv"))
    npairs = pairs[pairs["phenotype_id"] == H3_GENE]
    h3 = npairs.loc[npairs["PP.H3"].idxmax()] if len(npairs) else None
    ok(len(neg) == 1 and not bool(neg["coloc"].iloc[0]) and h3 is not None and float(h3["PP.H3"]) > 0.5
       and h3["hit_xqtl"] != h3["hit_gwas"] and float(npairs["PP.H4"].max()) < 0.5,
       f"SuSiE-coloc {H3_GENE}: not colocalized (max PP.H4={npairs['PP.H4'].max() if len(npairs) else float('nan'):.3f}); "
       f"distinct signals {h3['hit_xqtl'] if h3 is not None else 'NA'} vs {h3['hit_gwas'] if h3 is not None else 'NA'} "
       f"(PP.H3={h3['PP.H3'] if h3 is not None else float('nan'):.3f})")
    pp = ["PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4"]
    okrows = co[co["status"] == "ok"]
    ok(len(okrows) == len(co) and okrows[pp].notna().all().all()
       and np.allclose(okrows[pp].sum(axis=1), 1, atol=1e-6),
       f"coloc: {len(co)} tasks ok, posteriors sum to 1 ({okrows['method'].value_counts().to_dict()})")
    cb_path = os.path.join(cen, "coloc", "colocboost.tsv")
    cb = tsv(cb_path) if os.path.exists(cb_path) else pd.DataFrame(columns=["phenotype_id", "status"])
    cbp = cb[(cb["phenotype_id"] == COLOC_GENE) & (cb["status"] == "ok")]
    ok(len(cbp) >= 1 and cbp["colocalized_variables"].astype(str).str.contains(COLOC_LEAD, regex=False).any(),
       f"ColocBoost {COLOC_GENE}: colocalization set contains {COLOC_LEAD}")

    # ---- G x C -------------------------------------------------------------------------------------
    it = os.path.join(cen, "interactions")
    gx = pd.read_csv(os.path.join(it, "gxc_annotation.tsv.gz"), sep="\t")
    for t in truth[truth["interaction"].str.startswith("gxc:")].itertuples():
        ct = t.interaction.split(":", 1)[1]
        r = gx[(gx["scope"] == "ALL") & (gx["phenotype_id"] == t.phenotype_id) & (gx["variant_id"] == t.variant_id)]
        ok(len(r) == 1 and r["primary_cell_type"].iloc[0] == ct and r["primary_method"].iloc[0] == "ivw"
           and int(r["primary_k"].iloc[0]) >= 2,
           f"G x C {t.phenotype_id}: central primary cell type {r['primary_cell_type'].iloc[0] if len(r) else 'NA'} "
           f"(planted {ct}), k={r['primary_k'].iloc[0] if len(r) else 0}")

    # ---- GxE ---------------------------------------------------------------------------------------
    t1 = pd.read_csv(os.path.join(it, "gxe_tier1_meta.tsv.gz"), sep="\t")
    ok(set(t1["scope"]) <= SCOPES, f"tier-1 scope labels valid: {sorted(set(t1['scope']))}")
    single_ok = ((t1["n_contributions"] == 1) == (t1["scope"] == "single_study")).all()
    conc = t1[t1["scope"] == "shared_concordant_lead"]
    ok(single_ok and conc["p_meta"].notna().all() and t1.loc[t1["scope"] == "shared_discordant_lead", "p_meta"].isna().all(),
       "tier-1 scopes consistent with n_contributions; IVW only for concordant leads")
    t2 = pd.read_csv(os.path.join(it, "gxe_tier2_top.tsv.gz"), sep="\t")
    for t in truth[(truth["interaction"] == "gxe:GA") & (truth["modality"] == "expression")].itertuples():
        r = t1[(t1["exposure"] == "GA") & (t1["phenotype_id"] == t.phenotype_id)]
        ok(len(r) == 1 and r["scope"].iloc[0] != "single_study" and float(r["q_feature_acat"].iloc[0]) <= 0.05,
           f"GxE tier 1 {t.phenotype_id} x GA: scope {r['scope'].iloc[0] if len(r) else 'NA'}, "
           f"q={r['q_feature_acat'].iloc[0] if len(r) else float('nan'):.1e}")
        if len(r) and r["scope"].iloc[0] == "shared_concordant_lead":
            ok(r["meta_variant_id"].iloc[0] == t.variant_id, f"tier-1 concordant lead is the planted {t.variant_id}")
        s = t2[(t2["exposure"] == "GA") & (t2["phenotype_id"] == t.phenotype_id)]
        ok(len(s) == 1 and s["variant_id"].iloc[0] == t.variant_id and float(s["p_fe"].iloc[0]) <= 1e-4
           and int(s["k"].iloc[0]) == 4 and np.sign(s["beta_fe"].iloc[0]) == np.sign(t.interaction_effect),
           f"GxE tier 2 {t.phenotype_id} x GA: top variant {s['variant_id'].iloc[0] if len(s) else 'NA'}, "
           f"p={s['p_fe'].iloc[0] if len(s) else float('nan'):.1e}, k={s['k'].iloc[0] if len(s) else 0}")
    sc = tsv(os.path.join(it, "gxe_exposure_scales.tsv"))
    ok(bool(sc.loc[sc["exposure"] == "GA", "rescale"].iloc[0]) and float(sc.loc[sc["exposure"] == "GA", "pooled_sd"].iloc[0]) > 0,
       "GA rescaled to the pooled SD")

    # ---- TWAS inputs -------------------------------------------------------------------------------
    tw = os.path.join(cen, "twas")
    for anc in ("EUR", "EAS"):
        ss = tsv(os.path.join(tw, "otters", anc, "eQTLSumStats.txt"))
        an = tsv(os.path.join(tw, "otters", anc, "exp_anno.txt"))
        ok(list(ss.columns) == ["CHROM", "POS", "A1", "A2", "Zscore", "TargetID", "N"]
           and list(an.columns) == ["CHROM", "GeneStart", "GeneEnd", "TargetID"]
           and set(ss["TargetID"]) <= set(an["TargetID"]) and "expression:" + H3_GENE in set(an["TargetID"]),
           f"OTTERS {anc}: formats, every TargetID annotated")
        # weights are written for every locus SuShiE fitted (status ok or no_cs)
        n_fit = status[status["status"].isin(["ok", "no_cs"])]["locus_id"].nunique()
        w = glob.glob(os.path.join(tw, "fusion", anc, "WEIGHTS", "*.wgt.RDat"))
        posf = tsv(os.path.join(tw, "fusion", anc, f"{anc}.pos"))
        ok(len(w) == n_fit and len(posf) == len(w), f"FUSION {anc}: {len(w)} weight files for {n_fit} fitted loci")

    # ---- summary -----------------------------------------------------------------------------------
    sm = tsv(os.path.join(cen, "summary", "central_summary.tsv"))
    v = sm[(sm["section"] == "ingest") & (sm["item"] == "contributions")]["value"]
    ok(len(v) == 1 and int(v.iloc[0]) == 4, "summary: 4 contributions")
    en = tsv(os.path.join(cen, "summary", "pip_enrichment.tsv"))
    e = en[en["annotation"] == "planted_regions"]
    ok(len(e) == 1 and int(e["n_high_pip_in_class"].iloc[0]) >= 1 and float(e["odds_ratio"].iloc[0]) > 1,
       f"PIP enrichment in planted regions: OR={e['odds_ratio'].iloc[0] if len(e) else float('nan')}")

    print(f"check_central.py: {fails} failure(s)")
    return 1 if fails else 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("extract")
    e.add_argument("central_dir")
    e.add_argument("out_dir")
    e.add_argument("--truth")
    e.add_argument("--pairs", type=int, default=30)
    e.add_argument("--seed", type=int, default=1)
    c = sub.add_parser("check")
    c.add_argument("central_dir")
    c.add_argument("truth")
    a = ap.parse_args(argv)
    if a.cmd == "extract":
        extract(a)
        return 0
    return check(a)


if __name__ == "__main__":
    sys.exit(main())
