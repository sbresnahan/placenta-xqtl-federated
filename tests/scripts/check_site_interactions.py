#!/usr/bin/env python3
"""check_site_interactions.py EXPORT_COHORT_DIR ANC...  -- interaction checks (G x C and G x E, files G/H).

Per stratum, against the plants in tests/fixtures/make_site_fixture.py and the restricted test signals:
  G  * ENSGFIX0000011 (planted G x Cytotrophoblast) has primary cell type Cytotrophoblast, q <= 0.05;
     * test-set sources and untested reasons as constructed: ENSG00000100000 phenotype_absent,
       ENSGFIX0000003 chr22:20000001:A:C variant_absent, ENSGFIX0000006 one site and one mdacc pair,
       ENSGFIX0000011 both, clu_31 intron 1 both, ENSGFIX0000035 mdacc only.
  H  * tier 1: ENSGFIX0000021 x GA (planted, expression) q <= 0.05 with b_gi > 0;
     * tier 2: the planted variant chr22:20290948:A:G is the top variant of ENSGFIX0000021;
     * the splicing plant (clu_31 intron 1) is tested in tier 1 (main-effect / flag target only:
       splicing normalization attenuates a splicing interaction below detection at n = 60);
     * H_status: tier1, tier2 and sensitivity DONE for expression x GA;
     * H_sensitivity: maternal_age refit of ENSGFIX0000021 present and stable.
Only these pairs are asserted: other assignments at q <= 0.05 can occur by chance in 60 samples.
Exit 1 on any failure.
"""
import os
import sys

import pandas as pd

GXE_GENE, GXE_VAR = "ENSGFIX0000021", "chr22:20290948:A:G"
CLU31 = "ENSGFIX0000031__chr22:20762038:20762738:clu_31_+"


def read(path_noext):
    for ext, f in ((".parquet", pd.read_parquet), (".tsv.gz", lambda p: pd.read_csv(p, sep="\t"))):
        if os.path.exists(path_noext + ext):
            return f(path_noext + ext)
    raise FileNotFoundError(path_noext)


def main():
    root, ancs = sys.argv[1], sys.argv[2:]
    cohort = os.path.basename(root.rstrip("/"))
    fails = 0

    def check(ok, msg):
        nonlocal fails
        print(("PASS  " if ok else "FAIL  ") + msg)
        fails += 0 if ok else 1

    for anc in ancs:
        pre = os.path.join(root, anc, f"{cohort}_{anc}")
        # ---------------- G ----------------
        ann = pd.read_csv(pre + "_G_gxc_annotation.tsv.gz", sep="\t")
        r = ann[ann.phenotype_id == "ENSGFIX0000011"]
        ok = len(r) == 1 and r.primary_cell_type.iloc[0] == "Cytotrophoblast" and r.primary_q_GxC.iloc[0] <= 0.05
        q = r.primary_q_GxC.iloc[0] if len(r) else float("nan")
        check(ok, f"{anc} G x C ENSGFIX0000011 primary={r.primary_cell_type.tolist()} q={q:.1e}")
        g = read(pre + "_G_gxc")
        check(len(g) > 0 and {"cell_type", "beta_GxC", "q_GxC"} <= set(g.columns), f"{anc} G table: {len(g)} rows")

        def src(ph, var=None):
            x = ann[(ann.phenotype_id == ph) & ((ann.variant_id == var) if var else True)]
            return sorted(x.source), sorted(x.untested_reason.dropna())
        exp = {("ENSG00000100000", None): (["mdacc"], ["phenotype_absent"]),
               ("ENSGFIX0000003", "chr22:20000001:A:C"): (["mdacc"], ["variant_absent"]),
               ("ENSGFIX0000006", None): (["mdacc", "site"], []),
               ("ENSGFIX0000011", None): (["both"], []),
               (CLU31, None): (["both"], []),
               ("ENSGFIX0000035", None): (["mdacc"], [])}
        for (ph, var), want in exp.items():
            got = src(ph, var)
            check(got == want, f"{anc} test set {ph.split('__')[0]}{' ' + var if var else ''}: {got}")

        # ---------------- H ----------------
        t1 = read(pre + "_H_gxe_tier1")
        r = t1[(t1.modality == "expression") & (t1.exposure == "GA") & (t1.phenotype_id == GXE_GENE)]
        ok = len(r) == 1 and r.qval.iloc[0] <= 0.05 and r.b_gi.iloc[0] > 0
        desc = f"q={r.qval.iloc[0]:.1e} b_gi={r.b_gi.iloc[0]:+.2f} lead={r.variant_id.iloc[0]}" if len(r) else "missing"
        check(ok, f"{anc} G x E tier 1 {GXE_GENE} x GA: {desc}")
        check(((t1.modality == "splicing") & (t1.phenotype_id == CLU31)).any(), f"{anc} tier 1 tests clu_31 intron 1")
        t2 = read(os.path.join(root, anc, "H_gxe_tier2", f"{cohort}_{anc}_H_gxe_tier2_expression_GA"))
        top = t2[(t2.phenotype_id == GXE_GENE) & t2.top_variant]
        check(len(top) == 1 and top.variant_id.iloc[0] == GXE_VAR,
              f"{anc} tier 2 top variant of {GXE_GENE}: {top.variant_id.tolist()} (planted {GXE_VAR})")
        st = pd.read_csv(pre + "_H_status.tsv", sep="\t")
        s = st[(st.modality == "expression") & (st.exposure == "GA")].set_index("component").status.to_dict()
        check(all(s.get(c) == "DONE" for c in ("tier1", "tier2", "sensitivity")), f"{anc} H_status expression x GA: {s}")
        se = pd.read_csv(pre + "_H_sensitivity.tsv", sep="\t")
        r = se[(se.phenotype_id == GXE_GENE) & (se.covariates == "maternal_age")]
        check(len(r) == 1 and bool(r.stable.iloc[0]), f"{anc} sensitivity maternal_age {GXE_GENE}: stable={r.stable.tolist()}")
        check(os.path.exists(pre + "_H_exposures.tsv"), f"{anc} H_exposures present")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
