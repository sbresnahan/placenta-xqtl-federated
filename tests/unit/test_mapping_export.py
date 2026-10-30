"""Unit tests for the M3 helpers: leak scanner, export finalize, covariate pruning, k selection,
modality combination and parquet merging (no tensorQTL / R needed)."""
import gzip
import json
import os

import numpy as np
import pandas as pd
import pytest

from conftest import load

export = load("xqtl_export")
covariates = load("xqtl_covariates")
select_k = load("xqtl_select_k")
combine = load("xqtl_combine_modalities")
merge = load("xqtl_merge_parquet")


# ---------------- leak scanner ----------------
def _ids(tmp_path, ids):
    p = tmp_path / "metadata.tsv"
    pd.DataFrame({"rnaseq_id": ids, "genotype_id": [f"G{i}" for i in range(len(ids))]}).to_csv(p, sep="\t", index=False)
    return export.load_ids([str(p)])


def test_scanner_token_boundaries(tmp_path):
    sc = export.Scanner(_ids(tmp_path, ["S001", "HG00351", "12"]))
    assert sc.hit("S001") == "S001"
    assert sc.hit("total_reads_S001") == "S001"          # '_' is a boundary
    assert sc.hit("ENSGFIX0000006") is None              # embedded digits are not an ID
    assert sc.hit("S0010") is None
    assert sc.hit("12") == "12"                          # short/numeric IDs: exact match only
    assert sc.hit("chr12") is None


def test_scanner_skips_variant_ids(tmp_path):
    sc = export.Scanner(_ids(tmp_path, ["22"]))
    sc.strings(["chr22:100:A:G"], "x")
    assert sc.findings == []


def _clean_export(d):
    os.makedirs(d / "B_mapping")
    pd.DataFrame({"phenotype_id": ["ENSG1"], "variant_id": ["chr1:5:A:G"], "pval": [0.1]}).to_parquet(
        d / "B_mapping" / "C_EUR_B_permutation_expression.parquet", index=False)
    np.savez(d / "C_EUR_C_ld_chr1_0_3000000.npz", R=np.eye(2, dtype=np.float16))
    with open(d / "C_EUR_F_model_settings.json", "w") as fh:
        json.dump({"covariates": ["sex", "PC1", "ct_Hofbauer"]}, fh)


@pytest.mark.parametrize("mutation, expect", [
    (None, None),
    ("column", "S001"),
    ("value", "G1"),
    ("filename", "S002"),
    ("json", "S003"),
    ("npz", "npz must hold exactly one square matrix"),
    ("ext", "unexpected file type"),
])
def test_scan_dir_findings(tmp_path, mutation, expect):
    ids = _ids(tmp_path, ["S001", "S002", "S003"])
    d = tmp_path / "EUR"
    _clean_export(d)
    if mutation == "column":
        pd.DataFrame({"S001": [1.0]}).to_parquet(d / "C_EUR_I_summary.parquet", index=False)
    elif mutation == "value":
        pd.DataFrame({"item": ["G1"]}).to_csv(d / "C_EUR_I_summary.tsv.gz", sep="\t", index=False, compression="gzip")
    elif mutation == "filename":
        (d / "S002_notes.tsv").write_text("a\n1\n")
    elif mutation == "json":
        (d / "C_EUR_F_model_settings.json").write_text(json.dumps({"note": "dropped S003"}))
    elif mutation == "npz":
        np.savez(d / "C_EUR_C_ld_x.npz", R=np.eye(2), G=np.ones((10, 2)))
    elif mutation == "ext":
        (d / "dosages.csv").write_text("x\n")
    found = export.scan_dir(str(d), ids)
    if expect is None:
        assert found == []
    else:
        assert any(expect in f[1] for f in found), found


def test_finalize_deletes_export_on_leak(tmp_path):
    ids_file = tmp_path / "metadata.tsv"
    pd.DataFrame({"rnaseq_id": ["S001"], "genotype_id": ["HG1"]}).to_csv(ids_file, sep="\t", index=False)
    strata = tmp_path / "strata.tsv"
    pd.DataFrame({"ancestry": ["EUR"], "n_individuals": [60], "n_rnaseq": [60], "status": ["mapped"]}).to_csv(
        strata, sep="\t", index=False)
    d = tmp_path / "EUR"
    _clean_export(d)
    out = tmp_path / "export"
    export.main(["finalize", "--cohort", "C", "--strata", str(strata), "--ids", str(ids_file),
                 "--leak-report", str(tmp_path / "leak.tsv"), "--out", str(out), str(d)])
    man = pd.read_csv(out / "C" / "MANIFEST.tsv", sep="\t")
    assert {"C/EUR/C_EUR_C_ld_chr1_0_3000000.npz", "C/C_cohort_I_strata.tsv", "C/DATA_DICTIONARY.md"} <= set(man.path)
    # plant a sample ID: the export must fail and be deleted
    pd.DataFrame({"phenotype_id": ["HG1"]}).to_parquet(d / "C_EUR_E_features.parquet", index=False)
    with pytest.raises(SystemExit) as e:
        export.main(["finalize", "--cohort", "C", "--strata", str(strata), "--ids", str(ids_file),
                     "--leak-report", str(tmp_path / "leak.tsv"), "--out", str(out), str(d)])
    assert "LEAK SCAN FAILED" in str(e.value)
    assert not out.exists()
    assert "HG1" in (tmp_path / "leak.tsv").read_text()


# ---------------- covariates ----------------
def test_prune_correlated_keeps_higher_priority():
    rng = np.random.default_rng(1)
    x = rng.normal(size=50)
    cov = pd.DataFrame([x, x + rng.normal(scale=0.01, size=50), rng.normal(size=50)],
                       index=["PC1", "HCP1", "sex"])
    pri = {"sex": 0, "PC1": 1, "HCP1": 3}
    kept, dropped = covariates.prune_correlated(cov, pri, 0.9)
    assert list(kept.index) == ["PC1", "sex"]
    assert dropped[0]["covariate"] == "HCP1" and dropped[0]["correlated_with"] == "PC1"


# ---------------- k selection ----------------
def test_select_k_argmax_ties_to_smaller(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    evs = []
    for k, n in ((0, 1), (10, 5), (20, 5), (30, 4)):
        p = tmp_path / f"e{k}.tsv"
        pd.DataFrame([{"k": k, "scope": "chr1", "n_tested": 400, "n_egenes": n, "n_covariates": 8 + k}]).to_csv(
            p, sep="\t", index=False)
        evs.append(str(p))
    select_k.main(["--anc", "EUR", "--modality", "splicing", "--arm", "splicing", "--arm", "splicing_ungrouped",
                   "--grouped", "--n-samples", "80", "--nperm", "10000", "--window", "1000000", "--maf", "0.01",
                   "--fdr", "0.05", "--seed", "1", "--out-k", "k.txt"] + evs)
    assert (tmp_path / "k.txt").read_text().strip() == "10"
    g = json.load(open(tmp_path / "EUR_splicing.settings.json"))
    u = json.load(open(tmp_path / "EUR_splicing_ungrouped.settings.json"))
    assert g["grouped"] and not u["grouped"] and g["k_selected"] == u["k_selected"] == 10
    assert [r["k"] for r in g["hcp_grid"]] == [0, 10, 20, 30]


# ---------------- combined modalities ----------------
def _bed(path, ids, samples, start):
    df = pd.DataFrame({"#chr": "chr1", "start": start, "end": [s + 1 for s in start], "phenotype_id": ids})
    for i, s in enumerate(samples):
        df[s] = np.arange(len(ids)) + i
    df.to_csv(path, sep="\t", index=False, compression="gzip")


def test_combine_modalities(tmp_path):
    shared = [f"X{i}" for i in range(10)]
    _bed(tmp_path / "e.bed.gz", ["G1", "G2"], ["A"] + shared, [100, 200])
    _bed(tmp_path / "s.bed.gz", ["G1__intron1", "G1__intron2"], shared + ["D"], [100, 100])
    combine.main(["--out", str(tmp_path / "X"), f"expression={tmp_path / 'e.bed.gz'}",
                  f"splicing={tmp_path / 's.bed.gz'}"])
    with gzip.open(tmp_path / "X.bed.gz", "rt") as fh:
        bed = pd.read_csv(fh, sep="\t")
    assert list(bed.columns[4:]) == shared                          # sample intersection
    assert set(bed.phenotype_id) == {"expression__G1", "expression__G2", "splicing__G1__intron1", "splicing__G1__intron2"}
    grp = pd.read_csv(tmp_path / "X.phenotype_groups.txt", sep="\t", header=None, names=["p", "g"])
    assert grp.groupby("g").size().to_dict() == {"G1": 3, "G2": 1}


def test_merge_parquet(tmp_path):
    pd.DataFrame({"phenotype_id": ["a"], "variant_id": ["v"]}).to_parquet(tmp_path / "c1.parquet", index=False)
    pd.DataFrame(columns=["phenotype_id", "variant_id"]).to_parquet(tmp_path / "c2.parquet", index=False)
    merge.main([str(tmp_path / "o.parquet"), str(tmp_path / "c1.parquet"), str(tmp_path / "c2.parquet")])
    assert len(pd.read_parquet(tmp_path / "o.parquet")) == 1
    merge.main([str(tmp_path / "o2.parquet"), str(tmp_path / "c2.parquet")])
    assert len(pd.read_parquet(tmp_path / "o2.parquet")) == 0
