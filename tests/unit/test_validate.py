"""Unit tests for bin/xqtl_validate.py (cohort.yml + samplesheet + metadata validation)."""
import csv
import json
import os
import shutil

import pytest

from conftest import FIX, SCHEMA, load

V = load("xqtl_validate")


def base_params(tmp_path, **over):
    p = {
        "cohort_id": "TESTCOHORT",
        "mode": "site",
        "samplesheet": os.path.join(FIX, "samplesheet.csv"),
        "metadata": os.path.join(FIX, "metadata.tsv"),
        "restricted": os.path.join(FIX, "xqtl-restricted-test.tar.gz"),
        "site_refs": os.path.join(FIX, "site_refs"),
        "genotypes": {"imputed_dir": os.path.join(FIX, "genotypes")},
        "expected_restricted_version": "1",
        "max_cpus": 16,  # Nextflow params outside the schema are allowed
    }
    p.update(over)
    return p


def run(tmp_path, params, mode="site"):
    pj = tmp_path / "p.json"
    pj.write_text(json.dumps(params))
    oj, orep = tmp_path / "r.json", tmp_path / "r.tsv"
    rc = V.main(["--params", str(pj), "--schema", SCHEMA, "--mode", mode,
                 "--out-json", str(oj), "--out-report", str(orep), "--no-file-check"])
    rows = list(csv.DictReader(open(orep), delimiter="\t"))
    return rc, json.loads(oj.read_text()), rows


def test_fixture_is_valid_and_defaults_filled(tmp_path):
    rc, resolved, rows = run(tmp_path, base_params(tmp_path))
    assert rc == 0, rows
    assert resolved["qc"]["rsq_min"] == 0.4
    assert resolved["qc"]["maf_min"] == 0.01
    assert resolved["qc"]["mac_min"] == 5
    assert resolved["ancestry"]["min_n"] == 50
    assert resolved["mapping"]["fdr"] == 0.05
    assert resolved["gxe"]["maf_interaction"] == 0.05
    assert "Maternal" in resolved["gxc"]["exclude_cell_types"]


def test_missing_restricted_fails_with_email_hint(tmp_path, capsys):
    rc, _, rows = run(tmp_path, base_params(tmp_path, restricted=None))
    assert rc == 1
    msg = [r["detail"] for r in rows if r["check"] == "required.restricted"][0]
    assert "stbresnahan@mdanderson.org" in msg
    assert "stbresnahan@mdanderson.org" in capsys.readouterr().err


def test_wrong_restricted_version_fails(tmp_path):
    rc, _, rows = run(tmp_path, base_params(tmp_path, expected_restricted_version="2"))
    assert rc == 1
    assert any("does not match" in r["detail"] for r in rows)


def test_corrupt_restricted_fails(tmp_path):
    bad = tmp_path / "bad.tar.gz"
    bad.write_bytes(b"not an archive")
    rc, _, rows = run(tmp_path, base_params(tmp_path, restricted=str(bad)))
    assert rc == 1


def test_bad_cohort_id_and_threshold(tmp_path):
    rc, _, rows = run(tmp_path, base_params(tmp_path, cohort_id="1 bad id", qc={"rsq_min": 2}))
    assert rc == 1
    checks = {r["check"] for r in rows if r["status"] == "FAIL"}
    assert any(c.startswith("schema.cohort_id") for c in checks)
    assert any(c.startswith("schema.qc") for c in checks)


def test_ancestry_min_n_floor(tmp_path):
    rc, _, _ = run(tmp_path, base_params(tmp_path, ancestry={"min_n": 20}))
    assert rc == 1


def write_tsv(path, header, rows, sep="\t"):
    with open(path, "w") as fh:
        fh.write(sep.join(header) + "\n")
        for r in rows:
            fh.write(sep.join(r) + "\n")


def test_metadata_errors_are_all_reported(tmp_path):
    meta = tmp_path / "meta.tsv"
    write_tsv(meta, ["rnaseq_id", "genotype_id", "sex", "GA"],
              [["S01", "G01", "X", "39"], ["S02", "G02", "M", "12"], ["S03", "G03", "F", "abc"],
               ["S04", "G04", "F", "39"], ["S05", "G05", "F", "39"], ["S06", "G06", "F", "39"]])
    rc, _, rows = run(tmp_path, base_params(tmp_path, metadata=str(meta)))
    assert rc == 1
    fails = " ".join(r["detail"] for r in rows if r["status"] == "FAIL")
    assert "S01" in fails and "S02" in fails and "S03" in fails


def test_missing_required_metadata_column(tmp_path):
    meta = tmp_path / "meta.tsv"
    write_tsv(meta, ["rnaseq_id", "genotype_id", "GA"], [["S01", "G01", "39"]])
    rc, _, rows = run(tmp_path, base_params(tmp_path, metadata=str(meta)))
    assert rc == 1
    assert any("sex" in r["detail"] for r in rows if r["status"] == "FAIL")


def test_duplicate_sample_ids(tmp_path):
    ss = tmp_path / "ss.csv"
    shutil.copy(os.path.join(FIX, "samplesheet.csv"), ss)
    with open(ss, "a") as fh:
        fh.write("S01,a.fq.gz,b.fq.gz\n")
    rc, _, rows = run(tmp_path, base_params(tmp_path, samplesheet=str(ss)))
    assert rc == 1


def test_sample_without_metadata(tmp_path):
    ss = tmp_path / "ss.csv"
    shutil.copy(os.path.join(FIX, "samplesheet.csv"), ss)
    with open(ss, "a") as fh:
        fh.write("S99,a.fq.gz,b.fq.gz\n")
    rc, _, rows = run(tmp_path, base_params(tmp_path, samplesheet=str(ss)))
    assert rc == 1
    assert any("S99" in r["detail"] for r in rows if r["status"] == "FAIL")


def test_phenotyping_only_does_not_need_genotypes(tmp_path):
    p = base_params(tmp_path, mode="phenotyping_only")
    p.pop("genotypes")
    rc, _, rows = run(tmp_path, p, mode="phenotyping_only")
    assert rc == 0, rows


def test_imputation_prep_needs_raw_genotypes(tmp_path):
    rc, _, rows = run(tmp_path, {"cohort_id": "TC", "mode": "imputation_prep"}, mode="imputation_prep")
    assert rc == 1
    assert any(r["check"] == "required.genotypes.raw" for r in rows)


@pytest.mark.parametrize("mode", ["site", "phenotyping_only"])
def test_enabled_exposure_without_column_fails(tmp_path, mode):
    p = base_params(tmp_path, mode=mode,
                    gxe={"enabled": True, "exposures": {"ogtt": {"enabled": True, "type": "continuous"}}})
    rc, _, rows = run(tmp_path, p, mode=mode)
    if mode == "site":
        assert rc == 1
        assert any("ogtt" in r["detail"] for r in rows if r["status"] == "FAIL")


def test_user_exposure_extends_default_map(tmp_path):
    p = base_params(tmp_path, gxe={"exposures": {"ppBMI": {"enabled": True}}})
    rc, resolved, rows = run(tmp_path, p)
    ex = resolved["gxe"]["exposures"]
    assert ex["GA"]["enabled"] is True          # default kept
    assert ex["ppBMI"]["enabled"] is True       # user override
    assert ex["ppBMI"]["min"] == 14             # default fields of the same exposure kept
    assert rc == 0, rows


def test_sensitivity_default_off_and_missing_column_fails(tmp_path):
    rc, resolved, rows = run(tmp_path, base_params(tmp_path))
    sens = resolved["gxe"]["sensitivity"]
    assert not any(v["enabled"] for v in sens.values())
    assert rc == 0, rows
    p = base_params(tmp_path, gxe={"sensitivity": {"smoking": {"enabled": True}}})
    rc, resolved, rows = run(tmp_path, p)
    assert resolved["gxe"]["sensitivity"]["maternal_age"]["enabled"] is False   # defaults kept
    assert rc == 1
    assert any(r["check"] == "gxe.sensitivity.smoking" and r["status"] == "FAIL" for r in rows)


def _site_refs(tmp_path, read_length="101"):
    d = tmp_path / "refs"
    for f in V.SITE_REFS_FILES:
        (d / f).parent.mkdir(parents=True, exist_ok=True)
        (d / f).write_text("x")
    (d / "MANIFEST.tsv").write_text(f"key\tvalue\nread_length\t{read_length}\n")
    return str(d)


def test_site_refs_complete_and_read_length_warning(tmp_path):
    rep = V.Report()
    V.check_site_refs({"site_refs": _site_refs(tmp_path, "76"), "rnaseq": {"read_length": 101},
                       "mapping": {"modalities": ["expression", "rna_editing"]}}, rep)
    st = {c: s for c, s, _ in rep.rows}
    assert st["site_refs"] == "PASS"
    assert st["site_refs.read_length"] == "WARN"
    assert st["site_refs.editing"] == "WARN"


def test_site_refs_missing_or_incomplete_fails(tmp_path):
    rep = V.Report()
    V.check_site_refs({"site_refs": str(tmp_path / "nope")}, rep)
    assert rep.failed
    d = _site_refs(tmp_path)
    os.remove(os.path.join(d, "star", "SA"))
    rep = V.Report()
    V.check_site_refs({"site_refs": d}, rep)
    assert rep.failed and "star/SA" in rep.rows[-1][2]


def _ip_params(tmp_path, raw, build="GRCh38", **over):
    fa = tmp_path / "g.fa"
    fa.write_text(">chr22\nACGT\n")
    p = {"cohort_id": "TC", "mode": "imputation_prep", "genome_fasta": str(fa),
         "genotypes": {"raw": raw, "raw_build": build}}
    p.update(over)
    return p


def test_imputation_prep_plink_prefix_valid(tmp_path):
    raw = os.path.join(FIX, "imputation_prep", "kgp38_chr22")
    rc, resolved, rows = run(tmp_path, _ip_params(tmp_path, raw), mode="imputation_prep")
    assert rc == 0, rows
    assert resolved["imputation_prep"]["ref_match_min"] == 0.95
    assert any(r["check"] == "imputation_prep.raw" and r["status"] == "PASS" for r in rows)


def test_imputation_prep_grch37_needs_chain(tmp_path):
    raw = os.path.join(FIX, "imputation_prep", "kgp37_chr22.vcf.gz")
    p = _ip_params(tmp_path, raw, build="GRCh37", refs_dir=str(tmp_path / "norefs"))
    rc, _, rows = run(tmp_path, p, mode="imputation_prep")
    assert rc == 1
    assert any(r["check"] == "imputation_prep.liftover_chain" and r["status"] == "FAIL" for r in rows)
    p["liftover_chain"] = os.path.join(FIX, "imputation_prep", "hg19ToHg38.chr22.over.chain.gz")
    rc, _, rows = run(tmp_path, p, mode="imputation_prep")
    assert rc == 0, rows


def test_imputation_prep_bad_raw_and_wes_warning(tmp_path):
    p = _ip_params(tmp_path, str(tmp_path / "nothing"))
    p["genotypes"]["source"] = "wes"
    rc, _, rows = run(tmp_path, p, mode="imputation_prep")
    assert rc == 1
    assert any(r["check"] == "imputation_prep.raw" and r["status"] == "FAIL" for r in rows)
    assert any(r["check"] == "imputation_prep.source" and r["status"] == "WARN" for r in rows)
