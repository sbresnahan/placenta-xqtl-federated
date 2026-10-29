"""Tests for the MDACC packaging scripts, on a small synthetic mock of the atlas outputs.

Covers scripts/lib/atlas_restricted_tables.py, scripts/make_restricted_bundle.sh,
scripts/make_public_reference_record.sh (needs plink2 + samtools) and
scripts/pin_reference_record.py. All inputs are synthetic.
"""
import gzip
import hashlib
import importlib.util
import os
import shutil
import subprocess

import pandas as pd
import pytest

from conftest import REPO, load

R = load("xqtl_restricted")
spec = importlib.util.spec_from_file_location("art", os.path.join(REPO, "scripts", "lib", "atlas_restricted_tables.py"))
ART = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ART)


def _counts(path, junctions):
    with gzip.open(path, "wt") as fh:
        fh.write("S1 S2\n")
        for j in junctions:
            fh.write(f"{j} 3/10 4/12\n")


def mock_atlas(root):
    """EAS and EUR pooled splicing clusters plus tensorQTL tables for two modalities."""
    sp = root / "output"
    for anc, js in {
        "EAS": ["chr1:100:200:clu_1_+", "chr1:100:300:clu_1_+", "chr1:500:600:clu_2_-"],
        # clu_7 shares chr1:100:300 with EAS clu_1 -> same frozen cluster
        "EUR": ["chr1:100:300:clu_7_+", "chr1:150:300:clu_7_+", "chr2:10:20:clu_8_+"],
    }.items():
        d = sp / anc / "splicing" / "intermediate"
        d.mkdir(parents=True)
        _counts(d / "leafcutter_harmonized_perind_numers.counts.gz", js)
    q = root / "qtl_results"
    q.mkdir()
    # grouped splicing scan: the phenotype IDs are in the index, as tensorQTL writes them
    pd.DataFrame({"variant_id": ["chr1:150:A:G"], "qval": [0.01]},
                 index=pd.Index(["ENSG1__chr1:100:200:clu_1_+"], name=None)).to_parquet(q / "EAS_splicing_cisqtl.parquet")
    pd.DataFrame({"phenotype_id": ["ENSG1__chr1:150:300:clu_7_+", "ENSG2__chr2:10:20:clu_8_+"],
                  "variant_id": ["chr1:160:C:T", "chr2:15:G:A"], "qval": [0.30, 0.02]}).to_parquet(q / "EUR_splicing_cisqtl.parquet")
    pd.DataFrame({"phenotype_id": ["ENSG1__chr1:150:300:clu_7_+"], "variant_id": ["chr1:170:A:C"],
                  "rank": [2]}).to_parquet(q / "EUR_splicing_cisqtl_independent.parquet")
    pd.DataFrame({"phenotype_id": ["ENSG1", "ENSG3"], "variant_id": ["chr1:150:A:G", "chr3:1:A:T"],
                  "qval": [0.001, 0.5]}).to_parquet(q / "EAS_expression_cisqtl.parquet")
    pd.DataFrame({"phenotype_id": ["ENSG1", "ENSG1"], "variant_id": ["chr1:150:A:G", "chr1:900:G:C"],
                  "rank": [1, 2]}).to_parquet(q / "EAS_expression_cisqtl_independent.parquet")
    pd.DataFrame({"phenotype_id": ["ENSG1__ENST1"], "variant_id": ["chr1:150:A:G"],
                  "qval": [0.04]}).to_parquet(q / "EUR_isoforms_cisqtl.parquet")
    return sp, q


def test_frozen_clusters_merge_across_ancestries(tmp_path):
    sp, _ = mock_atlas(tmp_path)
    fc = ART.frozen_clusters({a: str(sp / a / "splicing/intermediate/leafcutter_harmonized_perind_numers.counts.gz")
                              for a in ("EAS", "EUR")})
    got = {(r.chrom, r.start, r.end, r.strand): r.cluster_id for r in fc.itertuples()}
    assert got[("chr1", 100, 200, "+")] == got[("chr1", 150, 300, "+")] == got[("chr1", 100, 300, "+")] == "clu_1"
    assert got[("chr1", 500, 600, "-")] == "clu_2"
    assert got[("chr2", 10, 20, "+")] == "clu_3"
    assert fc.loc[(fc.start == 100) & (fc.end == 300), "atlas_clusters"].item() == "EAS:clu_1;EUR:clu_7"


def test_signal_table(tmp_path):
    sp, q = mock_atlas(tmp_path)
    oc, os_ = tmp_path / "c.tsv.gz", tmp_path / "s.tsv.gz"
    rc = ART.main(["--ancestries", "EAS EUR", "--qtl-results", str(q), "--out-clusters", str(oc), "--out-signals", str(os_),
                   "--splicing-counts", str(sp / "{ANC}/splicing/intermediate/leafcutter_harmonized_perind_numers.counts.gz")])
    assert rc == 0
    s = pd.read_csv(os_, sep="\t")
    assert list(s.columns[:5]) == R.REQUIRED["signals/mdacc_signals.tsv.gz"]
    spl = s[s.modality == "splicing"].set_index("atlas_phenotype_id")
    # canonical IDs use the frozen cluster; EUR clu_7 -> clu_1, EUR clu_8 -> clu_3
    assert spl.loc["ENSG1__chr1:100:200:clu_1_+", "phenotype_id"] == "ENSG1__chr1:100:200:clu_1_+"
    assert spl.loc["ENSG1__chr1:150:300:clu_7_+", "phenotype_id"] == "ENSG1__chr1:150:300:clu_1_+"
    assert spl.loc["ENSG2__chr2:10:20:clu_8_+", "phenotype_id"] == "ENSG2__chr2:10:20:clu_3_+"
    # EUR lead with q = 0.30 is excluded; its independent rank-2 signal is kept
    eur1 = s[(s.ancestry == "EUR") & (s.atlas_phenotype_id == "ENSG1__chr1:150:300:clu_7_+")]
    assert eur1.variant_id.tolist() == ["chr1:170:A:C"] and eur1.signal_rank.tolist() == [2]
    expr = s[s.modality == "expression"]
    assert sorted(zip(expr.variant_id, expr.signal_rank)) == [("chr1:150:A:G", 1), ("chr1:900:G:C", 2)]
    assert "ENSG3" not in set(s.phenotype_id)              # q = 0.5
    assert set(s.modality) == {"splicing", "expression", "isoform_usage"}   # atlas 'isoforms' renamed
    # byte-reproducible output
    b1 = os_.read_bytes()
    ART.main(["--ancestries", "EAS EUR", "--qtl-results", str(q), "--out-clusters", str(oc), "--out-signals", str(os_),
              "--splicing-counts", str(sp / "{ANC}/splicing/intermediate/leafcutter_harmonized_perind_numers.counts.gz")])
    assert os_.read_bytes() == b1


def _env(tmp_path, sp, q):
    e = tmp_path / "mock.env"
    e.write_text(f'RESTRICTED_VERSION=7\nANCESTRIES="EAS EUR"\nQVALUE_MAX=0.05\n'
                 f'ATLAS_SPLICING_COUNTS="{sp}/{{ANC}}/splicing/intermediate/leafcutter_harmonized_perind_numers.counts.gz"\n'
                 f'ATLAS_QTL_RESULTS={q}\n')
    return e


def test_make_restricted_bundle_dry_run_and_build(tmp_path):
    sp, q = mock_atlas(tmp_path)
    env = _env(tmp_path, sp, q)
    sh = os.path.join(REPO, "scripts", "make_restricted_bundle.sh")
    out = tmp_path / "bundle"
    dry = subprocess.run(["bash", sh, "--env", str(env), "--out", str(out), "--dry-run"], capture_output=True, text=True)
    assert dry.returncode == 0 and not out.exists()
    assert "NOT READY" in dry.stdout and "EAS stability: MISSING" in dry.stdout   # mock has 3 modalities only
    # a full build refuses missing inputs ...
    assert subprocess.run(["bash", sh, "--env", str(env), "--out", str(out)], capture_output=True).returncode == 1
    # ... so give every modality a (signal-free) result table
    for anc in ("EAS", "EUR"):
        for mod in ("expression", "isoform_expression", "isoforms", "splicing", "intron_retention",
                    "alt_TSS", "alt_polyA", "RNA_editing", "stability"):
            p = q / f"{anc}_{mod}_cisqtl.parquet"
            if not p.exists():
                pd.DataFrame({"phenotype_id": ["X"], "variant_id": ["chr9:1:A:G"], "qval": [0.9]}).to_parquet(p)
    dry = subprocess.run(["bash", sh, "--env", str(env), "--out", str(out), "--dry-run"], capture_output=True, text=True)
    assert "dry run: READY" in dry.stdout
    r = subprocess.run(["bash", sh, "--env", str(env), "--out", str(out)], capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    arc = out / "xqtl-restricted-7.tar.gz"
    ok, msgs, members = R.verify(str(arc), "7")
    assert ok, msgs
    assert (out / "xqtl-restricted-7.tar.gz.sha256").read_text().split()[1] == arc.name


def test_make_public_reference_record(tmp_path):
    for t in ("plink2", "samtools"):
        if not shutil.which(t):
            pytest.skip(f"{t} not on PATH")
    a = tmp_path / "atlas"
    a.mkdir()
    for n in ("hplr.gtf", "norm.gtf", "class.txt", "rediportal.bed"):
        (a / n).write_text(f"synthetic {n}\n")
    tx = a / "txrevise"
    tx.mkdir()
    for g in ("grp_1", "grp_2"):
        for p in ("upstream", "contained", "downstream"):
            (tx / f"txrevise.{g}.{p}.gff3").write_text(f"##gff-version 3\n# {g} {p}\n")
    (a / "music.rds").write_bytes(b"synthetic rds")
    (a / "genome.fa").write_text(">chr1\nACGTACGTAC\n>chr2\nGGGCCC\n")
    # synthetic 1kGP-like panel: 80 samples x 400 variants on chr1
    subprocess.run(["plink2", "--dummy", "80", "400", "0.05", "acgt", "--seed", "1", "--make-pgen", "--out", str(a / "kg")],
                   check=True, capture_output=True)
    psam = (a / "kg.psam").read_text().splitlines()
    col = 1 if psam[0].startswith("#FID") else 0
    samples = [l.split("\t")[col] for l in psam[1:]]
    (a / "superpop.tsv").write_text("sample\tsuperpop\n" + "".join(f"{s}\t{'EUR' if i % 2 else 'EAS'}\n"
                                                                   for i, s in enumerate(samples[:60])))
    env = tmp_path / "rec.env"
    env.write_text("\n".join([
        f"ATLAS_HPLR_GTF={a}/hplr.gtf", f"ATLAS_NORMALIZED_GTF={a}/norm.gtf", f"ATLAS_CLASSIFICATION={a}/class.txt",
        f"ATLAS_TXREVISE_DIR={tx}", f"ATLAS_MUSIC_RDS={a}/music.rds", f"ATLAS_GENOME_FASTA={a}/genome.fa",
        f"KGP_PFILE={a}/kg", f"KGP_SUPERPOP={a}/superpop.tsv", "PANEL_MAF=0.01", 'PANEL_PRUNE="50 5 0.9"',
        f"ATLAS_REDIPORTAL_BED={a}/rediportal.bed", "INCLUDE_REDIPORTAL=no"]) + "\n")
    sh = os.path.join(REPO, "scripts", "make_public_reference_record.sh")
    out = tmp_path / "record"
    dry = subprocess.run(["bash", sh, "--env", str(env), "--out", str(out), "--dry-run"], capture_output=True, text=True)
    assert dry.returncode == 0 and "dry run: READY" in dry.stdout and not out.exists()
    r = subprocess.run(["bash", sh, "--env", str(env), "--out", str(out), "--threads", "2"], capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    files = pd.read_csv(out / "files.tsv", sep="\t")
    expect = {"HPLRv2.0.annotated.gtf.gz", "HPLRv2.0.annotated.PANTRY.gtf.gz", "HPLRv2_classification.annotated.txt.gz",
              "placenta_music_reference.rds", "1kGP_hg38_pca_panel.pgen", "1kGP_hg38_pca_panel.pvar.zst",
              "1kGP_hg38_pca_panel.psam", "1kg_superpop.tsv", "atlas_genome_contigs.tsv"} | {
        f"txrevise.{g}.{p}.gff3.gz" for g in ("grp_1", "grp_2") for p in ("upstream", "contained", "downstream")}
    assert set(files.name) == expect                      # no REDIportal unless INCLUDE_REDIPORTAL=yes
    assert len((out / "1kGP_hg38_pca_panel.psam").read_text().splitlines()) == 61   # header + 60 labelled
    contigs = pd.read_csv(out / "atlas_genome_contigs.tsv", sep="\t")
    assert contigs.name.tolist() == ["chr1", "chr2"] and contigs.length.tolist() == [10, 6]
    # every listed checksum is right
    for row in files.itertuples():
        assert hashlib.sha256((out / row.name).read_bytes()).hexdigest() == row.sha256

    # pin a copy of the sources table to this "record"
    src = tmp_path / "refs_sources.tsv"
    shutil.copy(os.path.join(REPO, "assets", "refs_sources.tsv"), src)
    pin = os.path.join(REPO, "scripts", "pin_reference_record.py")
    r = subprocess.run(["python3", pin, "--record", "123456", "--files", str(out / "files.tsv"), "--sources", str(src)],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    txt = src.read_text()
    assert "# ZENODO_RECORD=123456" in txt
    rows = pd.read_csv(src, sep="\t", comment="#")
    red = rows[rows.id == "rediportal_sites"].iloc[0]
    assert red.sha256 == "TO_BE_PINNED"                    # optional and not in the record
    gtf = rows[rows.id == "hplrv2_gtf"].iloc[0]
    assert gtf.sha256 == files.set_index("name").loc["HPLRv2.0.annotated.gtf.gz", "sha256"]
    # a record lacking a required file is refused
    files.iloc[1:].to_csv(tmp_path / "partial.tsv", sep="\t", index=False)
    r = subprocess.run(["python3", pin, "--record", "1", "--files", str(tmp_path / "partial.tsv"), "--sources", str(src)],
                       capture_output=True, text=True)
    assert r.returncode == 1 and "required files missing" in r.stderr
