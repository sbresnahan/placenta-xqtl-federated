"""Tests for scripts/fetch_references.sh with local file:// sources (synthetic data).

The build step needs bgzip + samtools (the image and the dev environment have both).
"""
import gzip
import hashlib
import os
import shutil
import subprocess

import pytest

from conftest import REPO, load

C = load("xqtl_refs_check")
SH = os.path.join(REPO, "scripts", "fetch_references.sh")
HDR = "id\trequired\turl\tsha256\tbytes\taction\tdest\tdescription\tsource\tlicense\n"


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def make_sources(tmp_path, extra_rows=""):
    src = tmp_path / "src"
    src.mkdir()
    with gzip.open(src / "genome.fna.gz", "wt") as fh:
        fh.write(">chr1\n" + "ACGT" * 50 + "\n>chr2\n" + "GGCC" * 30 + "\n")
    (src / "chain.txt").write_text("chain 1 chr1 10 + 0 10 chr1 10 + 0 10 1\n10\n")
    g, c = src / "genome.fna.gz", src / "chain.txt"
    rows = (f"genome\tyes\tfile://{g}\t{sha(g)}\t{g.stat().st_size}\tbgzip_fasta\tgenome/g.fa.gz\ttest genome\tsynthetic\tCC0\n"
            f"chain\tyes\tfile://{c}\t{sha(c)}\t{c.stat().st_size}\tcopy\tliftover/chain.txt\ttest chain\tsynthetic\tCC0\n"
            + extra_rows)
    t = tmp_path / "sources.tsv"
    t.write_text("# test sources\n# ZENODO_RECORD=TO_BE_PINNED\n" + HDR + rows)
    return t


def run(*args):
    return subprocess.run(["bash", SH, *map(str, args)], capture_output=True, text=True)


@pytest.fixture(autouse=True)
def tools():
    for t in ("bgzip", "samtools", "curl"):
        if not shutil.which(t):
            pytest.skip(f"{t} not on PATH")


def test_download_build_and_manifest(tmp_path):
    t = make_sources(tmp_path, "zen\toptional\thttps://zenodo.org/records/{ZENODO_RECORD}/files/x.bed.gz?download=1"
                               "\tTO_BE_PINNED\tTO_BE_PINNED\tcopy\tediting/x.bed.gz\topt\tz\tz\n")
    cache, out = tmp_path / "cache", tmp_path / "refs"
    r = run("download", "--sources", t, "--cache", cache)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "skip  zen" in r.stdout                        # optional + unpinned -> skipped, not fatal
    r = run("build", "--sources", t, "--cache", cache, "--out", out, "--threads", "1")
    assert r.returncode == 0, r.stdout + r.stderr
    paths = [l.split("\t")[0] for l in (out / "MANIFEST.tsv").read_text().splitlines()[1:]]
    assert paths == ["genome/g.fa.gz", "genome/g.fa.gz.fai", "genome/g.fa.gz.gzi", "liftover/chain.txt"]
    assert C.main(["--refs-dir", str(out)]) == 0
    # rebuilding gives byte-identical outputs (stable release manifest)
    m1 = (out / "MANIFEST.tsv").read_bytes()
    shutil.rmtree(out)
    run("build", "--sources", t, "--cache", cache, "--out", out, "--threads", "1")
    assert (out / "MANIFEST.tsv").read_bytes() == m1
    # expected-manifest comparison: strict fails on a difference, report only warns
    exp = tmp_path / "expected.tsv"
    exp.write_bytes(m1.replace(b"test chain", b"other chain"))
    assert run("build", "--sources", t, "--cache", cache, "--out", out, "--expected", exp,
               "--expected-mode", "strict").returncode == 1
    assert run("build", "--sources", t, "--cache", cache, "--out", out, "--expected", exp).returncode == 0
    exp.write_bytes(m1)
    r = run("build", "--sources", t, "--cache", cache, "--out", out, "--expected", exp, "--expected-mode", "strict")
    assert r.returncode == 0 and "matches the release manifest" in r.stdout


def test_checksum_mismatch_is_fatal(tmp_path):
    t = make_sources(tmp_path)
    t.write_text(t.read_text().replace("\tcopy\t", "\tcopy\t").replace(sha(tmp_path / "src" / "chain.txt"), "0" * 64))
    r = run("download", "--sources", t, "--cache", tmp_path / "cache")
    assert r.returncode == 1 and "checksum mismatch for chain" in r.stderr
    assert not list((tmp_path / "cache").glob("chain.txt*"))


def test_unpinned_required_rows(tmp_path):
    t = make_sources(tmp_path, "gtf\tyes\thttps://zenodo.org/records/{ZENODO_RECORD}/files/a.gtf.gz?download=1"
                               "\tTO_BE_PINNED\tTO_BE_PINNED\tcopy\tannotation/a.gtf.gz\tgtf\tz\tz\n")
    r = run("download", "--sources", t, "--cache", tmp_path / "cache")
    assert r.returncode == 1 and "not pinned yet: gtf" in r.stderr and "pin_reference_record.py" in r.stderr
    r = run("download", "--sources", t, "--cache", tmp_path / "cache", "--allow-unpinned")
    assert r.returncode == 0 and "development build only" in r.stdout


def test_real_sources_table_is_well_formed():
    rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(REPO, "assets", "refs_sources.tsv")) if not l.startswith("#")]
    assert rows[0] == HDR.rstrip("\n").split("\t")
    for r in rows[1:]:
        assert len(r) == 10, r
        assert r[1] in ("yes", "optional") and r[5] in ("copy", "bgzip_fasta")
        assert r[3] == "TO_BE_PINNED" or len(r[3]) == 64
    dests = [r[6] for r in rows[1:]]
    assert len(dests) == len(set(dests))
