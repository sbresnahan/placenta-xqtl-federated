"""Unit tests for bin/xqtl_restricted.py (restricted file: make / verify / extract)."""
import gzip
import io
import os
import tarfile

import pytest

from conftest import load

R = load("xqtl_restricted")


def gz(tmp_path, name, text):
    p = tmp_path / name
    p.write_bytes(gzip.compress(text.encode(), mtime=0))
    return str(p)


@pytest.fixture
def inputs(tmp_path):
    clu = gz(tmp_path, "c.tsv.gz", "chrom\tstart\tend\tstrand\tcluster_id\nchr1\t100\t200\t+\tclu_1_+\n")
    sig = gz(tmp_path, "s.tsv.gz",
             "phenotype_id\tvariant_id\tmodality\tancestry\tsignal_rank\nENSG1\tchr1:150:A:G\texpression\tEAS\t1\n")
    return clu, sig


def make(tmp_path, inputs, version="1", name="r.tar.gz"):
    out = str(tmp_path / name)
    assert R.main(["make", "--clusters", inputs[0], "--signals", inputs[1], "--version", version, "--out", out]) == 0
    return out


def test_roundtrip_and_extract(tmp_path, inputs):
    out = make(tmp_path, inputs)
    ok, msgs, members = R.verify(out, "1")
    assert ok, msgs
    assert R.main(["extract", out, str(tmp_path / "x"), "--expected-version", "1"]) == 0
    assert os.path.isfile(tmp_path / "x" / "signals" / "mdacc_signals.tsv.gz")
    assert os.path.isfile(tmp_path / "x" / "splicing" / "frozen_clusters.tsv.gz")


def test_deterministic(tmp_path, inputs):
    a = open(make(tmp_path, inputs, name="a.tar.gz"), "rb").read()
    b = open(make(tmp_path, inputs, name="b.tar.gz"), "rb").read()
    assert a == b


def test_version_mismatch(tmp_path, inputs):
    ok, msgs, _ = R.verify(make(tmp_path, inputs, version="1"), "2")
    assert not ok and any("does not match" in m for m in msgs)


def test_missing_file_message_has_email(tmp_path, capsys):
    assert R.main(["verify", str(tmp_path / "nope.tar.gz")]) != 0
    assert "stbresnahan@mdanderson.org" in capsys.readouterr().err


def rewrite(src, dst, mutate):
    with tarfile.open(src, "r:gz") as t:
        members = {m.name: t.extractfile(m).read() for m in t.getmembers() if m.isfile()}
    mutate(members)
    with tarfile.open(dst, "w:gz") as t:
        for n, b in members.items():
            ti = tarfile.TarInfo(n)
            ti.size = len(b)
            t.addfile(ti, io.BytesIO(b))


def test_tampered_member_detected(tmp_path, inputs):
    out = make(tmp_path, inputs)
    bad = str(tmp_path / "bad.tar.gz")

    def m(mem):
        mem["signals/mdacc_signals.tsv.gz"] = gzip.compress(
            b"phenotype_id\tvariant_id\tmodality\tancestry\tsignal_rank\nENSG2\tchr1:1:A:G\texpression\tEUR\t1\n")
    rewrite(out, bad, m)
    ok, msgs, _ = R.verify(bad, "1")
    assert not ok and any("checksum mismatch" in x for x in msgs)


def test_unsafe_path_rejected(tmp_path, inputs):
    out = make(tmp_path, inputs)
    bad = str(tmp_path / "evil.tar.gz")
    rewrite(out, bad, lambda mem: mem.__setitem__("../evil.txt", b"x"))
    ok, msgs, _ = R.verify(bad, "1")
    assert not ok


def test_bad_header_rejected(tmp_path, inputs):
    sig = gz(tmp_path, "bad.tsv.gz", "gene\tsnp\nENSG1\tchr1:1:A:G\n")
    out = str(tmp_path / "r.tar.gz")
    rc = R.main(["make", "--clusters", inputs[0], "--signals", sig, "--version", "1", "--out", out])
    if rc == 0:
        ok, msgs, _ = R.verify(out, "1")
        assert not ok
    else:
        assert rc != 0
