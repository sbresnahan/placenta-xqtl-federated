"""Unit tests for bin/xqtl_refs_check.py (embedded reference manifest check)."""
import hashlib
import shutil

from conftest import FIX, load

C = load("xqtl_refs_check")


def test_fixture_ok():
    assert C.main(["--refs-dir", f"{FIX}/refs_mini"]) == 0
    assert C.main(["--refs-dir", f"{FIX}/refs_mini", "--quick"]) == 0


def test_corruption_detected(tmp_path):
    d = tmp_path / "refs"
    shutil.copytree(f"{FIX}/refs_mini", d)
    p = d / "README.txt"
    b = p.read_bytes()
    p.write_bytes(b[:-2] + b"X\n")       # same size, different content
    assert C.main(["--refs-dir", str(d), "--quick"]) == 0   # quick = size only
    assert C.main(["--refs-dir", str(d)]) == 1


def test_missing_manifest(tmp_path, capsys):
    assert C.main(["--refs-dir", str(tmp_path)]) == 1
    assert "MANIFEST.tsv not found" in capsys.readouterr().out


def test_expected_manifest_mismatch(tmp_path):
    exp = tmp_path / "expected.tsv"
    exp.write_text("path\tsha256\tbytes\nREADME.txt\t" + hashlib.sha256(b"other").hexdigest() + "\t5\n")
    assert C.main(["--refs-dir", f"{FIX}/refs_mini", "--quick", "--expected", str(exp)]) == 1
