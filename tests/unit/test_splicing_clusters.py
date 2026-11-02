"""Unit tests for bin/xqtl_splicing_clusters.py and the atlas splicing-ID conversion."""
import gzip
import os
import sys

import pandas as pd

from conftest import REPO, load

S = load("xqtl_splicing_clusters")
sys.path.insert(0, os.path.join(REPO, "scripts", "lib"))
import atlas_restricted_tables as A  # noqa: E402


def wgz(p, text):
    with gzip.open(p, "wt") as fh:
        fh.write(text)
    return str(p)


def test_frozen_refined_and_merge(tmp_path):
    fz = wgz(tmp_path / "fz.tsv.gz", "chrom\tstart\tend\tstrand\tcluster_id\n"
             "chr1\t100\t200\t+\tclu_7_+\nchr1\t100\t300\t+\tclu_7_+\nchr2\t50\t90\t-\tclu_8_-\n")
    S.main(["frozen-refined", "--frozen", fz, "--chroms", "chr1", "--out", str(tmp_path / "ref.txt"),
            "--map", str(tmp_path / "map.tsv")])
    assert (tmp_path / "ref.txt").read_text() == "chr1:+ 100:200:0 100:300:0\n"
    # leafCutter numbers the refined clusters by line: clu_1_+ is frozen clu_7_+
    fn = wgz(tmp_path / "f.gz", "A B\nchr1:100:200:clu_1_+ 3 4\nchr1:100:300:clu_1_+ 1 0\n")
    loc = wgz(tmp_path / "l.gz", "A B\nchr1:100:200:clu_1_+ 3 4\nchr1:100:250:clu_1_+ 2 2\n"
              "chr1:500:600:clu_2_+ 9 9\nchr1:500:700:clu_2_+ 1 1\n")
    S.main(["merge", "--frozen-numers", fn, "--frozen-map", str(tmp_path / "map.tsv"), "--local-numers", loc,
            "--out", str(tmp_path / "o.gz"), "--summary", str(tmp_path / "s.tsv")])
    rows = gzip.open(tmp_path / "o.gz", "rt").read().splitlines()
    ids = [r.split()[0] for r in rows[1:]]
    # local cluster 1 shares a junction with the frozen set -> dropped; cluster 2 is novel
    assert ids == ["chr1:100:200:clu_7_+", "chr1:100:300:clu_7_+", "chr1:500:600:nclu_2_+", "chr1:500:700:nclu_2_+"]
    summ = pd.read_csv(tmp_path / "s.tsv", sep="\t").set_index("set")
    assert summ.loc["novel", "clusters"] == 1 and summ.loc["local_matched_to_frozen", "clusters"] == 1


def test_merge_with_no_local_clusters(tmp_path):
    fn = wgz(tmp_path / "f.gz", "A B\nchr1:100:200:clu_1_+ 3 4\n")
    (tmp_path / "map.tsv").write_text("line\tcluster_id\n1\tclu_7_+\n")
    loc = wgz(tmp_path / "l.gz", "A B\n")
    S.main(["merge", "--frozen-numers", fn, "--frozen-map", str(tmp_path / "map.tsv"), "--local-numers", loc,
            "--out", str(tmp_path / "o.gz"), "--summary", str(tmp_path / "s.tsv")])
    assert gzip.open(tmp_path / "o.gz", "rt").read().splitlines()[1:] == ["chr1:100:200:clu_7_+ 3 4"]


def test_atlas_splicing_ids_rewritten_to_frozen_clusters():
    jmap = {("chr22", 100, 200, "+"): "clu_7_+"}
    assert A.canonical_splicing("ENSG1.2__chr22_100_200_clu_33_+", jmap) == "ENSG1.2__chr22_100_200_clu_7_+"
    assert A.canonical_splicing("chr22:100:200:clu_3_+", jmap) == "chr22_100_200_clu_7_+"
