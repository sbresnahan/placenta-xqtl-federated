#!/usr/bin/env python3
"""Check the reference data embedded in the image (/opt/xqtl/refs) against its manifest.

refs/MANIFEST.tsv columns: path, sha256, bytes, description, source, license
  --quick      existence and byte size only (seconds; used by --mode smoke)
  (default)    full sha256 of every file (minutes; used at the start of site runs)
  --expected   a second manifest shipped with the pipeline code for this release
               (assets/refs_manifest.tsv); every row must match refs/MANIFEST.tsv exactly,
               which catches a pipeline run against the wrong image.
The last line of stdout is a one-line summary. Exit 0 = all OK, 1 = problems.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import os
import sys

COLS = ["path", "sha256", "bytes", "description", "source", "license"]


def read_manifest(path: str, allow_empty: bool = False) -> list[dict]:
    with open(path, newline="") as fh:
        reader = csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t")
        rows = list(reader)
        cols = reader.fieldnames or []
    if [c for c in COLS[:3] if c not in cols] or (not rows and not allow_empty):
        raise ValueError(f"{path}: needs columns {COLS[:3]} and at least one row")
    return rows


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refs-dir", required=True)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--expected", default=None)
    a = ap.parse_args(argv)

    man = os.path.join(a.refs_dir, "MANIFEST.tsv")
    if not os.path.isfile(man):
        print(f"FAIL\t{man} not found: is this the placenta-xqtl-federated image (or is --refs_dir wrong)?")
        return 1
    rows = read_manifest(man)
    problems = []
    for r in rows:
        p = os.path.join(a.refs_dir, r["path"])
        if not os.path.isfile(p):
            problems.append(f"missing {r['path']}")
            continue
        if str(os.path.getsize(p)) != str(r["bytes"]):
            problems.append(f"size {r['path']}: {os.path.getsize(p)} != {r['bytes']}")
            continue
        if not a.quick and sha256(p) != r["sha256"]:
            problems.append(f"sha256 {r['path']}")
    if a.expected:
        exp = {r["path"]: r for r in read_manifest(a.expected, allow_empty=True)}
        if not exp:
            print(f"NOTE\t{a.expected} has no entries (development build): release pins not checked")
        got = {r["path"]: r for r in rows}
        for k, r in exp.items():
            if k not in got:
                problems.append(f"expected-not-in-image {k}")
            elif got[k]["sha256"] != r["sha256"]:
                problems.append(f"image-differs-from-release {k}")
    for pr in problems:
        print(f"PROBLEM\t{pr}")
    level = "quick (size)" if a.quick else "full (sha256)"
    if problems:
        print(f"FAIL\t{len(problems)} problem(s) in {len(rows)} reference files [{level}]")
        return 1
    print(f"OK\t{len(rows)} reference files verified [{level}]")
    return 0


if __name__ == "__main__":
    sys.exit(main())
