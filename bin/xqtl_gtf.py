"""xqtl_gtf.py -- minimal GTF reader (replaces gtfparse.read_gtf for the ported atlas/PANTRY scripts).

read_gtf(path, features=None) returns a pandas DataFrame with the columns gtfparse produces
(seqname, source, feature, start, end, score, strand, frame) plus one column per attribute
key (first occurrence; repeated keys such as `tag` are joined with ","). Plain or gzipped input.
"""
from __future__ import annotations

import gzip
import re

import pandas as pd

_ATTR = re.compile(r'\s*([^\s;]+)\s+(?:"([^"]*)"|([^;\s]+))\s*;?')
COLS = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame"]


def _open(path):
    path = str(path)
    with open(path, "rb") as fh:
        magic = fh.read(2)
    return gzip.open(path, "rt") if magic == b"\x1f\x8b" else open(path)


def parse_attributes(s: str) -> dict:
    out: dict = {}
    for m in _ATTR.finditer(s):
        k = m.group(1)
        v = m.group(2) if m.group(2) is not None else m.group(3)
        out[k] = f"{out[k]},{v}" if k in out else v
    return out


def read_gtf(path, features=None) -> pd.DataFrame:
    feats = set(features) if features else None
    rows = []
    with _open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or (feats is not None and c[2] not in feats):
                continue
            rec = dict(zip(COLS, c[:8]))
            rec.update(parse_attributes(c[8]))
            rows.append(rec)
    df = pd.DataFrame(rows)
    if df.empty:
        return pd.DataFrame(columns=COLS)
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    return df
