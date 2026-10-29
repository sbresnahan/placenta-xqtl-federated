#!/usr/bin/env python3
"""Pin the public MDACC Zenodo reference record in assets/refs_sources.tsv.

After uploading the files written by scripts/make_public_reference_record.sh to Zenodo:

    scripts/pin_reference_record.py --record 1234567 --files /path/record/files.tsv

sets the '# ZENODO_RECORD=' line and the sha256/bytes of every Zenodo row whose file is in
files.tsv (name, sha256, bytes). Required rows missing from the record are an error;
optional rows (REDIportal sites, the genome provenance table) stay unpinned and are skipped
by the image build. Then commit assets/refs_sources.tsv.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--record", required=True, help="Zenodo record id (digits)")
    ap.add_argument("--files", required=True, help="files.tsv from make_public_reference_record.sh")
    ap.add_argument("--sources", default=str(HERE / "assets" / "refs_sources.tsv"))
    a = ap.parse_args(argv)
    if not re.fullmatch(r"\d+", a.record):
        print("--record must be the numeric Zenodo record id", file=sys.stderr)
        return 2

    with open(a.files, newline="") as fh:
        files = {r["name"]: r for r in csv.DictReader(fh, delimiter="\t")}

    lines = Path(a.sources).read_text().splitlines()
    out, missing, pinned = [], [], 0
    header = None
    for line in lines:
        if line.startswith("#"):
            out.append(re.sub(r"^# *ZENODO_RECORD=\S*", f"# ZENODO_RECORD={a.record}", line))
            continue
        f = line.split("\t")
        if header is None:
            header = f
            out.append(line)
            continue
        r = dict(zip(header, f))
        if "{ZENODO_RECORD}" in r["url"]:
            name = r["url"].split("/files/", 1)[1].split("?", 1)[0]
            if name in files:
                r["sha256"], r["bytes"] = files[name]["sha256"], files[name]["bytes"]
                pinned += 1
            else:
                r["sha256"] = r["bytes"] = "TO_BE_PINNED"
                if r["required"] == "yes":
                    missing.append(name)
        out.append("\t".join(r[h] for h in header))
    if missing:
        print("required files missing from the record: " + ", ".join(missing), file=sys.stderr)
        return 1
    Path(a.sources).write_text("\n".join(out) + "\n")
    print(f"pinned {pinned} file(s) from Zenodo record {a.record} in {a.sources}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
