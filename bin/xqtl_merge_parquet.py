#!/usr/bin/env python3
"""xqtl_merge_parquet.py OUT.parquet IN.parquet...  -- concatenate parquet tables (e.g. independent-signal
chunks) in the order given; empty inputs are skipped and an empty table is written if all are empty."""
import sys

import pandas as pd


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) < 2:
        sys.exit(__doc__)
    out, ins = argv[0], argv[1:]
    parts = [d for d in (pd.read_parquet(p) for p in ins) if len(d)]
    df = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=["phenotype_id", "variant_id"])
    df.to_parquet(out, index=False)
    print(f"{out}: {len(df)} rows from {len(ins)} file(s)")


if __name__ == "__main__":
    main()
