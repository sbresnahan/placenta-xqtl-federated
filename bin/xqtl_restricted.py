#!/usr/bin/env python3
"""Create, verify and unpack the restricted file xqtl-restricted-<version>.tar.gz.

The restricted file is sent by MDACC to approved sites by email. It contains:
  RESTRICTED_VERSION                     format/content version (one line)
  MANIFEST.sha256                        "<sha256>  <member path>" for every other member
  splicing/frozen_clusters.tsv.gz        shared leafCutter reference clusters
                                         (chrom, start, end, strand, cluster_id); junction
                                         coordinates only, no per-sample data
  signals/mdacc_signals.tsv.gz           MDACC xQTL signals tested at every site for G x C and
                                         GxE tier 2 (phenotype_id, variant_id, modality,
                                         ancestry, signal_rank)
  README.txt

Subcommands
  make     build the archive from the two tables (MDACC only; scripts/make_restricted_bundle.sh)
  verify   check version, members and checksums; exit 1 with instructions if anything is wrong
  extract  verify, then unpack to a directory
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import os
import sys
import tarfile

REQUIRED = {
    "splicing/frozen_clusters.tsv.gz": ["chrom", "start", "end", "strand", "cluster_id"],
    "signals/mdacc_signals.tsv.gz": ["phenotype_id", "variant_id", "modality", "ancestry", "signal_rank"],
}
HINT = (
    "If your site has not received the restricted file, or it fails verification, request it by "
    "email from the scientific contact (Sean T. Bresnahan, stbresnahan@mdanderson.org), stating "
    "your cohort, the pipeline version and the approved data-sharing agreement. See docs/SETUP.md."
)


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def read_members(path: str) -> dict[str, bytes]:
    out = {}
    with tarfile.open(path, "r:gz") as tf:
        for m in tf.getmembers():
            if not m.isfile():
                continue
            name = m.name[2:] if m.name.startswith("./") else m.name
            if name.startswith("/") or ".." in name.split("/"):
                raise ValueError(f"unsafe member path in archive: {m.name}")
            out[name] = tf.extractfile(m).read()
    return out


def verify(path: str, expected_version: str | None) -> tuple[bool, list[str], dict[str, bytes]]:
    msgs: list[str] = []
    if not os.path.exists(path):
        return False, [f"restricted file not found: {path}"], {}
    try:
        members = read_members(path)
    except (tarfile.TarError, OSError, ValueError) as e:
        return False, [f"cannot read {path} as a .tar.gz archive: {e}"], {}
    for need in ["RESTRICTED_VERSION", "MANIFEST.sha256", *REQUIRED]:
        if need not in members:
            msgs.append(f"missing member: {need}")
    if msgs:
        return False, msgs, members
    version = members["RESTRICTED_VERSION"].decode().strip()
    if expected_version and version != str(expected_version):
        msgs.append(f"restricted file version {version} does not match the version this pipeline "
                    f"release expects ({expected_version})")
    for line in members["MANIFEST.sha256"].decode().splitlines():
        if not line.strip():
            continue
        digest, name = line.split(None, 1)
        name = name.strip()
        if name not in members:
            msgs.append(f"manifest lists {name} but it is not in the archive")
        elif sha256_bytes(members[name]) != digest:
            msgs.append(f"checksum mismatch for {name}")
    for name, cols in REQUIRED.items():
        try:
            header = gzip.decompress(members[name]).decode().split("\n", 1)[0].split("\t")
        except OSError as e:
            msgs.append(f"{name} is not gzip-compressed text: {e}")
            continue
        if header[: len(cols)] != cols:
            msgs.append(f"{name} header {header} does not start with {cols}")
    return not msgs, msgs, members


def cmd_verify(a) -> int:
    ok, msgs, members = verify(a.archive, a.expected_version)
    if not ok:
        for m in msgs:
            print(f"[FAIL] {m}", file=sys.stderr)
        print("\n" + HINT, file=sys.stderr)
        return 1
    v = members["RESTRICTED_VERSION"].decode().strip()
    n_cl = gzip.decompress(members["splicing/frozen_clusters.tsv.gz"]).count(b"\n") - 1
    n_sig = gzip.decompress(members["signals/mdacc_signals.tsv.gz"]).count(b"\n") - 1
    print(f"restricted file OK: version {v}; {n_cl} frozen splicing junctions; {n_sig} MDACC signals")
    return 0


def cmd_extract(a) -> int:
    ok, msgs, members = verify(a.archive, a.expected_version)
    if not ok:
        for m in msgs:
            print(f"[FAIL] {m}", file=sys.stderr)
        print("\n" + HINT, file=sys.stderr)
        return 1
    for name, data in members.items():
        dest = os.path.join(a.outdir, name)
        os.makedirs(os.path.dirname(dest) or ".", exist_ok=True)
        with open(dest, "wb") as fh:
            fh.write(data)
    print(f"extracted {len(members)} files to {a.outdir}")
    return 0


def cmd_make(a) -> int:
    payload = {}
    for name, src in (("splicing/frozen_clusters.tsv.gz", a.clusters), ("signals/mdacc_signals.tsv.gz", a.signals)):
        with open(src, "rb") as fh:
            raw = fh.read()
        if not src.endswith(".gz"):
            raw = gzip.compress(raw, mtime=0)
        header = gzip.decompress(raw).decode().split("\n", 1)[0].split("\t")
        if header[: len(REQUIRED[name])] != REQUIRED[name]:
            print(f"{src}: header {header} must start with {REQUIRED[name]}", file=sys.stderr)
            return 1
        payload[name] = raw
    payload["README.txt"] = (
        f"placenta-xqtl-federated restricted file, version {a.version}\n"
        "Prepared by MD Anderson Cancer Center (Bhattacharya Lab).\n"
        "For use by approved participating cohorts only, under their data-sharing agreement.\n"
        "Do not redistribute. Place the archive at the path given by 'restricted' in cohort.yml.\n"
    ).encode()
    manifest = "".join(f"{sha256_bytes(d)}  {n}\n" for n, d in sorted(payload.items()))
    payload["MANIFEST.sha256"] = manifest.encode()
    payload["RESTRICTED_VERSION"] = f"{a.version}\n".encode()
    # Byte-reproducible: fixed member order/mtime/owner, and the gzip layer is written
    # separately with mtime 0 and no embedded file name.
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w", format=tarfile.PAX_FORMAT) as tf:
        for name in sorted(payload):
            info = tarfile.TarInfo(name)
            info.size = len(payload[name])
            info.mtime = 0
            info.mode = 0o644
            info.uid = info.gid = 0
            info.uname = info.gname = ""
            tf.addfile(info, io.BytesIO(payload[name]))
    with open(a.out, "wb") as fh:
        fh.write(gzip.compress(buf.getvalue(), compresslevel=9, mtime=0))
    print(f"wrote {a.out} (sha256 {hashlib.sha256(open(a.out, 'rb').read()).hexdigest()})")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    v = sub.add_parser("verify")
    v.add_argument("archive")
    v.add_argument("--expected-version", default=None)
    v.set_defaults(func=cmd_verify)
    e = sub.add_parser("extract")
    e.add_argument("archive")
    e.add_argument("outdir")
    e.add_argument("--expected-version", default=None)
    e.set_defaults(func=cmd_extract)
    m = sub.add_parser("make")
    m.add_argument("--clusters", required=True)
    m.add_argument("--signals", required=True)
    m.add_argument("--version", required=True)
    m.add_argument("--out", required=True)
    m.set_defaults(func=cmd_make)
    a = ap.parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
