#!/usr/bin/env python3
"""Reciprocal best hits between the anchor set and one proteome (MMseqs2 easy-rbh).

usage: rbh.py --anchors data/genes/anchors.faa --proteome data/proteomes/<key>.faa --out results/rbh/<key>.tsv \
              --params config/params.yaml --threads 8
Output columns: gene, anchor, target, pident, evalue, bits, query_coverage, target_coverage
"""
import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import load_yaml, log, write_tsv  # noqa: E402

FIELDS = ["gene", "anchor", "target", "pident", "evalue", "bits", "query_coverage", "target_coverage"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--anchors", required=True)
    ap.add_argument("--proteome", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--params", required=True)
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()
    p = load_yaml(args.params)["rbh"]
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    raw = out.with_suffix(".raw.tsv")
    with tempfile.TemporaryDirectory() as tmp:
        cmd = ["mmseqs", "easy-rbh", args.anchors, args.proteome, str(raw), tmp, "--threads", str(args.threads),
               "-s", str(p["sensitivity"]), "-e", "10", "--format-output",
               "query,target,pident,alnlen,qstart,qend,qlen,tstart,tend,tlen,evalue,bits"]
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    rows = []
    for line in open(raw):
        q, t, pident, alnlen, qstart, qend, qlen, tstart, tend, tlen, evalue, bits = line.rstrip("\n").split("\t")
        rows.append({"gene": q.split("|")[0], "anchor": q, "target": t, "pident": pident, "evalue": evalue,
                     "bits": bits, "query_coverage": f"{(int(qend) - int(qstart) + 1) / int(qlen):.3f}",
                     "target_coverage": f"{(int(tend) - int(tstart) + 1) / int(tlen):.3f}"})
    raw.unlink()
    write_tsv(out, rows, FIELDS)
    log(f"{len(rows)} reciprocal best hits -> {out}")


if __name__ == "__main__":
    main()
