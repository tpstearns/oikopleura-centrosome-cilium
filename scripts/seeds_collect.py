#!/usr/bin/env python3
"""Collect ortholog seed sequences per gene from the seed species via reciprocal best hits.

For each seed species proteome, run MMseqs2 easy-rbh (anchors vs proteome); for each gene with an RBH
passing thresholds, add the orthologous sequence to data/seeds/<gene>.faa (anchor sequence always included).

usage: seeds_collect.py --anchors data/genes/anchors.faa --params config/params.yaml --outdir data/seeds \
                        --threads 8 --proteomes Hsap=data/proteomes/Hsap.faa Mmus=data/proteomes/Mmus.faa ...
"""
import argparse
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import die, load_yaml, log  # noqa: E402
from Bio import SeqIO  # noqa: E402


def run_rbh(query: Path, target: Path, out: Path, threads: int, sensitivity: float) -> list[list[str]]:
    with tempfile.TemporaryDirectory() as tmp:
        cmd = ["mmseqs", "easy-rbh", str(query), str(target), str(out), tmp, "--threads", str(threads),
               "-s", str(sensitivity), "--format-output",
               "query,target,pident,alnlen,qstart,qend,qlen,tstart,tend,tlen,evalue,bits"]
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return [line.rstrip("\n").split("\t") for line in open(out)]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--anchors", required=True)
    ap.add_argument("--params", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--proteomes", nargs="+", required=True, help="KEY=path.faa ...")
    args = ap.parse_args()
    p = load_yaml(args.params)["rbh"]
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    anchors = {rec.id: rec for rec in SeqIO.parse(args.anchors, "fasta")}
    seeds: dict[str, list] = defaultdict(list)
    for aid, rec in anchors.items():
        gene = aid.split("|")[0]
        rec.id, rec.description = f"{gene}|anchor|{aid.split('|')[1]}", ""
        seeds[gene].append(rec)

    for item in args.proteomes:
        key, path = item.split("=", 1)
        prot = Path(path)
        tmp_out = outdir / f"_rbh_{key}.tsv"
        rows = run_rbh(Path(args.anchors), prot, tmp_out, args.threads, p["sensitivity"])
        index = SeqIO.index(str(prot), "fasta")
        n = 0
        for q, t, pident, alnlen, qstart, qend, qlen, tstart, tend, tlen, evalue, bits in rows:
            if float(evalue) > p["evalue_max"]:
                continue
            if (int(qend) - int(qstart) + 1) / int(qlen) < p["min_query_coverage"]:
                continue
            gene = q.split("|")[0]
            rec = index[t]
            rec.id, rec.description = f"{gene}|{key}|{t.split('|')[1] if '|' in t else t}", ""
            seeds[gene].append(rec)
            n += 1
        tmp_out.unlink()
        log(f"[{key}] {n} RBH seeds added")

    for gene, recs in seeds.items():
        SeqIO.write(recs, outdir / f"{gene}.faa", "fasta")
    log(f"wrote seeds for {len(seeds)} genes to {outdir}")


if __name__ == "__main__":
    main()
