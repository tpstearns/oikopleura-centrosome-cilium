#!/usr/bin/env python3
"""Spliced protein-to-genome alignment (miniprot) of anchors and seeds against one genome.

This is the search that lets a gene be called E3 (present in the genome, missing from the annotation).
Queries are all seed sequences (anchor + orthologs) for every gene, so a divergent Oikopleura gene has
several chances to be found.

usage: search_genomes.py --genome data/genomes/<key>.fna --seeds data/seeds --out results/miniprot/<key>.tsv \
                         --params config/params.yaml --threads 8
Output columns: gene, query, contig, start, end, strand, identity, positive, query_coverage, rank
"""
import argparse
import re
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import load_yaml, log, write_tsv  # noqa: E402
from Bio import SeqIO  # noqa: E402

FIELDS = ["gene", "query", "contig", "start", "end", "strand", "identity", "positive", "query_coverage", "rank"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", required=True)
    ap.add_argument("--seeds", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--params", required=True)
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()
    p = load_yaml(args.params)["miniprot"]
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as tmp:
        query = Path(tmp) / "queries.faa"
        qlen = {}
        with open(query, "w") as fh:
            for seedfile in sorted(Path(args.seeds).glob("*.faa")):
                for rec in SeqIO.parse(seedfile, "fasta"):
                    qlen[rec.id] = len(rec.seq)
                    SeqIO.write(rec, fh, "fasta")
        idx = Path(tmp) / "genome.mpi"
        subprocess.run(["miniprot", "-t", str(args.threads), "-d", str(idx), args.genome], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        gff = Path(tmp) / "aln.gff"
        with open(gff, "w") as fh:
            subprocess.run(["miniprot", "-t", str(args.threads), "--gff", "-G", str(p["max_intron_bp"]),
                            "-N", str(p["outn"]), str(idx), str(query)], check=True, stdout=fh,
                           stderr=subprocess.DEVNULL)
        rows = []
        attr_re = re.compile(r"([A-Za-z_]+)=([^;]+)")
        for line in open(gff):
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "mRNA":
                continue
            attrs = dict(attr_re.findall(f[8]))
            tgt = attrs.get("Target", "").split()
            if len(tgt) != 3:
                continue
            q, tstart, tend = tgt[0], int(tgt[1]), int(tgt[2])
            cov = (tend - tstart + 1) / qlen[q] if q in qlen else 0.0
            rows.append({"gene": q.split("|")[0], "query": q, "contig": f[0], "start": f[3], "end": f[4],
                         "strand": f[6], "identity": attrs.get("Identity", ""), "positive": attrs.get("Positive", ""),
                         "query_coverage": f"{cov:.3f}", "rank": attrs.get("Rank", "")})
    write_tsv(out, rows, FIELDS)
    log(f"{len(rows)} spliced alignments -> {out}")


if __name__ == "__main__":
    main()
