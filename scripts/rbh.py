#!/usr/bin/env python3
"""Reciprocal best hits between the anchor genes and one proteome, reciprocal against the FULL anchor-species proteome.

For each anchor species (human for most genes, Xenopus for TUBZ1, ...), MMseqs2 easy-rbh runs the whole reference
proteome of that species against the target proteome; a target protein is the RBH of a listed gene only when its
best hit in the anchor species is that gene's anchor. Paralogs that are not in the gene list therefore cannot be
reported as the listed gene.

usage: rbh.py --anchors data/genes/anchors.faa --anchors-tsv data/genes/anchors.tsv \
              --anchor-proteomes Hsap=data/proteomes/Hsap.faa Xtro=data/proteomes/Xtro.faa \
              --proteome data/proteomes/<key>.faa --out results/rbh/<key>.tsv --params config/params.yaml --threads 8
Output columns: gene, anchor, target, pident, evalue, bits, query_coverage, target_coverage
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import anchor_table, die, load_yaml, log, rbh_full_proteome, write_tsv  # noqa: E402

FIELDS = ["gene", "anchor", "target", "pident", "evalue", "bits", "query_coverage", "target_coverage"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--anchors", required=True)
    ap.add_argument("--anchors-tsv", required=True)
    ap.add_argument("--anchor-proteomes", nargs="+", required=True, help="KEY=path.faa for every anchor species")
    ap.add_argument("--proteome", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--params", required=True)
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()
    p = load_yaml(args.params)["rbh"]
    prot_by_species = dict(item.split("=", 1) for item in args.anchor_proteomes)
    anchors = anchor_table(args.anchors, args.anchors_tsv)
    rows = []
    for sp in sorted({a["species"] for a in anchors.values()}):
        if sp not in prot_by_species:
            die(f"no proteome given for anchor species {sp}")
        wanted = {a["accession"]: g for g, a in anchors.items() if a["species"] == sp}
        hits = rbh_full_proteome(prot_by_species[sp], args.proteome, wanted, args.threads, p["sensitivity"], evalue=10.0)
        for h in hits:
            rows.append({"gene": h["gene"], "anchor": h["anchor"], "target": h["target"], "pident": h["pident"],
                         "evalue": h["evalue"], "bits": h["bits"],
                         "query_coverage": f"{(int(h['qend']) - int(h['qstart']) + 1) / int(h['qlen']):.3f}",
                         "target_coverage": f"{(int(h['tend']) - int(h['tstart']) + 1) / int(h['tlen']):.3f}"})
        log(f"  [{sp} proteome as reciprocal set] {len(hits)} of {len(wanted)} anchors have an RBH")
    write_tsv(args.out, rows, FIELDS)
    log(f"{len(rows)} reciprocal best hits -> {args.out}")


if __name__ == "__main__":
    main()
