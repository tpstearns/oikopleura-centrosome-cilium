#!/usr/bin/env python3
"""Collect ortholog seed sequences per gene from the seed species via reciprocal best hits.

For each seed species proteome, run MMseqs2 easy-rbh of the FULL anchor-species proteome (human; Xenopus for
TUBZ1) against it, so a seed is the reciprocal best hit of the anchor among all anchor-species proteins and
paralogs outside the gene list cannot seed a profile. Seeds passing thresholds go to data/seeds/<gene>.faa
(anchor sequence always included).

usage: seeds_collect.py --anchors data/genes/anchors.faa --anchors-tsv data/genes/anchors.tsv --params config/params.yaml --outdir data/seeds \
                        --exclude config/seed_exclude.tsv --threads 8 \
                        --proteomes Hsap=data/proteomes/Hsap.faa Mmus=data/proteomes/Mmus.faa ...
--exclude lists curated (gene, species, accession, reason) RBH seeds to drop: spurious coiled-coil/repeat matches
that would corrupt a profile. Seed files are rewritten only when their content changes, so only the affected
profiles are rebuilt.
"""
import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import anchor_table, die, load_yaml, log, rbh_full_proteome  # noqa: E402
from Bio import SeqIO  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--anchors", required=True)
    ap.add_argument("--anchors-tsv", required=True)
    ap.add_argument("--params", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--proteomes", nargs="+", required=True, help="KEY=path.faa ...")
    ap.add_argument("--exclude", help="TSV with columns gene, species, accession, reason")
    args = ap.parse_args()
    excluded: set[tuple[str, str, str]] = set()
    if args.exclude:
        with open(args.exclude, newline="") as fh:
            excluded = {(r["gene"], r["species"], r["accession"]) for r in csv.DictReader(fh, delimiter="\t")}
    n_excluded = 0
    p = load_yaml(args.params)["rbh"]
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    anchors = {rec.id: rec for rec in SeqIO.parse(args.anchors, "fasta")}
    table = anchor_table(args.anchors, args.anchors_tsv)
    prot_by_key = dict(item.split("=", 1) for item in args.proteomes)
    seeds: dict[str, list] = defaultdict(list)
    for aid, rec in anchors.items():
        gene = aid.split("|")[0]
        rec.id, rec.description = f"{gene}|anchor|{aid.split('|')[1]}", ""
        seeds[gene].append(rec)

    for key, path in prot_by_key.items():
        prot = Path(path)
        index = SeqIO.index(str(prot), "fasta")
        n = 0
        for sp in sorted({a["species"] for a in table.values()}):
            if sp not in prot_by_key:
                die(f"anchor species {sp} is not among the seed proteomes; add it to seed_species")
            wanted = {a["accession"]: g for g, a in table.items() if a["species"] == sp}
            hits = rbh_full_proteome(prot_by_key[sp], prot, wanted, args.threads, p["sensitivity"])
            for h in hits:
                if float(h["evalue"]) > p["evalue_max"]:
                    continue
                if (int(h["qend"]) - int(h["qstart"]) + 1) / int(h["qlen"]) < p["min_query_coverage"]:
                    continue
                gene, t = h["gene"], h["target"]
                rec = index[t]
                if any(str(rec.seq) == str(r.seq) for r in seeds[gene]):
                    continue  # identical to a seed already present (e.g. the anchor itself in the anchor species)
                acc = t.split("|")[1] if "|" in t else t
                if (gene, key, acc) in excluded:
                    n_excluded += 1
                    log(f"  excluded curated seed {gene}|{key}|{acc}")
                    continue
                rec.id, rec.description = f"{gene}|{key}|{acc}", ""
                seeds[gene].append(rec)
                n += 1
        log(f"[{key}] {n} RBH seeds added")

    n_written = 0
    for gene, recs in seeds.items():
        path = outdir / f"{gene}.faa"
        from io import StringIO
        buf = StringIO()
        SeqIO.write(recs, buf, "fasta")
        if not path.exists() or path.read_text() != buf.getvalue():
            path.write_text(buf.getvalue())
            n_written += 1
    log(f"seeds for {len(seeds)} genes in {outdir} ({n_written} files written/updated, {n_excluded} curated exclusions)")


if __name__ == "__main__":
    main()
