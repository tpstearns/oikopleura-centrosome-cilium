#!/usr/bin/env python3
"""Pfam domain architecture of candidate hits versus their anchors.

Candidates are the targets listed in the hmmsearch and RBH tables for one species. Each candidate and each
anchor is scanned against Pfam-A (hmmsearch with the Pfam profiles over the small candidate set), and the
fraction of the anchor's Pfam families found in the candidate is reported.

usage: pfam_scan.py --pfam data/pfam/Pfam-A.hmm --anchors data/genes/anchors.faa --proteome data/proteomes/<key>.faa \
                    --hmmsearch results/hmmsearch/<key>.tsv --rbh results/rbh/<key>.tsv --out results/pfam/<key>.tsv \
                    --params config/params.yaml --threads 8
Output columns: gene, target, anchor_domains, target_domains, shared_fraction
"""
import argparse
import sys
from collections import defaultdict
from pathlib import Path

import pyhmmer

sys.path.insert(0, str(Path(__file__).parent))
from common import load_yaml, log, read_tsv, text, write_tsv  # noqa: E402

FIELDS = ["gene", "target", "anchor_domains", "target_domains", "shared_fraction"]


def domains_for(seqs, pfam_path: str, threads: int, evalue: float) -> dict[str, set[str]]:
    """Map sequence name -> set of Pfam family accessions (PFxxxxx) with at least one domain under evalue."""
    doms: dict[str, set[str]] = defaultdict(set)
    with pyhmmer.plan7.HMMFile(pfam_path) as hf:
        for hits in pyhmmer.hmmsearch(hf, seqs, cpus=threads, E=evalue, domE=evalue):
            fam = text(hits.query.accession or hits.query.name).split(".")[0]
            for hit in hits:
                if any(d.i_evalue <= evalue for d in hit.domains if d.i_evalue is not None):
                    doms[text(hit.name)].add(fam)
    return doms


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pfam", required=True)
    ap.add_argument("--anchors", required=True)
    ap.add_argument("--proteome", required=True)
    ap.add_argument("--hmmsearch", required=True)
    ap.add_argument("--rbh", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--params", required=True)
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()
    p = load_yaml(args.params)
    ev = float(p["pfam"]["evalue_max"])
    loose = float(p["hmmsearch"]["evalue_max"]) * 1e3  # scan a generous candidate set

    pairs: set[tuple[str, str]] = set()
    for r in read_tsv(args.hmmsearch):
        if float(r["evalue"]) <= loose:
            pairs.add((r["gene"], r["target"]))
    for r in read_tsv(args.rbh):
        pairs.add((r["gene"], r["target"]))
    targets = {t for _, t in pairs}

    alphabet = pyhmmer.easel.Alphabet.amino()
    with pyhmmer.easel.SequenceFile(args.anchors, digital=True, alphabet=alphabet) as sf:
        anchor_seqs = sf.read_block()
    with pyhmmer.easel.SequenceFile(args.proteome, digital=True, alphabet=alphabet) as sf:
        cand = pyhmmer.easel.DigitalSequenceBlock(alphabet, [s for s in sf if text(s.name) in targets])
    log(f"scanning {len(anchor_seqs)} anchors and {len(cand)} candidates against Pfam-A")
    adoms = domains_for(anchor_seqs, args.pfam, args.threads, ev)
    cdoms = domains_for(cand, args.pfam, args.threads, ev)
    anchor_by_gene = {text(s.name).split("|")[0]: text(s.name) for s in anchor_seqs}

    rows = []
    for gene, target in sorted(pairs):
        a = adoms.get(anchor_by_gene.get(gene, ""), set())
        c = cdoms.get(target, set())
        shared = (len(a & c) / len(a)) if a else (1.0 if not c else 0.0)
        rows.append({"gene": gene, "target": target, "anchor_domains": ";".join(sorted(a)),
                     "target_domains": ";".join(sorted(c)), "shared_fraction": f"{shared:.2f}"})
    write_tsv(args.out, rows, FIELDS)
    log(f"wrote {len(rows)} gene-target architecture comparisons to {args.out}")


if __name__ == "__main__":
    main()
