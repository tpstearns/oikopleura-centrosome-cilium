#!/usr/bin/env python3
"""Combine evidence for one species into tiers E1–E4 per gene (definitions in docs/evidence_tiers.md).

usage: tiers.py --species KEY --genes config/gene_modules.tsv --params config/params.yaml \
                --hmmsearch results/hmmsearch/KEY.tsv --rbh results/rbh/KEY.tsv --out results/tiers/KEY.tsv \
                [--pfam results/pfam/KEY.tsv] [--miniprot results/miniprot/KEY.tsv] [--transcripts results/tx/KEY.tsv] \
                [--foldseek results/foldseek/KEY.tsv]
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import load_genes, load_yaml, log, read_tsv, write_tsv  # noqa: E402

FIELDS = ["gene", "module", "species", "tier", "best_target", "hmm_evalue", "hmm_profile_coverage", "rbh_target",
          "rbh_evalue", "pfam_shared_fraction", "genome_hit", "genome_identity", "genome_coverage", "transcript_hit",
          "foldseek_hit", "notes"]


def fnum(x: str, default: float) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", required=True)
    ap.add_argument("--genes", required=True)
    ap.add_argument("--params", required=True)
    ap.add_argument("--hmmsearch", required=True)
    ap.add_argument("--rbh", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--pfam")
    ap.add_argument("--miniprot")
    ap.add_argument("--transcripts")
    ap.add_argument("--foldseek")
    args = ap.parse_args()
    P = load_yaml(args.params)
    ph, pr, pp, pm, pt, pf = (P["hmmsearch"], P["rbh"], P["pfam"], P["miniprot"], P["transcripts"], P["foldseek"])
    require_pfam = bool(P.get("tiers", {}).get("require_pfam_for_e1", False))

    # A target protein is credited only to the profile that scores it best ("owner"); hits from other profiles
    # are cross-hits within a family (tubulins, cyclins, E2Fs) and do not count as presence of that gene.
    hits = read_tsv(args.hmmsearch)
    owner: dict[str, str] = {}
    owner_bits: dict[str, float] = {}
    for r in hits:
        if r["target"] not in owner_bits or fnum(r["bits"], 0) > owner_bits[r["target"]]:
            owner[r["target"]], owner_bits[r["target"]] = r["gene"], fnum(r["bits"], 0)
    hmm_best: dict[str, dict] = {}
    cross_best: dict[str, dict] = {}
    for r in hits:
        ok = (fnum(r["evalue"], 1e9) <= ph["evalue_max"] and fnum(r["best_domain_ievalue"], 1e9) <= ph["domain_evalue_max"]
              and fnum(r["profile_coverage"], 0) >= ph["min_profile_coverage"])
        if not ok:
            continue
        if owner[r["target"]] != r["gene"]:
            if r["gene"] not in cross_best or fnum(r["bits"], 0) > fnum(cross_best[r["gene"]]["bits"], 0):
                cross_best[r["gene"]] = r
            continue
        if r["gene"] not in hmm_best or fnum(r["bits"], 0) > fnum(hmm_best[r["gene"]]["bits"], 0):
            hmm_best[r["gene"]] = r

    rbh_best: dict[str, dict] = {}
    for r in read_tsv(args.rbh):
        if fnum(r["evalue"], 1e9) <= pr["evalue_max"] and fnum(r["query_coverage"], 0) >= pr["min_query_coverage"]:
            rbh_best[r["gene"]] = r

    pfam: dict[tuple[str, str], float] = {}
    if args.pfam:
        for r in read_tsv(args.pfam):
            pfam[(r["gene"], r["target"])] = fnum(r["shared_fraction"], 0)

    genome_best: dict[str, dict] = {}
    if args.miniprot:
        for r in read_tsv(args.miniprot):
            if fnum(r["identity"], 0) >= pm["min_identity"] and fnum(r["query_coverage"], 0) >= pm["min_target_coverage"]:
                g = r["gene"]
                if g not in genome_best or fnum(r["identity"], 0) * fnum(r["query_coverage"], 0) > \
                        fnum(genome_best[g]["identity"], 0) * fnum(genome_best[g]["query_coverage"], 0):
                    genome_best[g] = r

    tx_best: dict[str, dict] = {}
    if args.transcripts:
        for r in read_tsv(args.transcripts):
            if fnum(r["evalue"], 1e9) <= pt["evalue_max"] and fnum(r["query_coverage"], 0) >= pt["min_query_coverage"]:
                tx_best[r["gene"]] = r

    fs_best: dict[str, dict] = {}
    if args.foldseek:
        for r in read_tsv(args.foldseek):
            if fnum(r["evalue"], 1e9) <= pf["evalue_max"]:
                fs_best[r["gene"]] = r

    rows = []
    for g in load_genes(args.genes):
        gene = g["gene"]
        h, rb = hmm_best.get(gene), rbh_best.get(gene)
        target = (rb or h or {}).get("target", "")
        shared = pfam.get((gene, target))
        pfam_ok = (shared is None and not require_pfam) or (shared is not None and shared >= pp["min_shared_domain_fraction"])
        notes = []
        if h and rb and h["target"] != rb["target"]:
            notes.append(f"hmm best target {h['target']} differs from RBH target {rb['target']}")
        if not h and gene in cross_best:
            x = cross_best[gene]
            notes.append(f"profile cross-hit only: {x['target']} is best explained by {owner[x['target']]}")
        if rb and owner.get(rb["target"], gene) != gene:
            notes.append(f"RBH target {rb['target']} is best explained by profile {owner[rb['target']]}")
        if shared is None and require_pfam:
            notes.append("no Pfam comparison available")
        if h and rb and pfam_ok:
            tier = "E1"
        elif h or rb or gene in fs_best:
            tier = "E2"
            if h and rb and not pfam_ok:
                notes.append("RBH+HMM but Pfam architecture incomplete")
        elif gene in genome_best or gene in tx_best:
            tier = "E3"
        else:
            tier = "E4"
        gb, tb, fb = genome_best.get(gene), tx_best.get(gene), fs_best.get(gene)
        rows.append({"gene": gene, "module": g["module"], "species": args.species, "tier": tier, "best_target": target,
                     "hmm_evalue": h["evalue"] if h else "", "hmm_profile_coverage": h["profile_coverage"] if h else "",
                     "rbh_target": rb["target"] if rb else "", "rbh_evalue": rb["evalue"] if rb else "",
                     "pfam_shared_fraction": "" if shared is None else f"{shared:.2f}",
                     "genome_hit": f"{gb['contig']}:{gb['start']}-{gb['end']}({gb['strand']})" if gb else "",
                     "genome_identity": gb["identity"] if gb else "", "genome_coverage": gb["query_coverage"] if gb else "",
                     "transcript_hit": tb["target"] if tb else "", "foldseek_hit": fb["target"] if fb else "",
                     "notes": "; ".join(notes)})
    write_tsv(args.out, rows, FIELDS)
    counts = {t: sum(1 for r in rows if r["tier"] == t) for t in ("E1", "E2", "E3", "E4")}
    log(f"[{args.species}] {counts}")


if __name__ == "__main__":
    main()
