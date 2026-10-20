#!/usr/bin/env python3
"""Assemble the gene × species evidence matrix and module summaries from per-species tier tables.

usage: matrix.py --genes config/gene_modules.tsv --tiers results/tiers/*.tsv --out results/matrix.tsv \
                 --summary results/module_summary.tsv --md results/matrix_summary.md
"""
import argparse
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import load_genes, log, read_tsv, write_tsv  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", required=True)
    ap.add_argument("--tiers", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--md", required=True)
    args = ap.parse_args()

    genes = load_genes(args.genes)
    species, cells = [], {}
    for path in args.tiers:
        rows = read_tsv(path)
        if not rows:
            continue
        sp = rows[0]["species"]
        species.append(sp)
        for r in rows:
            cells[(r["gene"], sp)] = r["tier"]
    species = sorted(set(species))

    matrix_rows = [{"module": g["module"], "gene": g["gene"], **{sp: cells.get((g["gene"], sp), "") for sp in species}}
                   for g in genes]
    write_tsv(args.out, matrix_rows, ["module", "gene"] + species)

    summary = []
    by_module: dict[str, list] = defaultdict(list)
    for g in genes:
        by_module[g["module"]].append(g["gene"])
    for module, gl in by_module.items():
        for sp in species:
            c = Counter(cells.get((gene, sp), "NA") for gene in gl)
            summary.append({"module": module, "species": sp, "n_genes": len(gl), "E1": c["E1"], "E2": c["E2"],
                            "E3": c["E3"], "E4": c["E4"], "NA": c["NA"],
                            "present_fraction": f"{(c['E1'] + c['E2'] + c['E3']) / len(gl):.2f}"})
    write_tsv(args.summary, summary, ["module", "species", "n_genes", "E1", "E2", "E3", "E4", "NA", "present_fraction"])

    lines = ["# Module retention summary", "", "Present = E1+E2+E3; absent = E4 (see docs/evidence_tiers.md).", "",
             "| module | " + " | ".join(species) + " |", "|---|" + "---|" * len(species)]
    for module, gl in by_module.items():
        vals = []
        for sp in species:
            c = Counter(cells.get((gene, sp), "NA") for gene in gl)
            vals.append(f"{c['E1'] + c['E2'] + c['E3']}/{len(gl)} (E4 {c['E4']})")
        lines.append(f"| {module} | " + " | ".join(vals) + " |")
    Path(args.md).write_text("\n".join(lines) + "\n")
    log(f"matrix {len(genes)} genes × {len(species)} species -> {args.out}")


if __name__ == "__main__":
    main()
