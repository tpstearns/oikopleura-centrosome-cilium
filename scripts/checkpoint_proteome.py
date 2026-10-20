#!/usr/bin/env python3
"""Proteome-level checkpoint: per gene x species, what the profile search and RBH found, before genome searches.

This is NOT the evidence-tier matrix. It classifies proteome evidence only:
  HMM+RBH      profile hit passing thresholds and reciprocal best hit to the anchor
  HMM          profile hit only      RBH   reciprocal best hit only
  HMM-partial  strong profile hit (E <= partial_evalue) whose profile coverage is below the threshold
  none         no proteome hit at thresholds: genome search pending, never reported as absent
A profile hit counts only if this gene's profile is the best-scoring profile for that target protein ("owner");
otherwise it is a cross-hit from a paralog's profile and is reported as such (tubulins, cyclins, E2Fs...).
A near-miss column shows the best raw profile hit (any E <= 10) so divergent candidates are visible.

usage: checkpoint_proteome.py --genes config/gene_modules.tsv --params config/params.yaml \
         --species Odio_OKI Crob ... [--hmmsearch-dir results/hmmsearch] [--rbh-dir results/rbh] \
         [--miniprot-dir results/miniprot] --detail A_tubulins I_multiciliogenesis --out results/x.tsv --md results/x.md
"""
import argparse
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import die, load_genes, load_yaml, log, read_tsv, write_tsv  # noqa: E402

FIELDS = ["species", "module", "gene", "proteome_class", "hmm_target", "hmm_evalue", "hmm_profile_cov", "hmm_target_cov",
          "rbh_target", "rbh_evalue", "rbh_query_cov", "rbh_pident", "best_raw_target", "best_raw_evalue", "best_raw_profile_cov",
          "genome_hit", "genome_identity", "genome_query_cov", "note"]


def f(x, default=float("inf")):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes", required=True)
    ap.add_argument("--params", required=True)
    ap.add_argument("--species", nargs="+", required=True)
    ap.add_argument("--hmmsearch-dir", default="results/hmmsearch")
    ap.add_argument("--rbh-dir", default="results/rbh")
    ap.add_argument("--miniprot-dir", default="results/miniprot")
    ap.add_argument("--detail", nargs="*", default=[])
    ap.add_argument("--partial-evalue", type=float, default=1e-10,
                    help="E-value for the HMM-partial class (strong hit, profile coverage below threshold)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--md", required=True)
    args = ap.parse_args()
    P = load_yaml(args.params)
    ph, pr, pm = P["hmmsearch"], P["rbh"], P["miniprot"]
    genes = load_genes(args.genes)

    rows = []
    for sp in args.species:
        hpath, rpath = Path(args.hmmsearch_dir) / f"{sp}.tsv", Path(args.rbh_dir) / f"{sp}.tsv"
        if not hpath.exists() or not rpath.exists():
            die(f"[{sp}] missing {hpath} or {rpath}")
        hits = read_tsv(hpath)
        owner: dict[str, str] = {}   # target -> gene whose profile scores it best
        owner_bits: dict[str, float] = {}
        for r in hits:
            if r["target"] not in owner_bits or f(r["bits"], 0) > owner_bits[r["target"]]:
                owner[r["target"]], owner_bits[r["target"]] = r["gene"], f(r["bits"], 0)
        best_pass: dict[str, dict] = {}
        best_partial: dict[str, dict] = {}
        best_raw: dict[str, dict] = {}
        cross: dict[str, dict] = {}   # gene -> best passing hit that belongs to another profile
        for r in hits:
            g = r["gene"]
            if g not in best_raw or f(r["evalue"]) < f(best_raw[g]["evalue"]):
                best_raw[g] = r
            strong = f(r["evalue"]) <= ph["evalue_max"] and f(r["best_domain_ievalue"]) <= ph["domain_evalue_max"]
            covered = f(r["profile_coverage"], 0) >= ph["min_profile_coverage"]
            if strong and covered and owner[r["target"]] != g:
                if g not in cross or f(r["bits"], 0) > f(cross[g]["bits"], 0):
                    cross[g] = r
                continue
            if strong and covered and (g not in best_pass or f(r["bits"], 0) > f(best_pass[g]["bits"], 0)):
                best_pass[g] = r
            elif (f(r["evalue"]) <= args.partial_evalue and owner[r["target"]] == g
                  and (g not in best_partial or f(r["bits"], 0) > f(best_partial[g]["bits"], 0))):
                best_partial[g] = r
        rbh: dict[str, dict] = {}
        for r in read_tsv(rpath):
            if f(r["evalue"]) <= pr["evalue_max"] and f(r["query_coverage"], 0) >= pr["min_query_coverage"]:
                rbh[r["gene"]] = r
        genome: dict[str, dict] = {}
        mpath = Path(args.miniprot_dir) / f"{sp}.tsv"
        if mpath.exists():
            for r in read_tsv(mpath):
                if f(r["identity"], 0) >= pm["min_identity"] and f(r["query_coverage"], 0) >= pm["min_target_coverage"]:
                    g = r["gene"]
                    score = f(r["identity"], 0) * f(r["query_coverage"], 0)
                    if g not in genome or score > f(genome[g]["identity"], 0) * f(genome[g]["query_coverage"], 0):
                        genome[g] = r
        for g in genes:
            gene = g["gene"]
            h, rb, raw, gm = best_pass.get(gene), rbh.get(gene), best_raw.get(gene), genome.get(gene)
            part, x = best_partial.get(gene), cross.get(gene)
            cls = "HMM+RBH" if h and rb else "HMM" if h else "RBH" if rb else "HMM-partial" if part else "none"
            notes = []
            if cls == "HMM-partial":
                notes.append(f"partial: E={part['evalue']} profile cov={part['profile_coverage']}")
            if x and not h:
                notes.append(f"cross-hit only: {x['target']} is best explained by {owner[x['target']]} (E={x['evalue']})")
            if h and rb and h["target"] != rb["target"]:
                notes.append(f"HMM best {h['target']} != RBH {rb['target']}")
            if rb and owner.get(rb["target"], gene) != gene:
                notes.append(f"RBH target {rb['target']} is best explained by profile {owner[rb['target']]}")
            if cls == "none" and raw and f(raw["evalue"]) <= 1e-3:
                notes.append(f"near-miss E={raw['evalue']} cov={raw['profile_coverage']}")
            if cls == "none" and gm:
                notes.append("genome hit (miniprot)")
            rows.append({"species": sp, "module": g["module"], "gene": gene, "proteome_class": cls,
                         "hmm_target": h["target"] if h else "", "hmm_evalue": h["evalue"] if h else "",
                         "hmm_profile_cov": h["profile_coverage"] if h else "", "hmm_target_cov": h["target_coverage"] if h else "",
                         "rbh_target": rb["target"] if rb else "", "rbh_evalue": rb["evalue"] if rb else "",
                         "rbh_query_cov": rb["query_coverage"] if rb else "", "rbh_pident": rb["pident"] if rb else "",
                         "best_raw_target": raw["target"] if raw else "", "best_raw_evalue": raw["evalue"] if raw else "",
                         "best_raw_profile_cov": raw["profile_coverage"] if raw else "",
                         "genome_hit": f"{gm['contig']}:{gm['start']}-{gm['end']}" if gm else "",
                         "genome_identity": gm["identity"] if gm else "", "genome_query_cov": gm["query_coverage"] if gm else "",
                         "note": "; ".join(notes)})
    write_tsv(args.out, rows, FIELDS)

    # ---- markdown
    sps = args.species
    by_mod: dict[str, list] = defaultdict(list)
    for g in genes:
        by_mod[g["module"]].append(g["gene"])
    cell = {(r["species"], r["gene"]): r for r in rows}
    have_genome = {sp for sp in sps if (Path(args.miniprot_dir) / f"{sp}.tsv").exists()}
    md = ["# Proteome-search checkpoint", "",
          "Proteome evidence only (hmmsearch with the per-gene profiles + MMseqs2 RBH to the human anchors, thresholds from "
          "`config/params.yaml`). `none` means no proteome hit at thresholds and is **not** an absence call: the genome "
          "search decides. Table: `" + args.out + "`.", "",
          "Classes: **HMM+RBH** both; **HMM** profile only; **RBH** reciprocal best hit only (reciprocal against the full "
          "anchor-species proteome); **HMM-partial** profile hit with E ≤ " + f"{args.partial_evalue:g}" + " but profile coverage below "
          f"{ph['min_profile_coverage']}; **none** no proteome hit"
          + (" (genome hit noted where miniprot found one)" if have_genome else "") + ".", "",
          "## Genes with proteome evidence: HMM+RBH / HMM or RBH only / HMM-partial, out of module size", "",
          "| module | n | " + " | ".join(sps) + " |", "|---|---|" + "---|" * len(sps)]
    for m, gl in by_mod.items():
        vals = []
        for sp in sps:
            c = Counter(cell[(sp, g)]["proteome_class"] for g in gl)
            vals.append(f"{c['HMM+RBH']} / {c['HMM'] + c['RBH']} / {c['HMM-partial']}")
        md.append(f"| {m} | {len(gl)} | " + " | ".join(vals) + " |")
    md += ["", "Totals: genes with HMM/RBH evidence at thresholds (+ HMM-partial) / 272: " + ", ".join(
        f"{sp} {sum(1 for g in genes if cell[(sp, g['gene'])]['proteome_class'] in ('HMM+RBH', 'HMM', 'RBH'))}"
        f" (+{sum(1 for g in genes if cell[(sp, g['gene'])]['proteome_class'] == 'HMM-partial')})" for sp in sps) + "."]

    sym = {"HMM+RBH": "●", "HMM": "◐", "RBH": "◑", "HMM-partial": "◔", "none": "○"}
    for m in args.detail:
        if m not in by_mod:
            die(f"module {m} not in gene list")
        md += ["", f"## {m}: per-gene detail", "", "● HMM+RBH  ◐ HMM only  ◑ RBH only  ◔ HMM-partial (strong E, low profile coverage)  ○ no proteome hit (genome search pending)"
               + ("; G = miniprot genome hit at thresholds" if have_genome else ""), "",
               "| gene | " + " | ".join(sps) + " |", "|---|" + "---|" * len(sps)]
        for g in by_mod[m]:
            vals = []
            for sp in sps:
                r = cell[(sp, g)]
                s = sym[r["proteome_class"]]
                if r["proteome_class"] in ("HMM+RBH", "HMM", "RBH"):
                    ev = r["hmm_evalue"] or r["rbh_evalue"]
                    s += f" {ev}"
                elif r["proteome_class"] == "HMM-partial":
                    s += " " + r["note"].split("partial: ")[1].split(";")[0]
                elif r["genome_hit"]:
                    s += " G"
                elif r["best_raw_evalue"] and f(r["best_raw_evalue"]) <= 1e-3:
                    s += f" ({r['best_raw_evalue']}, cov {r['best_raw_profile_cov']})"
                vals.append(s)
            md.append(f"| {g} | " + " | ".join(vals) + " |")
        notes = [(g, sp, cell[(sp, g)]["note"]) for g in by_mod[m] for sp in sps if cell[(sp, g)]["note"]]
        if notes:
            md += ["", "Notes:"] + [f"- {g} / {sp}: {n}" for g, sp, n in notes]

    nohit_all = [g["gene"] for g in genes if all(cell[(sp, g["gene"])]["proteome_class"] in ("none", "HMM-partial") for sp in sps)]
    md += ["", f"## Genes with no proteome hit at thresholds in any of these {len(sps)} species ({len(nohit_all)}; HMM-partial allowed)", "",
           "Candidates for the genome search; listed by module.", ""]
    for m, gl in by_mod.items():
        miss = [g for g in gl if g in nohit_all]
        if miss:
            md.append(f"- {m}: " + ", ".join(miss))
    Path(args.md).write_text("\n".join(md) + "\n")
    log(f"checkpoint: {len(rows)} cells -> {args.out}, {args.md}")


if __name__ == "__main__":
    main()
