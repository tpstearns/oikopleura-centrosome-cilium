#!/usr/bin/env python3
"""Tubulin-family phylogeny to place the Oikopleura tubulins against reference alpha/beta/gamma/delta/epsilon/zeta clades.

Reference set: the per-gene seed files for TUBA1A, TUBB, TUBG1, TUBD1, TUBE1, TUBZ1 (anchor + RBH orthologs), extra
reference FASTAs given with --extra (label taken from the file name), and, for the species in --owner-species, every
protein that the TUBD1/TUBE1/TUBZ1 profile scores best among the tubulin profiles (labelled "<species>_<family>?").
Query set: every protein of the --query-species proteomes hit by any tubulin profile at E <= --query-evalue.
Pipeline: MAFFT L-INS-i -> column trimming (>= --min-occupancy non-gap) -> IQ-TREE (ModelFinder, UFBoot, SH-aLRT) ->
clade report (for each reference family: MRCA support, monophyly, and which query sequences fall inside).

usage: tubulin_tree.py --profiles data/profiles --seeds data/seeds --proteomes-dir data/proteomes \\
         --query-species Odio_OKI Odio_OdB3 --owner-species Bflo Clep Cint Pmam Bsch Scla Tdem \\
         --extra Xla_zeta=path.faa --outdir results/tubulin_tree --threads 12
"""
import argparse
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import pyhmmer
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, str(Path(__file__).parent))
from common import die, log, text  # noqa: E402

FAMILIES = {"TUBA1A": "alpha", "TUBB": "beta", "TUBG1": "gamma", "TUBD1": "delta", "TUBE1": "epsilon", "TUBZ1": "zeta"}


def clean(label: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", label)


def owner_hits(profiles_dir: Path, proteome: Path, evalue: float, threads: int):
    """target -> (owning family gene, bits, evalue) over the six family profiles (+TUBG2 as a gamma alias)."""
    alphabet = pyhmmer.easel.Alphabet.amino()
    hmms = []
    for g in list(FAMILIES) + ["TUBG2"]:
        with pyhmmer.plan7.HMMFile(profiles_dir / f"{g}.hmm") as hf:
            hmms.append(next(hf))
    with pyhmmer.easel.SequenceFile(str(proteome), digital=True, alphabet=alphabet) as sf:
        seqs = sf.read_block()
    best = {}
    for hits in pyhmmer.hmmsearch(hmms, seqs, cpus=threads, E=evalue):
        g = text(hits.query.name)
        g = "TUBG1" if g == "TUBG2" else g
        for h in hits:
            n = text(h.name)
            if n not in best or h.score > best[n][1]:
                best[n] = (g, h.score, h.evalue)
    return best


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--profiles", default="data/profiles")
    ap.add_argument("--seeds", default="data/seeds")
    ap.add_argument("--proteomes-dir", default="data/proteomes")
    ap.add_argument("--query-species", nargs="+", required=True)
    ap.add_argument("--owner-species", nargs="*", default=[])
    ap.add_argument("--extra", nargs="*", default=[], help="LABEL=fasta with extra reference sequences")
    ap.add_argument("--query-evalue", type=float, default=1e-30)
    ap.add_argument("--min-length", type=int, default=300)
    ap.add_argument("--min-occupancy", type=float, default=0.5)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    profiles, seeds, pdir = Path(args.profiles), Path(args.seeds), Path(args.proteomes_dir)

    records, seen, labels = [], {}, {}   # labels: label -> (kind, family, species)

    def add(label, seq, kind, family, species):
        label = clean(label)
        if len(seq) < args.min_length:
            return
        s = str(seq)
        if s in seen:
            labels[seen[s]] = labels[seen[s]] + (label,)  # record duplicate names
            return
        seen[s] = label
        labels[label] = (kind, family, species)
        records.append(SeqRecord(Seq(s), id=label, description=""))

    # 1. references from seed files
    for gene, fam in FAMILIES.items():
        for rec in SeqIO.parse(seeds / f"{gene}.faa", "fasta"):
            _, sp, acc = rec.id.split("|")[:3]
            sp = "Hsap" if sp == "anchor" and gene != "TUBZ1" else ("Xtro" if sp == "anchor" else sp)
            add(f"REF_{fam}_{sp}_{acc}", rec.seq, "ref", fam, sp)
    # 2. extra references
    for item in args.extra:
        label, path = item.split("=", 1)
        fam = label.split("_")[-1]
        for rec in SeqIO.parse(path, "fasta"):
            acc = rec.id.split("|")[1] if "|" in rec.id else rec.id
            add(f"REF_{fam}_{label.split('_')[0]}_{acc}", rec.seq, "ref", fam, label.split("_")[0])
    # 3. owner-assigned delta/epsilon/zeta proteins from other panel species (candidate references)
    for sp in args.owner_species:
        best = owner_hits(profiles, pdir / f"{sp}.faa", 1e-20, args.threads)
        index = SeqIO.index(str(pdir / f"{sp}.faa"), "fasta")
        for n, (g, bits, ev) in best.items():
            if g in ("TUBD1", "TUBE1", "TUBZ1"):
                acc = n.split("|")[1] if "|" in n else n
                add(f"CAND_{FAMILIES[g]}_{sp}_{acc}", index[n].seq, "cand", FAMILIES[g], sp)
    # 4. queries: all tubulin-family proteins of the query species
    for sp in args.query_species:
        best = owner_hits(profiles, pdir / f"{sp}.faa", args.query_evalue, args.threads)
        index = SeqIO.index(str(pdir / f"{sp}.faa"), "fasta")
        for n, (g, bits, ev) in best.items():
            acc = n.split("|")[1] if "|" in n else n
            add(f"OD_{sp}_{acc}_best-{FAMILIES[g]}", index[n].seq, "query", "", sp)
    fasta = out / "tubulins.faa"
    SeqIO.write(records, fasta, "fasta")
    log(f"{len(records)} unique sequences (refs {sum(1 for v in labels.values() if v[0] == 'ref')}, "
        f"candidates {sum(1 for v in labels.values() if v[0] == 'cand')}, queries {sum(1 for v in labels.values() if v[0] == 'query')})")

    # alignment and trimming
    aln = out / "tubulins.linsi.faa"
    with open(aln, "w") as fh:
        subprocess.run(["mafft-linsi", "--thread", str(args.threads), "--quiet", str(fasta)], check=True, stdout=fh)
    msa = list(SeqIO.parse(aln, "fasta"))
    L = len(msa[0].seq)
    keep = [i for i in range(L) if sum(1 for r in msa if r.seq[i] != "-") / len(msa) >= args.min_occupancy]
    trimmed = out / "tubulins.trim.faa"
    SeqIO.write([SeqRecord(Seq("".join(str(r.seq)[i] for i in keep)), id=r.id, description="") for r in msa], trimmed, "fasta")
    log(f"alignment {L} columns -> {len(keep)} kept at occupancy >= {args.min_occupancy}")

    # tree
    prefix = out / "tubulins"
    cmd = ["iqtree3", "-s", str(trimmed), "-m", "MFP", "-B", "1000", "-alrt", "1000", "-T", str(args.threads),
           "--prefix", str(prefix), "-redo", "--quiet"]
    log("  " + " ".join(cmd))
    subprocess.run(cmd, check=True)
    treefile = Path(f"{prefix}.treefile")

    # clade report
    from ete3 import Tree
    t = Tree(treefile.read_text(), format=1)
    t.set_outgroup(t.get_midpoint_outgroup())
    model = ""
    for line in open(f"{prefix}.iqtree"):
        if line.startswith("Best-fit model"):
            model = line.strip()
    lines = ["# Tubulin-family tree", "", f"{len(records)} sequences, {len(keep)} trimmed alignment columns; IQ-TREE 3, {model}; "
             "support = SH-aLRT % / ultrafast bootstrap %; midpoint-rooted for display.", "",
             "Reference labels: REF_<family>_<species>_<accession> (seed-set orthologs; TUBZ1 anchor is X. tropicalis), "
             "CAND_<family>_<species> (proteins owned by the delta/epsilon/zeta profile in other panel proteomes), "
             "OD_<lineage>_<accession>_best-<family> (every Oikopleura tubulin-family protein; suffix = best-scoring profile).", "",
             "## Reference clades and the Oikopleura sequences inside them", "",
             "| family | reference seqs | MRCA support (aLRT/UFBoot) | monophyletic (allowing CAND/OD inside) | CAND inside | OD inside | foreign REF inside |",
             "|---|---|---|---|---|---|---|"]
    leaf_kind = {l: labels[l][0] for l in t.get_leaf_names() if l in labels}
    leaf_fam = {l: labels[l][1] for l in t.get_leaf_names() if l in labels}
    od_placement = {}
    for fam in FAMILIES.values():
        refs = [l for l in t.get_leaf_names() if leaf_kind.get(l) == "ref" and leaf_fam.get(l) == fam]
        if len(refs) < 2:
            lines.append(f"| {fam} | {len(refs)} | n/a | n/a | | | |")
            continue
        mrca = t.get_common_ancestor(refs)
        inside = mrca.get_leaf_names()
        foreign = [l for l in inside if leaf_kind.get(l) == "ref" and leaf_fam.get(l) != fam]
        cands = [l for l in inside if leaf_kind.get(l) == "cand"]
        ods = [l for l in inside if leaf_kind.get(l) == "query"]
        for l in ods:
            od_placement.setdefault(l, []).append(fam)
        sup = mrca.name if mrca.name else "root"
        lines.append(f"| {fam} | {len(refs)} | {sup} | {'yes' if not foreign else 'no'} | {len(cands)}: {', '.join(c.replace('CAND_', '') for c in cands) or '—'} | "
                     f"{len(ods)}: {', '.join(o.replace('OD_', '') for o in ods) or '—'} | {', '.join(foreign) or '—'} |")
    lines += ["", "## Every Oikopleura sequence: best profile vs tree placement", "", "| sequence | best profile | inside reference clade(s) |", "|---|---|---|"]
    for l in sorted(l for l in t.get_leaf_names() if leaf_kind.get(l) == "query"):
        lines.append(f"| {l.replace('OD_', '')} | {l.split('best-')[-1]} | {', '.join(od_placement.get(l, [])) or 'none (outside all reference clades)'} |")
    dups = {k: v[3:] for k, v in labels.items() if len(v) > 3}
    if dups:
        lines += ["", "Identical sequences collapsed: " + "; ".join(f"{k} = {', '.join(v)}" for k, v in dups.items())]
    lines += ["", f"Files: `{treefile.name}` (Newick), `{trimmed.name}`, `{aln.name}`, `{Path(str(prefix) + '.iqtree').name}` (IQ-TREE report), `tubulin_tree.pdf`."]
    (out / "tubulin_tree.md").write_text("\n".join(lines) + "\n")

    # figure (Bio.Phylo passes the label *string* to label_colors)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from Bio import Phylo
        tree = Phylo.read(str(treefile), "newick")
        tree.root_at_midpoint()
        tree.ladderize()
        n = len(tree.get_terminals())
        fig = plt.figure(figsize=(12, 0.22 * n + 2))
        ax = fig.add_subplot(1, 1, 1)
        colors = {"alpha": "#999999", "beta": "#bbbbbb", "gamma": "#7f7f7f", "delta": "#1f77b4", "epsilon": "#2ca02c", "zeta": "#d62728"}

        def lab(cl):
            if not cl.name or not cl.is_terminal():
                return None
            return cl.name.replace("REF_", "").replace("CAND_", "cand:").replace("OD_", "Od:")

        def col(name):
            if name.startswith("Od:"):
                return "#ff7f0e"
            for fam, c in colors.items():
                if name.startswith((f"{fam}_", f"cand:{fam}_")):
                    return c
            return "black"
        Phylo.draw(tree, axes=ax, do_show=False, label_func=lab, label_colors=col, branch_labels=lambda c: None)
        ax.set_title(f"Tubulin family tree ({n} sequences; IQ-TREE, {model.split(': ')[-1] if ': ' in model else ''}). "
                     "Oikopleura orange; delta blue, epsilon green, zeta red; alpha/beta/gamma grey")
        ax.set_xlabel("substitutions per site")
        fig.tight_layout()
        fig.savefig(out / "tubulin_tree.pdf")
        Path("results/figures").mkdir(parents=True, exist_ok=True)
        fig.savefig("results/figures/tubulin_tree.pdf")
        fig.savefig("results/figures/tubulin_tree.png", dpi=110)
    except Exception as exc:  # noqa: BLE001
        log(f"figure not drawn: {exc}")
    log(f"wrote {out / 'tubulin_tree.md'}")


if __name__ == "__main__":
    main()
