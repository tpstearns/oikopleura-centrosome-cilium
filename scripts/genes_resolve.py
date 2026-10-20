#!/usr/bin/env python3
"""Resolve gene symbols in config/gene_modules.tsv to UniProt accessions and anchor sequences.

For each gene: use `anchor_accession` if given; otherwise search UniProt (reviewed entries) in the anchor
species by exact gene symbol, then by each alias. Unresolved genes are written to the report and the run
fails, so nothing downstream proceeds on a guessed anchor.

usage: genes_resolve.py config/gene_modules.tsv data/manifest/species.tsv data/genes/anchors.tsv data/genes/anchors.faa
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import die, load_genes, log, read_tsv, uniprot_fetch_fasta, uniprot_search_tsv, write_tsv  # noqa: E402

FIELDS = ["module", "gene", "anchor", "taxid", "accession", "entry", "gene_primary", "protein_name", "length",
          "resolved_by"]


def search_symbol(symbol: str, taxid: str) -> list[dict]:
    q = f"(gene_exact:{symbol}) AND (organism_id:{taxid}) AND (reviewed:true)"
    return uniprot_search_tsv(q, "accession,id,gene_primary,gene_synonym,protein_name,length", size=10)


def main(genes_tsv: str, species_tsv: str, out_tsv: str, out_faa: str) -> None:
    taxid_by_key = {r["key"]: r["taxid"] for r in read_tsv(species_tsv) if r["taxid"]}
    rows, unresolved = [], []
    for g in load_genes(genes_tsv):
        anchor = g.get("anchor") or "Hsap"
        taxid = taxid_by_key.get(anchor)
        if not taxid:
            die(f"anchor species '{anchor}' for gene {g['gene']} is not resolved in {species_tsv}")
        base = {"module": g["module"], "gene": g["gene"], "anchor": anchor, "taxid": taxid}
        if g.get("anchor_accession"):
            hits = uniprot_search_tsv(f"(accession:{g['anchor_accession']})",
                                      "accession,id,gene_primary,protein_name,length", size=1)
            if not hits:
                unresolved.append((g["gene"], f"accession {g['anchor_accession']} not found"))
                continue
            h = hits[0]
            rows.append({**base, "accession": h["Entry"], "entry": h["Entry Name"], "gene_primary": h.get("Gene Names (primary)", ""),
                         "protein_name": h.get("Protein names", ""), "length": h.get("Length", ""), "resolved_by": "accession"})
            continue
        candidates = [g["gene"]] + [a for a in (g.get("aliases") or "").split("|") if a]
        found = None
        for sym in candidates:
            hits = search_symbol(sym, taxid)
            # prefer a hit whose primary gene name equals the symbol we asked for
            exact = [h for h in hits if h.get("Gene Names (primary)", "").upper() == sym.upper()]
            pick = exact or hits
            if pick:
                found = (pick[0], sym)
                break
        if found is None:
            unresolved.append((g["gene"], f"no reviewed {anchor} entry for {candidates}"))
            continue
        h, sym = found
        rows.append({**base, "accession": h["Entry"], "entry": h["Entry Name"], "gene_primary": h.get("Gene Names (primary)", ""),
                     "protein_name": h.get("Protein names", ""), "length": h.get("Length", ""), "resolved_by": f"symbol:{sym}"})
        log(f"  {g['gene']:<10} -> {h['Entry']} ({h.get('Gene Names (primary)', '')}, {h.get('Length', '')} aa)")

    write_tsv(out_tsv, rows, FIELDS)
    fasta = uniprot_fetch_fasta([r["accession"] for r in rows])
    # Rename FASTA headers to "<gene>|<accession>" so every downstream file carries the gene symbol.
    acc2gene = {r["accession"]: r["gene"] for r in rows}
    out_lines = []
    for line in fasta.splitlines():
        if line.startswith(">"):
            acc = line.split("|")[1] if "|" in line else line[1:].split()[0]
            out_lines.append(f">{acc2gene.get(acc, 'UNKNOWN')}|{acc}")
        else:
            out_lines.append(line)
    Path(out_faa).parent.mkdir(parents=True, exist_ok=True)
    Path(out_faa).write_text("\n".join(out_lines) + "\n")
    log(f"resolved {len(rows)} genes; {len(unresolved)} unresolved")
    if unresolved:
        report = Path(out_tsv).with_name("unresolved.tsv")
        write_tsv(report, [{"gene": g, "reason": r} for g, r in unresolved], ["gene", "reason"])
        die(f"{len(unresolved)} genes unresolved; see {report}. Add anchor_accession in config/gene_modules.tsv.")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        die(__doc__)
    main(*sys.argv[1:])
