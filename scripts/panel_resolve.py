#!/usr/bin/env python3
"""Validate the species panel against NCBI Taxonomy and resolve UniProt proteomes.

Writes data/manifest/species.tsv. Fails on any name/taxid mismatch. Entries whose name is TO-RESOLVE
are written with status 'unresolved' and skipped from downstream runs.

usage: panel_resolve.py config/species_panel.yaml data/manifest/species.tsv
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import die, load_species, log, ncbi_taxid_for_name, uniprot_reference_proteome, write_tsv  # noqa: E402

FIELDS = ["key", "name", "ncbi_name", "taxid", "role", "lineage", "proteome_source", "upid", "proteome_type",
          "protein_count", "genome_source", "status"]


def main(panel_path: str, out_path: str) -> None:
    rows = []
    for sp in load_species(panel_path):
        key, name = sp["key"], sp["name"]
        row = {"key": key, "name": name, "role": sp.get("role", ""), "lineage": sp.get("lineage", ""),
               "proteome_source": (sp.get("proteome") or {}).get("source", ""),
               "genome_source": (sp.get("genome") or {}).get("source", ""), "status": sp.get("status", "")}
        if name == "TO-RESOLVE":
            row.update({"ncbi_name": "", "taxid": "", "upid": "", "proteome_type": "", "protein_count": "",
                        "status": "unresolved"})
            rows.append(row)
            log(f"[{key}] placeholder entry; skipped")
            continue
        taxid, sci = ncbi_taxid_for_name(name)
        if "taxid" in sp and int(sp["taxid"]) != taxid:
            die(f"[{key}] taxid mismatch: panel says {sp['taxid']}, NCBI says {taxid} for '{name}'")
        if sci != name:
            die(f"[{key}] NCBI scientific name '{sci}' differs from panel name '{name}'")
        row.update({"ncbi_name": sci, "taxid": taxid, "upid": "", "proteome_type": "", "protein_count": ""})
        if row["proteome_source"] == "uniprot_reference":
            prot = uniprot_reference_proteome(taxid)
            if prot is None:
                die(f"[{key}] no UniProt proteome for taxid {taxid}; change proteome source in the panel")
            row.update({"upid": prot["id"], "proteome_type": prot.get("proteomeType", ""),
                        "protein_count": prot.get("proteinCount", "")})
        log(f"[{key}] {sci} taxid={taxid} proteome={row['upid'] or row['proteome_source']}")
        rows.append(row)
    write_tsv(out_path, rows, FIELDS)
    log(f"wrote {out_path} ({len(rows)} species)")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        die(__doc__)
    main(sys.argv[1], sys.argv[2])
