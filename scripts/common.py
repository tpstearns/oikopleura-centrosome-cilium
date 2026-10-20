"""Shared helpers: config loading, manifests, checksums, HTTP with retries, UniProt/NCBI queries.

All functions fail loudly. Nothing here guesses identifiers.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import sys
import time
from pathlib import Path
from typing import Any, Iterable

import requests
import yaml

UNIPROT_REST = "https://rest.uniprot.org"
UNIPROT_FTP_REF = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes"
NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ZENODO_API = "https://zenodo.org/api/records"
USER_AGENT = "oikopleura-centrosome-cilium/0.1 (research pipeline)"


def text(x) -> str:
    """pyhmmer returns bytes for names/descriptions in versions <0.12 and str from 0.12; normalise to str."""
    if x is None:
        return ""
    return x.decode() if isinstance(x, (bytes, bytearray)) else str(x)


def set_text(obj, attr: str, value: str) -> None:
    """Assign a name/description on a pyhmmer object, matching the type the installed version expects."""
    try:
        setattr(obj, attr, value)
    except TypeError:
        setattr(obj, attr, value.encode())


def die(msg: str, code: int = 1) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(code)


def log(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)


def load_yaml(path: str | Path) -> dict:
    with open(path) as fh:
        return yaml.safe_load(fh)


def load_species(panel_path: str | Path) -> list[dict]:
    panel = load_yaml(panel_path)
    species = panel.get("species", [])
    keys = [s["key"] for s in species]
    dups = {k for k in keys if keys.count(k) > 1}
    if dups:
        die(f"duplicate species keys in panel: {sorted(dups)}")
    return species


def species_by_key(panel_path: str | Path, key: str) -> dict:
    for s in load_species(panel_path):
        if s["key"] == key:
            return s
    die(f"species key not in panel: {key}")


def load_genes(tsv_path: str | Path) -> list[dict]:
    with open(tsv_path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    genes = [r["gene"] for r in rows]
    dups = {g for g in genes if genes.count(g) > 1}
    if dups:
        die(f"duplicate gene symbols in gene list: {sorted(dups)}")
    return rows


def sha256_of(path: str | Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def http_get(url: str, params: dict | None = None, stream: bool = False, retries: int = 4,
             timeout: int = 120, headers: dict | None = None) -> requests.Response:
    hdrs = {"User-Agent": USER_AGENT}
    if headers:
        hdrs.update(headers)
    last = None
    for attempt in range(retries):
        try:
            r = requests.get(url, params=params, stream=stream, timeout=timeout, headers=hdrs)
            if r.status_code in (429, 500, 502, 503, 504):
                raise requests.HTTPError(f"{r.status_code} from {url}")
            r.raise_for_status()
            return r
        except Exception as exc:  # noqa: BLE001
            last = exc
            wait = 2 ** attempt
            log(f"  retry {attempt + 1}/{retries} after error: {exc} (sleep {wait}s)")
            time.sleep(wait)
    die(f"GET failed after {retries} attempts: {url} ({last})")


def download(url: str, dest: str | Path, params: dict | None = None) -> Path:
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    r = http_get(url, params=params, stream=True, timeout=600)
    with open(tmp, "wb") as fh:
        for block in r.iter_content(1 << 20):
            fh.write(block)
    tmp.replace(dest)
    return dest


def record_file(manifest: str | Path, path: str | Path, source: str, note: str = "") -> None:
    """Append (or replace) a row in data/manifest/files.tsv with the file's checksum."""
    manifest = Path(manifest)
    manifest.parent.mkdir(parents=True, exist_ok=True)
    path = Path(path)
    row = {
        "path": str(path),
        "source": source,
        "sha256": sha256_of(path),
        "bytes": str(path.stat().st_size),
        "date": time.strftime("%Y-%m-%d"),
        "note": note,
    }
    rows: list[dict] = []
    if manifest.exists():
        with open(manifest, newline="") as fh:
            rows = [r for r in csv.DictReader(fh, delimiter="\t") if r["path"] != row["path"]]
    rows.append(row)
    with open(manifest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(row.keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def write_tsv(path: str | Path, rows: Iterable[dict], fieldnames: list[str]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def read_tsv(path: str | Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ---------------------------------------------------------------- NCBI Taxonomy

NCBI_MIN_INTERVAL = 0.4  # seconds between E-utilities calls (NCBI allows 3 requests/s without an API key)
_last_ncbi_call = 0.0


def _ncbi_throttle() -> None:
    global _last_ncbi_call
    wait = NCBI_MIN_INTERVAL - (time.monotonic() - _last_ncbi_call)
    if wait > 0:
        time.sleep(wait)
    _last_ncbi_call = time.monotonic()


def ncbi_taxid_for_name(name: str) -> tuple[int, str]:
    """Resolve a scientific name to an NCBI taxon ID; returns (taxid, scientific_name). Exact-name match only."""
    _ncbi_throttle()
    r = http_get(f"{NCBI_EUTILS}/esearch.fcgi",
                 params={"db": "taxonomy", "term": f'"{name}"[Scientific Name]', "retmode": "json"})
    ids = r.json().get("esearchresult", {}).get("idlist", [])
    if len(ids) != 1:
        die(f"NCBI Taxonomy: expected exactly one match for '{name}', got {ids}")
    taxid = int(ids[0])
    _ncbi_throttle()
    r2 = http_get(f"{NCBI_EUTILS}/efetch.fcgi", params={"db": "taxonomy", "id": taxid, "retmode": "xml"})
    import xml.etree.ElementTree as ET
    root = ET.fromstring(r2.text)
    sci = root.findtext(".//Taxon/ScientificName")
    if sci is None:
        die(f"NCBI Taxonomy: no ScientificName returned for taxid {taxid}")
    return taxid, sci


# ---------------------------------------------------------------- UniProt

def uniprot_reference_proteome(taxid: int) -> dict | None:
    """Return the UniProt proteome record for a taxon ID, preferring reference proteomes. None if none exists."""
    r = http_get(f"{UNIPROT_REST}/proteomes/search",
                 params={"query": f"(organism_id:{taxid})", "format": "json", "size": 50})
    results = r.json().get("results", [])
    if not results:
        return None
    # proteomeType values returned by the current REST API (checked 2026-10-20): "Reference proteome",
    # "Non Reference proteome", "Redundant proteome", "Excluded". Older labels kept for safety.
    rank = {"Reference and representative proteome": 0, "Reference proteome": 0, "Representative proteome": 1,
            "Non Reference proteome": 2, "Other proteome": 2, "Redundant proteome": 3, "Excluded": 4}
    unknown = [p.get("proteomeType", "") for p in results if p.get("proteomeType", "") not in rank]
    if unknown:
        die(f"UniProt returned an unrecognised proteomeType {sorted(set(unknown))} for taxid {taxid}; update rank map")
    results.sort(key=lambda p: (rank[p["proteomeType"]], -int(p.get("proteinCount", 0) or 0)))
    return results[0]


def uniprot_search_tsv(query: str, fields: str, size: int = 25) -> list[dict]:
    r = http_get(f"{UNIPROT_REST}/uniprotkb/search",
                 params={"query": query, "fields": fields, "format": "tsv", "size": size})
    lines = r.text.strip().splitlines()
    if not lines:
        return []
    reader = csv.DictReader(lines, delimiter="\t")
    return list(reader)


def uniprot_fetch_fasta(accessions: list[str]) -> str:
    out = []
    for i in range(0, len(accessions), 100):
        chunk = accessions[i:i + 100]
        q = " OR ".join(f"accession:{a}" for a in chunk)
        r = http_get(f"{UNIPROT_REST}/uniprotkb/stream", params={"query": q, "format": "fasta"})
        out.append(r.text)
    return "".join(out)


# ---------------------------------------------------------------- Zenodo

def zenodo_files(record: int) -> list[dict]:
    r = http_get(f"{ZENODO_API}/{record}")
    files = r.json().get("files", [])
    if not files:
        die(f"Zenodo record {record}: no files listed")
    return [{"name": f["key"], "url": f["links"]["self"], "size": f.get("size"), "checksum": f.get("checksum")}
            for f in files]
