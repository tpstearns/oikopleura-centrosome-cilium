#!/usr/bin/env python3
"""Download proteomes, genomes, and annotations for one species, recording checksums in the manifest.

usage:
  fetch_data.py proteome --panel config/species_panel.yaml --species-tsv data/manifest/species.tsv \
                         --key Hsap --out data/proteomes/Hsap.faa --manifest data/manifest/files.tsv
  fetch_data.py genome   ... --out data/genomes/Odio_OKI.fna
  fetch_data.py pfam     --out data/pfam/Pfam-A.hmm --manifest data/manifest/files.tsv

Sources handled: uniprot_reference (reference proteome FASTA, one sequence per gene), zenodo (record + file
globs), zenodo_archive (record + archive file name + exact member path inside the tar; the archive is kept under
data/downloads/), ncbi_datasets (GCA_/GCF_ accession via the `datasets` CLI), url. Any other source value stops
the run.
"""
import argparse
import fnmatch
import gzip
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import (UNIPROT_FTP_REF, die, download, log, md5_of, read_tsv, record_file,  # noqa: E402
                    species_by_key, zenodo_files)

PFAM_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"


def gunzip_to(src: Path, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(src, "rb") as fin, open(dest, "wb") as fout:
        shutil.copyfileobj(fin, fout)


def fetch_uniprot_reference(upid: str, taxid: str, out: Path, manifest: Path) -> None:
    # Reference proteomes are grouped by superkingdom on the FTP site; try Eukaryota first, then the others.
    last_err = None
    for kingdom in ("Eukaryota", "Bacteria", "Archaea", "Viruses"):
        url = f"{UNIPROT_FTP_REF}/{kingdom}/{upid}/{upid}_{taxid}.fasta.gz"
        try:
            gz = download(url, out.with_suffix(".fasta.gz"))
        except SystemExit as exc:  # download() dies on failure; keep trying other kingdoms
            last_err = exc
            continue
        gunzip_to(gz, out)
        record_file(manifest, out, url, note=f"UniProt reference proteome {upid}, one sequence per gene")
        gz.unlink()
        return
    die(f"could not download reference proteome {upid} from UniProt FTP ({last_err})")


def fetch_zenodo(record: int, globs: list[str], out: Path, manifest: Path) -> None:
    files = zenodo_files(record)
    matches = [f for f in files if any(fnmatch.fnmatch(f["name"], g) for g in globs)]
    if len(matches) != 1:
        names = [f["name"] for f in files]
        die(f"Zenodo record {record}: expected exactly one file matching {globs}, got "
            f"{[m['name'] for m in matches]}. Files in record: {names}. Narrow the globs in the panel.")
    f = matches[0]
    dest = out.parent / f["name"]
    download(f["url"], dest)
    if f["name"].endswith(".gz"):
        gunzip_to(dest, out)
        dest.unlink()
    else:
        dest.replace(out)
    record_file(manifest, out, f"zenodo:{record}/{f['name']}", note=f"zenodo md5 {f.get('checksum')}")


def zenodo_archive(record: int, archive: str, manifest: Path) -> Path:
    """Download one archive from a Zenodo record into data/downloads (once) and record it; return its path."""
    files = {f["name"]: f for f in zenodo_files(record)}
    if archive not in files:
        die(f"Zenodo record {record}: no file named {archive!r}; files: {sorted(files)}")
    f = files[archive]
    dest = Path("data/downloads") / f"zenodo_{record}" / archive
    if dest.exists() and dest.stat().st_size == int(f["size"] or 0):
        log(f"  archive already present: {dest}")
    else:
        log(f"  downloading {archive} ({int(f['size'] or 0) / 1e6:.0f} MB) from Zenodo record {record}")
        download(f["url"], dest)
        expected = (f.get("checksum") or "").replace("md5:", "")
        if expected and md5_of(dest) != expected:
            die(f"{dest}: md5 {md5_of(dest)} does not match Zenodo checksum {expected}")
        record_file(manifest, dest, f"zenodo:{record}/{archive}", note=f"zenodo md5 {f.get('checksum')} verified")
    return dest


def fetch_zenodo_archive(record: int, archive: str, member: str, out: Path, manifest: Path) -> None:
    """Extract one member (exact path inside the tar) from a Zenodo archive into `out`."""
    import tarfile
    arc = zenodo_archive(record, archive, manifest)
    out.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(arc) as tf:
        try:
            info = tf.getmember(member)
        except KeyError:
            names = [n for n in tf.getnames() if not n.endswith("/")]
            die(f"{arc}: no member {member!r}. Members ({len(names)}): {names[:60]}{' ...' if len(names) > 60 else ''}")
        if not info.isfile():
            die(f"{arc}: member {member!r} is not a regular file")
        src = tf.extractfile(info)
        if member.endswith(".gz"):
            with gzip.open(src, "rb") as fin, open(out, "wb") as fout:
                shutil.copyfileobj(fin, fout)
        else:
            with open(out, "wb") as fout:
                shutil.copyfileobj(src, fout)
    record_file(manifest, out, f"zenodo:{record}/{archive}!{member}")


def fetch_ncbi_datasets(accession: str, kind: str, out: Path, manifest: Path) -> None:
    include = {"genome": "genome", "proteome": "protein", "annotation": "gff3"}[kind]
    with tempfile.TemporaryDirectory() as tmp:
        zip_path = Path(tmp) / "ncbi.zip"
        cmd = ["datasets", "download", "genome", "accession", accession, "--include", include,
               "--filename", str(zip_path)]
        log("  " + " ".join(cmd))
        subprocess.run(cmd, check=True)
        subprocess.run(["unzip", "-q", "-o", str(zip_path), "-d", tmp], check=True)
        suffix = {"genome": "*.fna", "proteome": "protein.faa", "annotation": "genomic.gff"}[kind]
        found = list(Path(tmp).rglob(suffix))
        if len(found) != 1:
            die(f"NCBI Datasets {accession}: expected one {suffix}, found {found}")
        out.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(found[0], out)
    record_file(manifest, out, f"ncbi_datasets:{accession}:{include}")


def fetch_url(url: str, out: Path, manifest: Path) -> None:
    dest = download(url, out if not url.endswith(".gz") else out.with_suffix(out.suffix + ".gz"))
    if url.endswith(".gz"):
        gunzip_to(dest, out)
        dest.unlink()
    record_file(manifest, out, url)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("kind", choices=["proteome", "genome", "annotation", "pfam"])
    ap.add_argument("--panel")
    ap.add_argument("--species-tsv")
    ap.add_argument("--key")
    ap.add_argument("--out", required=True)
    ap.add_argument("--manifest", required=True)
    args = ap.parse_args()
    out, manifest = Path(args.out), Path(args.manifest)

    if args.kind == "pfam":
        gz = download(PFAM_URL, out.with_suffix(".hmm.gz"))
        gunzip_to(gz, out)
        gz.unlink()
        record_file(manifest, out, PFAM_URL)
        return

    if not (args.panel and args.species_tsv and args.key):
        die("--panel, --species-tsv and --key are required for proteome/genome/annotation")
    sp = species_by_key(args.panel, args.key)
    srow = next((r for r in read_tsv(args.species_tsv) if r["key"] == args.key), None)
    if srow is None or srow["status"] == "unresolved":
        die(f"[{args.key}] not resolved in {args.species_tsv}; run panel_resolve first")
    spec = sp.get(args.kind) or {}
    source = spec.get("source")
    log(f"[{args.key}] {args.kind} via {source}")
    if source == "uniprot_reference":
        if not srow["upid"]:
            die(f"[{args.key}] no UniProt proteome id in species.tsv")
        fetch_uniprot_reference(srow["upid"], srow.get("proteome_taxid") or srow["taxid"], out, manifest)
    elif source == "zenodo":
        fetch_zenodo(int(spec["record"]), spec.get("files", ["*"]), out, manifest)
    elif source == "ncbi_datasets":
        acc = spec.get("accession", "")
        if not acc or acc == "TO-RESOLVE" or not acc.startswith(("GCA_", "GCF_")):
            die(f"[{args.key}] ncbi_datasets accession not set: {acc!r}")
        fetch_ncbi_datasets(acc, args.kind, out, manifest)
    elif source == "zenodo_archive":
        for k in ("record", "archive", "member"):
            if not spec.get(k):
                die(f"[{args.key}] zenodo_archive source needs {k}")
        fetch_zenodo_archive(int(spec["record"]), spec["archive"], spec["member"], out, manifest)
    elif source == "url":
        fetch_url(spec["url"], out, manifest)
    else:
        die(f"[{args.key}] {args.kind} source '{source}' is not fetchable; resolve it in the panel first")


if __name__ == "__main__":
    main()
