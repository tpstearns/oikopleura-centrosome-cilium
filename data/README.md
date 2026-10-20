# data/

Everything here except `manifest/` is git-ignored. Downloaded files are recorded in `manifest/files.tsv` (path, source URL or accession, SHA-256, size, date). Species-level resolution (taxon ID, proteome ID, assembly accession) is in `manifest/species.tsv`, written by `scripts/panel_resolve.py`.

Expected subdirectories after `snakemake fetch`:

```
data/proteomes/<species_key>.faa      one protein per gene where the source provides it
data/genomes/<species_key>.fna        assembly FASTA
data/annotations/<species_key>.gff    annotation when available
data/transcripts/<species_key>.fna    transcript assemblies when available (OikoBase, Osaka)
data/genes/anchors.faa                anchor sequences for the gene set
data/seeds/<gene>.faa                 ortholog seeds per gene
data/profiles/<gene>.hmm              per-gene HMMs
data/pfam/Pfam-A.hmm                  Pfam-A (download separately; see scripts/fetch_data.py --pfam)
```
