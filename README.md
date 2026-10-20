# Centrosome and cilium gene loss and retention in *Oikopleura dioica* and related tunicates

Computational project (T. Stearns, Rockefeller University). Goal: an in-depth, evidence-tiered analysis of which centrosome, centriole, and cilium genes have been lost and which retained in *Oikopleura dioica* (Od), its cryptic lineages, other larvaceans, and other tunicates, placed against the recurrent patterns of centriole/cilium reduction and loss across eukaryotes. Output: a paper and a reusable resource (curated gene modules, module HMMs, presence/absence matrix with evidence tiers, reproducible pipeline).

The working plan lives in `docs/plan.md`; the verified bibliography in `docs/bibliography.md`; decisions in `docs/decisions.md`. `CLAUDE.md` carries the standing instructions for Claude Code sessions in this repo.

## Layout

```
config/
  species_panel.yaml   # species, taxon IDs, assembly/proteome sources (verified vs. to-verify flagged)
  gene_modules.tsv     # curated gene set by module; edit freely (TS)
  params.yaml          # thresholds for searches and evidence tiers
workflow/Snakefile     # pipeline wiring
scripts/               # one script per stage; all take explicit paths
docs/                  # plan, bibliography, evidence-tier definitions, decision log
data/                  # downloaded data (git-ignored) + manifests with checksums (committed)
results/               # small outputs committed (matrix, summaries, per-gene evidence)
```

## Pipeline stages

1. `panel_resolve` — validate species names against NCBI Taxonomy and resolve UniProt reference proteomes by taxon ID; write `data/manifest/species.tsv`. Fails closed on any mismatch.
2. `fetch` — download proteomes (UniProt reference proteomes), genomes/annotations (NCBI Datasets, ENA, Zenodo), with SHA-256 recorded in `data/manifest/files.tsv`.
3. `genes_resolve` — map gene symbols in `config/gene_modules.tsv` to UniProt accessions and sequences for the anchor species (human by default; non-human anchors where a gene is absent in human). Fails on unresolved symbols.
4. `seeds` + `profiles` — collect ortholog seeds per gene from the reference backbone (reciprocal best hits, MMseqs2), align (MAFFT), build per-gene HMMs (HMMER/pyhmmer).
5. `search_proteomes` — hmmsearch every profile against every proteome; `rbh` — reciprocal best hits against anchors; `pfam` — domain architecture of candidate hits.
6. `search_genomes` — spliced protein-to-genome alignment (miniprot) of anchors and seeds against each genome, to catch genes missing from annotations.
7. `foldseek` (optional) — structure-aware confirmation with AlphaFold DB models or ProstT5 3Di prediction.
8. `tiers` + `matrix` — combine evidence into tiers E1–E4 per gene × species (see `docs/evidence_tiers.md`), build the matrix and module summaries.

## Setup (macOS arm64 or Linux)

```
mamba env create -f environment.yml
mamba activate oiko-cc
snakemake -n            # dry run; shows the stages and what each needs
snakemake --cores 12 panel_manifest
```

Large data stays under `data/` and out of git; manifests with checksums are committed so any run can be reproduced.

## Status

See `docs/plan.md` §Status and `docs/decisions.md`.
