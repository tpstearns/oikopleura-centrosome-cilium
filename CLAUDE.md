# Instructions for Claude Code sessions in this repo

## What this is
Comparative-genomic analysis of centrosome/centriole/cilium genes in *Oikopleura dioica* (Od), its cryptic lineages, other larvaceans, other tunicates, and a reference backbone of chordates and non-chordate eukaryotes. Output: a paper and a resource (gene modules, HMMs, evidence-tiered matrix, pipeline). PI: Tim Stearns (senior centrosome/cilium biologist; expert reader; lead with the answer, no preamble, no hedging beyond what changes a decision).

Read `docs/plan.md` first, then `docs/decisions.md`. Update both when something durable changes (decisions, scope, thresholds, data sources). Keep `docs/plan.md` §Status current at the end of each session.

## Hard rules
1. **No fabricated references, DOIs, accessions, or quotations.** Every citation and accession is verified against a live source before it enters `docs/bibliography.md`, `config/species_panel.yaml`, or any manuscript text. If it cannot be verified, say so and mark it `to-verify`.
2. **Absence must be earned.** A gene is called absent (tier E4) only after the proteome search, the genome search (miniprot), and any transcript search all fail at the thresholds in `config/params.yaml`, and the sensitivity controls for that species pass. False absences from divergent sequence and incomplete annotation are the main risk in this organism.
3. **Large data never goes into git.** `data/` is git-ignored except manifests. Record every downloaded file with URL and SHA-256 in `data/manifest/files.tsv`.
4. **Every call is traceable.** Each matrix cell links to its evidence (hit tables under `results/`). Do not hand-edit result files; change code or config and rerun.
5. **Scripts take explicit paths and fail loudly.** No silent fallbacks, no guessed taxon IDs or accessions in code.
6. **Commit messages describe the scientific change**, not the mechanics ("add miniprot genome search for Okinawa lineages", not "update script").

## Working conventions
- Environment: `mamba env create -f environment.yml && mamba activate oiko-cc`. On Apple silicon, if a bioconda package lacks an osx-arm64 build, use `CONDA_SUBDIR=osx-64` for that package or a Docker fallback; note it in `docs/decisions.md`.
- Pipeline: Snakemake (`workflow/Snakefile`). Run stages by target name; keep `snakemake -n` clean.
- Thresholds live in `config/params.yaml`, never inline in scripts.
- Gene set: `config/gene_modules.tsv` is curated by TS. Do not add or remove genes without recording why in `docs/decisions.md`.
- Species: `config/species_panel.yaml`. Entries flagged `status: to-verify` must be verified (taxon ID, assembly accession, URL) before use.
- Write results summaries to `results/` as TSV + a short Markdown summary. Figures go to `results/figures/`.

## Current priorities (update as work proceeds)
1. Verify and complete `config/species_panel.yaml`: pull accessions for the Osaka, Bergen, Barcelona, Aomori, Kume assemblies from Plessy et al. 2024 Supplementary Table 1; OdB3 (Genoscope/OikoBase) download location; the 2025 haplotype-resolved Okinawa assembly (ENA PRJEB91365); larvacean genomes from Naville et al. 2019; ascidian proteomes from ANISEED/NCBI.
2. Run `panel_manifest` and `fetch`; confirm checksums and proteome sizes look sane (Od proteomes ~16–19k proteins).
3. Run `genes_resolve`; send TS the list of unresolved or ambiguous symbols.
4. Build profiles; run proteome searches on the Od lineages and *Ciona robusta* first; report the ZED-module and multiciliogenesis-module results as the first checkpoint.
5. Genome searches (miniprot) for every E4 candidate before any absence is reported.
