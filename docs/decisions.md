# Decision log

Newest first. One line per decision with the reason; longer rationale in the relevant doc.

- 2026-10-20 — Environment: micromamba 2.9.0 (Homebrew), root prefix `~/micromamba`, env `oiko-cc` from `environment.yml`. Every package, including all bioconda tools (mmseqs2 18, hmmer 3.4, pyhmmer 0.12, miniprot 0.18, foldseek 10, orthofinder 3.1, diamond 2.2, iqtree 3.1, mafft 7.526, ncbi-datasets-cli 18), solved natively for osx-arm64; no `CONDA_SUBDIR=osx-64` or Docker fallback was needed. Activate with `eval "$(micromamba shell hook -s zsh)" && micromamba activate oiko-cc` (MAMBA_ROOT_PREFIX=~/micromamba). Snakemake 9.27; `snakemake -n` clean (347 jobs).
- 2026-10-20 — `panel_resolve` throttles NCBI E-utilities to 3 req/s (429s seen on the first run) and reports every panel failure in one pass instead of stopping at the first. Reason: faster iteration on the panel.
- 2026-10-20 — UniProt proteome selection fixed: the REST API now labels proteomes "Reference proteome"/"Non Reference proteome"; the old rank map ignored the type and picked the largest proteome (yeast resolved to a 12,520-protein non-reference strain). Reference proteomes are preferred; unknown labels fail the run. Panel entries may set `proteome_taxid` for reference proteomes attached to strain-level taxids (yeast S288c 559292 → UP000002311; Tetrahymena SB210 312017 → UP000009168; both checked live).
- 2026-10-20 — UniProt reference proteome FASTA files (`*_<taxid>.fasta.gz`, canonical sequences, release 2026_03) are the proteome source for the backbone and control species; isoform (`_additional`) files are not used.

- 2026-10-20 — Compute moves to Claude Code on TS's MacBook (M5 Max, 128 GB); this repo is the shared backbone; the claude.ai Project keeps planning, literature, and manuscript. Reason: the cloud sandbox has 2 CPUs, 7 GB RAM, a network allowlist that blocks ENA/NCBI/UniProt/Zenodo/AFDB, and no persistence across sessions.
- 2026-10-20 — Scope confirmed by TS: strictly computational; outputs are a paper and resources. Framing: what centrosome/cilium genes Od and related species have lost and retained, against the eukaryote-wide pattern of centriole/cilium reduction and loss.
- 2026-10-20 — Recommended panel: Od cryptic lineages as the core, other larvaceans and ascidians/thaliaceans as the comparative set, chordate backbone and non-chordate controls. Breadth still to be confirmed by TS.
- 2026-10-20 — Evidence tiers E1–E4 adopted (docs/evidence_tiers.md); genome-level search mandatory before any absence call in a core species.
- 2026-10-20 — Gene set is human-anchored with explicit non-human anchors where human lacks the gene (zeta-tubulin → Xenopus). Draft list in config/gene_modules.tsv awaits TS curation.
- 2026-10-20 — Correction recorded: Naville, Henriet et al. larvacean genome-size paper is Current Biology 2019 (doi:10.1016/j.cub.2019.01.080), not 2026.
