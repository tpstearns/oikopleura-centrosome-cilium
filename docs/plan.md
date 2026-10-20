# Plan of work

**Project description (TS, 2026-10-20).** An unusual feature of *Oikopleura* is that the genome is extremely reduced compared to other chordates. Our task is to do an in-depth analysis of the genes associated with proteins of the centrosome and cilium and determine what has been lost and what has been retained in Od and related species. This is of interest because the centrosome and particularly the cilium are subject to rapid evolutionary change in many groups, including complete loss in several branches of the eukaryotic tree. Goal: publish a paper. Strictly computational; outputs are a paper and resources.

## Objective
Determine which centrosome, centriole, and cilium genes *O. dioica* and related tunicates have lost and retained relative to other chordates, place those losses against the recurrent patterns of centriole/cilium reduction and loss across eukaryotes, and publish the analysis with its resource. The distinctive position: *Oikopleura* is not a cilia-loss lineage. It keeps sperm flagella, motile ciliated epithelia (gill ring, ciliated funnel, gut), and monociliated sensory cells (Nishida et al. 2021), so it tests reduction without loss: what a chordate must keep to build centrioles and cilia when gene loss is cheap.

## Species panel
- Core: *O. dioica* cryptic lineages (Okinawa OKI2018_I69 and the 2025 haplotype-resolved assembly; Osaka; Bergen; Barcelona) plus the legacy OdB3 annotation. Differences among them are loss in progress without morphological change.
- Related larvaceans: whichever of the Naville et al. 2019 genomes are public (to verify).
- Ascidians and thaliaceans via ANISEED/NCBI.
- Reference backbone: human, mouse, zebrafish, *Xenopus*, amphioxus; *Chlamydomonas*, *Tetrahymena* (complete toolkit); *Drosophila*, *C. elegans* (reduced/divergent); *Arabidopsis*, yeast (centriole-free).
- Eukaryote-wide backdrop imported from Hodges 2010, Carvalho-Santos 2010, Nevers 2017, Dobbelaere 2023 rather than recomputed.

## Gene set
~200 human-anchored genes by module (config/gene_modules.tsv): tubulins incl. ZED; duplication/cartwheel; distal end; PCM/γ-TuRC; appendages/basal foot; transition zone; IFT/BBSome/motors; motile axoneme and assembly factors; multiciliogenesis; ciliary signaling; linkers/satellites/IDPs; tubulin PTM enzymes; controls. Non-human anchors where human lacks the gene. TS curates.

## Pipeline (designed around false absences)
1. Orthogroups across the panel (OrthoFinder) for paralog structure.
2. Per-gene HMMs from vertebrate and invertebrate seeds (pyhmmer) against proteomes; miniprot against genomes for unannotated genes; transcripts where available.
3. RBH (MMseqs2) and Pfam architecture to confirm identity.
4. Structure-aware pass: Foldseek against AlphaFold models where they exist; ProstT5 3Di prediction for new assemblies.
5. Evidence tiers E1–E4 per cell (docs/evidence_tiers.md); sensitivity controls per species.
6. Loss timing by parsimony on the chordate/tunicate tree.

## Depth beyond presence/absence
- Retained core: conservation of known assembly interfaces (SAS-6 head/coiled-coil, STIL–PLK4, CEP135–SAS-6, CPAP–tubulin), domain truncations, paralog collapse; branch-length ratios of retained centriole genes versus the genome-wide distribution.
- Module-level coelimination: are losses clustered by module, and do they match the modules that go first elsewhere in eukaryotes (BBSome, ZED, PCM-specific proteins), or does compaction cut into the conserved core?
- Cilium types against modules: ZED and basal-foot status vs. predicted centriole architecture; multiciliogenesis genes vs. whether gill-ring, funnel, and gut cells are multiciliated (ultrastructure literature); sensory-cilium modules vs. monociliated sensory cells.
- Expression support from OikoBase developmental transcriptomes.

## Outputs
Manuscript (GBE / MBE / Current Biology depending on the biology); evidence-tiered matrix with per-gene evidence; curated module HMMs; pipeline (this repo); Zenodo deposit; browsable matrix page.

## Phases
1. Landscape and data inventory (in progress; novelty check clean so far).
2. Pipeline and matrix (bulk of the work; Claude Code on TS's MacBook).
3. Depth analyses and figures.
4. Drafting in TS's voice with verified references.

Division of labor: TS curates the gene list, adjudicates ambiguous calls, owns interpretation; Claude runs the pipeline, analyses, and drafts.

## Status (2026-10-20, evening)
- Environment: `oiko-cc` built natively on TS's MacBook (osx-arm64, no fallbacks); `snakemake -n` clean.
- Panel (`config/species_panel.yaml`, 37 entries): verified and fetched — 20 proteomes (6 tunicates, 6 chordate backbone, 6 non-chordate controls, plus Ciona type B and Thalia) and 6 larvacean genomes (Naville 2019; genome-only, no annotations exist). Still to-verify: Halocynthia, Molgula ×2, Botrylloides (ANISEED downloads not scriptable), Salpa thompsoni (no annotation). *O. dioica* core: Okinawa OKI2018 (OIST + RefSeq annotations, genome), OdB3 = Bergen (GenBank annotation, genome), and the four Plessy 2024 Zenodo genomes (Osaka OSKA2016v1.9, Barcelona Bar2_p4, Kume KUM-M3-7f, Aomori AOM-5-5f; lengths match Table 1 exactly) are fetched; their Augustus proteomes are fetched, as are Augustus proteomes for the six larvacean genomes (uncurated; miniprot remains mandatory). PRJEB91365 is reads only. Summary: `results/panel_summary.md`.
- Genes: 272/272 symbols resolved to anchors (`results/anchors.tsv`); TUBZ1 anchored on *X. tropicalis* tubz1 by accession.
- Seeds and profiles built for all 272 genes (`results/profiles_summary.md`): RBH recovery falls from 266/272 (mouse) to 162 (Ciona), 138 (Chlamydomonas), 85 (worm); eleven profiles rest on ≤3 seeds (TUBG2 anchor-only).
- Next: proteome searches on Od lineages and *C. robusta* (ZED and multiciliogenesis modules as first checkpoint); genome-only species need a tiers path without proteome evidence (miniprot-only), to add to the Snakefile.
- GitHub: commits are local; no GitHub credential on this Mac yet (no gh, SSH key, or keychain entry), so pushes wait for TS.
- Novelty: Hodges 2010 (45 genomes) does not include *Oikopleura*; Dobbelaere 2023 (147 genomes) tunicate inclusion to check (Dataset EV1A); no *Oikopleura*-focused survey found.
- Open for TS: panel breadth (four more RefSeq-annotated ascidians available); gene list curation.
