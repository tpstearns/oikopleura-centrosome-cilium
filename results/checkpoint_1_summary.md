# Checkpoint 1: proteome searches on the *O. dioica* lineages and *C. robusta* (2026-10-20)

Scope: hmmsearch with the 272 per-gene profiles and MMseqs2 reciprocal best hits (reciprocal against the full human proteome) on seven *O. dioica* proteomes (Okinawa OIST and RefSeq annotations, Osaka, Barcelona, Kume, Aomori, OdB3/Bergen) and *Ciona robusta* (KH). miniprot genome searches have been run on the six Od genomes and are used here only to annotate proteome misses. **Nothing below is an absence call**; E4 requires the full tier logic, the genome search, and the sensitivity controls (docs/evidence_tiers.md). Machine-readable tables: `checkpoint_proteome.tsv` (all cells), `checkpoint_proteome.md` (module counts + A, I, G detail), `checkpoint_OKI_best_human_hit.tsv`, `checkpoint_paralog_collisions.tsv`.

## Headline numbers
Genes with HMM and/or RBH evidence at the current thresholds, out of 272: Okinawa 168 (OIST) / 163 (RefSeq), Osaka 164, Barcelona 163, Kume 170, Aomori 170, OdB3 166; *Ciona* 172. A further 31–39 genes per Od lineage (57 in *Ciona*) have a strong profile hit (E ≤ 1e-10) with profile coverage below 0.5 ("HMM-partial"). The seven Od lineages agree closely; differences among them are dominated by annotation (the Augustus proteomes miss some genes that miniprot finds in the genome, e.g. ε-tubulin in Osaka/Barcelona/Aomori, δ-tubulin in Barcelona).

## ZED module (module A)
| gene | *O. dioica* (7 proteomes, 6 genomes) | *C. robusta* |
|---|---|---|
| TUBD1 (δ) | present: HMM+RBH in 6/7 proteomes (Barcelona Augustus set lacks the protein; miniprot finds the locus, PAR:9,498,952–9,500,282, 34% identity to human over 99%) | present (HMM+RBH) |
| TUBE1 (ε) | present: HMM+RBH in Okinawa (both annotations), Kume, OdB3; Osaka/Barcelona/Aomori proteomes carry only a β-tubulin cross-hit, but miniprot finds the ε locus in all six genomes (49–52% identity to *Ciona* ε over 94%) | present (HMM+RBH) |
| TUBZ1 (ζ) | **no distinct protein and no distinct locus.** The ζ profile's best proteome hit is the δ-tubulin (Okinawa, Kume, OdB3) or a β-tubulin (Osaka, Barcelona, Aomori) in every lineage; miniprot aligns the ζ queries (Xenopus anchor plus amphioxus, *Ciona*, *Chlamydomonas*, *Tetrahymena* seeds) only inside the ε locus at 29–31% identity over 39% of the query. | present: a third protein distinct from δ and ε, HMM+RBH, E = 3e-175 |
| TUBG1/TUBG2 | one γ-tubulin (reciprocal to TUBG2 in the full human proteome, owned by the TUBG1 profile); TUBG1/2 is a vertebrate duplication | one γ-tubulin |

Preliminary reading: *O. dioica* keeps δ- and ε-tubulin and has lost ζ-tubulin, while *Ciona* keeps all three. This is the pattern Turk et al. 2015 predicted for the ZED module (ζ can be lost while δ/ε persist). It needs (i) a tubulin tree placing the Od δ and ε candidates and the *Ciona* ζ with the vertebrate δ/ε/ζ clades, and (ii) a six-frame/loosened miniprot pass before ζ is reported as E4.

## Multiciliogenesis module (module I)
- **Vertebrate-specific regulators behave as expected and are uninformative about Od**: MCIDAS, GMNC and DEUP1 have no hit in any Od proteome *or in Ciona*; their profiles are built from vertebrate seeds only (Mmus, Xtro, Drer), i.e. no invertebrate ortholog was found anywhere in the seed panel. DEUP1 is the CEP63 paralog; CEP63 itself is present in Od.
- **CCDC78** (deuterosome/centriole-amplification protein): HMM+RBH in *Ciona* (E = 4e-229), no proteome hit and no miniprot hit in any Od lineage. Candidate loss; first gene for the full absence workflow.
- **FOXJ1**: present. The forkhead domain hits at E ≈ 1e-38 in every lineage but covers only 33–40% of the 421-aa profile, so it is "HMM-partial" at the current 50% coverage threshold. Threshold decision below.
- **RFX2/RFX3**: one RFX2/3-type gene in Od (credited to the RFX3 profile; the RFX2 reading is a cross-hit). RFX2/3 is a vertebrate duplication.
- **E2F4/E2F5**: one E2F4/5-type gene (owned by E2F4). GMNN present. CCNO present (HMM+RBH in 6/7). MYB present.
- **CDC20B**: the Od and *Ciona* hits are CDC20/FZR1 (best human hit outside the gene list); CDC20B itself is a vertebrate paralog. The CDC20B profile now rests on mouse and Xenopus seeds only.
- **TP73** (p53 family): no profile hit in Od, but a Pfam P53-domain search finds exactly one protein per Od proteome at E ≈ 1e-5–1e-6 (versus 1e-52–1e-58 for the two *Ciona* p53-family proteins). Od has a single, extremely divergent p53-family gene; the vertebrate-seeded TP73 profile cannot see it. A good test case for the Foldseek/structure-aware pass.

## IFT and BBSome (module G), noted in passing
IFT-A and IFT-B are intact in all Od lineages (IFT25 and IFT43 are the only "none"; both are short proteins, IFT43 has near-misses at E ≈ 1e-5). The BBSome core (BBS1, 2, 4, 5, 7, 8, ARL6, LZTFL1) is present; BBS9 is partial; BBIP1 (92 aa) is not found. The chaperonin-like BBS6/MKKS, BBS10 and BBS12 are not found: their profile "hits" are CCT subunits by best human hit, and the same holds in *Ciona*, consistent with these being vertebrate-specific. Retrograde dynein-2 is present (DYNC2H1 strong; the WD40 intermediate chains DYNC2I1/I2 and DYNC2LI1 are partial/HMM-only). KIF3A and KIF3B present; KIF3C (vertebrate paralog of KIF3B) reads as a cross-hit.

## Methodological findings that change the pipeline (all implemented, in docs/decisions.md)
1. **Reciprocity must be against the full anchor-species proteome.** The first RBH implementation was reciprocal only among the 272 anchors, so paralogs outside the list (CDC20 for CDC20B; fly fzy seeding the CDC20B profile) passed as the listed gene. Fixed in `rbh.py` and `seeds_collect.py`; seeds, profiles and searches were rebuilt. Okinawa RBH count went from 210 to 155; of the HMM-only genes that lost RBH, 27 have a best human hit outside the gene list and 6 a listed paralog (`checkpoint_OKI_best_human_hit.tsv`).
2. **A target protein is credited only to the profile that scores it best** (owner rule), now in the checkpoint and in `tiers.py`. Without it β-tubulins registered as ε, the δ-tubulin as ζ, CCTs as BBS10.
3. **Spurious RBH seeds** (coiled-coil/repeat matches) are excluded by a curated list with reasons (`config/seed_exclude.tsv`).

## Decisions needed from TS
1. **Paralog groups.** Vertebrate duplications (TUBG1/2, RFX2/3, E2F4/5, KIF3B/C, VASH1/2, PKD2/PKD2L1, TCTN1/3, TTLL subfamilies, DNAH7/12, DNAH9/11, AGBL, RPGRIP1/1L, CEP250/CROCC, GLI1/2/3, SHH/IHH/DHH, PTCH1/2, TUBGCP2/3/4/5, CEP57/57L1, CEP63/DEUP1, MZT2A/B, BBS6/10/12) each resolve to one Od gene. The matrix should score these as groups (one presence for the group, with the vertebrate paralog count recorded) or every such pair will inflate the loss count. `checkpoint_paralog_collisions.tsv` lists the 23 collisions detected automatically; a `group` column in `config/gene_modules.tsv` is the proposed mechanism, and group-level profiles (seeds pooled across the group) would also cure the vertebrate bias of profiles such as TP73 and CDC20B.
2. **Profile-coverage threshold.** 50% profile coverage fails single-domain transcription factors (FOXJ1 at E 1e-38) and long coiled-coil proteins (CEP192, PCM1, CDK5RAP2, CEP128, NPHP3 reach E ≤ 1e-40 in *Ciona* with 10–40% coverage). Proposal: keep 50% for E1; admit E ≤ 1e-10 with ≥ 10% coverage as E2 ("homolog present, architecture not established"). Only `config/params.yaml` and the decision log change.
3. **Best-human-hit as a standard column.** The cheap reverse search (every candidate vs the whole human proteome) resolved most ambiguities above; it should run in the RBH stage and feed `tiers.py` (a candidate whose best human hit is an unlisted paralog is not presence of the listed gene).
4. Whether to add CDC20, FZR1, CCNA/CCNB, E2F1-3, RFX1, TUBA/TUBB isotypes and the CCT subunits as explicit decoys in the gene list (alternative to 3).

## Next steps (in order)
Tubulin tree for the ZED call; CCDC78 and ζ-tubulin through the full absence workflow (looser miniprot, six-frame, transcript evidence where available); implement best-human-hit in the RBH stage; apply TS's decisions on groups and thresholds; then run the remaining 25 proteomes and the larvacean genomes and build the first matrix.
