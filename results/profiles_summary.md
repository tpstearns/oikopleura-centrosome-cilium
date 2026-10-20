# Seeds and per-gene profiles (2026-10-20)

Seeds: MMseqs2 reciprocal best hits (E ≤ 1e-5, ≥50% anchor coverage, -s 7.5) of each anchor against the ten seed proteomes (human, mouse, Xenopus, zebrafish, amphioxus, Ciona robusta, Chlamydomonas, Tetrahymena, fly, worm); the human anchor is always included, identical sequences are not added twice, and five curated spurious RBHs are dropped (`config/seed_exclude.tsv`). Profiles: MAFFT --auto alignment, hmmbuild via pyhmmer. Table: `results/profiles_summary.tsv`.

272 profiles. Seeds per gene: 1: 1, 2: 4, 3: 6, 4: 23, 5: 37, 6: 47, 7: 30, 8: 38, 9: 44, 10: 40, 11: 2.

## RBH recovery per seed species (genes with an RBH / 272)

| species | genes with RBH |
|---|---|
| Mmus | 266 |
| Xtro | 260 |
| Drer | 251 |
| Bflo | 223 |
| Crob | 162 |
| Crei | 136 |
| Tthe | 135 |
| Dmel | 123 |
| Cele | 83 |

## Mean seeds per gene by module

| module | genes | mean seeds | genes with ≤3 seeds |
|---|---|---|---|
| A_tubulins | 7 | 7.0 | TUBG2 |
| B_duplication_cartwheel | 24 | 5.8 | CEP57L1, CEP295 |
| C_distal_end_length | 11 | 6.7 | — |
| D_PCM_gTuRC | 18 | 5.9 | MZT2A |
| E_appendages_basal_foot | 19 | 5.5 | — |
| F_transition_zone | 24 | 7.2 | RPGRIP1 |
| G_IFT_BBSome_motors | 48 | 8.9 | — |
| H_motile_axoneme | 51 | 8.2 | DNAH12 |
| I_multiciliogenesis | 14 | 6.1 | DEUP1 |
| J_ciliary_signaling | 21 | 5.8 | PTCH2, GLI1 |
| K_linkers_satellites_IDPs | 8 | 5.4 | ALMS1 |
| L_tubulin_PTM | 23 | 5.9 | AGBL3 |
| M_controls | 4 | 10.0 | — |

## Weak profiles (≤3 seeds; hits need RBH or architecture support, E2 at best otherwise)

| gene | module | seeds | species | profile length |
|---|---|---|---|---|
| TUBG2 | A_tubulins | 1 | — | 451 |
| CEP57L1 | B_duplication_cartwheel | 3 | Mmus, Xtro | 455 |
| CEP295 | B_duplication_cartwheel | 2 | Mmus | 2621 |
| MZT2A | D_PCM_gTuRC | 2 | Crob | 167 |
| RPGRIP1 | F_transition_zone | 3 | Mmus, Xtro | 1274 |
| DNAH12 | H_motile_axoneme | 3 | Mmus, Xtro | 3081 |
| DEUP1 | I_multiciliogenesis | 3 | Mmus, Xtro | 603 |
| PTCH2 | J_ciliary_signaling | 2 | Mmus | 1203 |
| GLI1 | J_ciliary_signaling | 2 | Mmus | 1113 |
| ALMS1 | K_linkers_satellites_IDPs | 3 | Mmus, Xtro | 3448 |
| AGBL3 | L_tubulin_PTM | 3 | Mmus, Xtro | 1004 |

## Curated seed exclusions

| gene | species | accession | reason |
|---|---|---|---|
| AKAP9 | Dmel | A0A0B4KG35 | Strn-Mlck (stretchin-MLCK, 8255 aa): coiled-coil RBH, not an AKAP9 ortholog |
| GLI1 | Crei | A0A2K3CT35 | pherophorin-domain cell-wall glycoprotein (3003 aa); Chlamydomonas has no GLI/Hedgehog pathway |
| PTCH2 | Crei | A0A2K3E1H5 | SSD/RND-transporter family protein (PF00873+PF12349), not a Patched Hedgehog receptor; Chlamydomonas lacks Hedgehog signalling |
| RSPH1 | Cele | Q22667 | junctophilin-1: MORN-repeat match only; C. elegans has no motile cilia or radial spokes |
| TCHP | Cele | A0A1C3NSL9 | ajm-1 apical junction molecule: coiled-coil RBH, not a trichoplein ortholog |

Notes: TUBG2 is anchor-only because every species' RBH collapses onto TUBG1. Within-species duplicates (CCNO/Bflo, DNAAF6/Mmus, EEF1A1/Crei, IFT70A/Mmus, ODAD1/Xtro, ODAD2/Xtro, TTLL3/Xtro) are MMseqs2 score ties between near-identical paralogs or isoforms. Length outliers kept after inspection: DNAAF4, FAM161A, MKS1 (B9 domain), RSPH1 (MORN), INPP5E in Chlamydomonas; zebrafish hyls1 long isoform; amphioxus/Ciona SVBP-domain models; a Tetrahymena TTL-family protein for TTL (may be a TTLL; flagged for the PTM module).
