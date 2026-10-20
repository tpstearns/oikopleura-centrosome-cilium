# Seeds and per-gene profiles (2026-10-20, rebuilt with full-proteome reciprocity)

Seeds: MMseqs2 reciprocal best hits (E ≤ 1e-5, ≥50% anchor coverage, -s 7.5) of the whole anchor-species proteome against each of the ten seed proteomes (human, mouse, Xenopus, zebrafish, amphioxus, Ciona robusta, Chlamydomonas, Tetrahymena, fly, worm), kept where the reciprocal best is the anchor itself; the human anchor is always included, identical sequences are not added twice, and curated spurious RBHs are dropped (`config/seed_exclude.tsv`). Profiles: MAFFT --auto alignment, hmmbuild via pyhmmer. Table: `results/profiles_summary.tsv`.

272 profiles. Seeds per gene: 1: 1, 2: 6, 3: 10, 4: 28, 5: 43, 6: 48, 7: 31, 8: 38, 9: 37, 10: 30.

## RBH recovery per seed species (genes with an RBH / 272)

| species | genes with RBH |
|---|---|
| Mmus | 265 |
| Xtro | 256 |
| Drer | 244 |
| Bflo | 209 |
| Crob | 149 |
| Crei | 121 |
| Tthe | 121 |
| Dmel | 106 |
| Cele | 65 |

## Mean seeds per gene by module

| module | genes | mean seeds | genes with ≤3 seeds |
|---|---|---|---|
| A_tubulins | 7 | 6.1 | TUBB, TUBG2 |
| B_duplication_cartwheel | 24 | 5.6 | CEP57L1, CEP295 |
| C_distal_end_length | 11 | 6.4 | — |
| D_PCM_gTuRC | 18 | 5.7 | MZT2A, CDK5RAP2 |
| E_appendages_basal_foot | 19 | 5.4 | CCDC68 |
| F_transition_zone | 24 | 7.1 | RPGRIP1 |
| G_IFT_BBSome_motors | 48 | 8.4 | IFT70A, TULP3 |
| H_motile_axoneme | 51 | 7.8 | DNAH12 |
| I_multiciliogenesis | 14 | 4.9 | CDC20B, DEUP1 |
| J_ciliary_signaling | 21 | 5.1 | PTCH2, GLI1 |
| K_linkers_satellites_IDPs | 8 | 5.2 | ALMS1 |
| L_tubulin_PTM | 23 | 5.9 | AGBL3 |
| M_controls | 4 | 8.5 | — |

## Weak profiles (≤3 seeds; hits need RBH or architecture support, E2 at best otherwise)

| gene | module | seeds | species | profile length |
|---|---|---|---|---|
| TUBB | A_tubulins | 3 | Xtro, Drer | 444 |
| TUBG2 | A_tubulins | 1 | — | 451 |
| CEP57L1 | B_duplication_cartwheel | 3 | Mmus, Xtro | 455 |
| CEP295 | B_duplication_cartwheel | 2 | Mmus | 2621 |
| MZT2A | D_PCM_gTuRC | 2 | Crob | 167 |
| CDK5RAP2 | D_PCM_gTuRC | 2 | Mmus | 1911 |
| CCDC68 | E_appendages_basal_foot | 3 | Mmus, Xtro | 335 |
| RPGRIP1 | F_transition_zone | 3 | Mmus, Xtro | 1274 |
| IFT70A | G_IFT_BBSome_motors | 3 | Crei, Dmel | 650 |
| TULP3 | G_IFT_BBSome_motors | 2 | Mmus | 462 |
| DNAH12 | H_motile_axoneme | 3 | Mmus, Xtro | 3081 |
| CDC20B | I_multiciliogenesis | 3 | Mmus, Xtro | 519 |
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

Notes: TUBG2 is anchor-only because every species' γ-tubulin is reciprocal to TUBG1. Full-proteome reciprocity removed paralog seeds such as fly/worm/Chlamydomonas CDC20 (fzy) from CDC20B, whose profile now rests on mouse and Xenopus only; profiles of vertebrate-expanded families (TP73, E2F5, RFX2, KIF3C, TTLLs) are correspondingly vertebrate-biased, which is the paralog-group problem recorded in docs/decisions.md.
