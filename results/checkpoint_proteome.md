# Proteome-search checkpoint

Proteome evidence only (hmmsearch with the per-gene profiles + MMseqs2 RBH to the human anchors, thresholds from `config/params.yaml`). `none` means no proteome hit at thresholds and is **not** an absence call: the genome search decides. Table: `results/checkpoint_proteome.tsv`.

Classes: **HMM+RBH** both; **HMM** profile only; **RBH** reciprocal best hit only (reciprocal against the full anchor-species proteome); **HMM-partial** profile hit with E ≤ 1e-10 but profile coverage below 0.5; **none** no proteome hit (genome hit noted where miniprot found one).

## Genes with proteome evidence: HMM+RBH / HMM or RBH only / HMM-partial, out of module size

| module | n | Odio_OKI | Odio_OKI_refseq | Odio_OSA | Odio_BCN | Odio_KUM | Odio_AOM | Odio_OdB3 | Crob |
|---|---|---|---|---|---|---|---|---|---|
| A_tubulins | 7 | 2 / 4 / 0 | 2 / 4 / 0 | 1 / 4 / 1 | 0 / 4 / 1 | 2 / 4 / 0 | 1 / 4 / 1 | 3 / 3 / 0 | 5 / 1 / 0 |
| B_duplication_cartwheel | 24 | 7 / 4 / 3 | 6 / 2 / 4 | 5 / 5 / 3 | 7 / 3 / 2 | 9 / 3 / 1 | 8 / 1 / 2 | 7 / 4 / 2 | 9 / 3 / 8 |
| C_distal_end_length | 11 | 2 / 4 / 3 | 3 / 3 / 3 | 1 / 5 / 3 | 2 / 4 / 3 | 4 / 2 / 3 | 1 / 3 / 5 | 4 / 2 / 3 | 7 / 1 / 2 |
| D_PCM_gTuRC | 18 | 3 / 4 / 5 | 4 / 3 / 3 | 4 / 3 / 4 | 4 / 4 / 4 | 4 / 3 / 5 | 4 / 4 / 2 | 4 / 4 / 4 | 9 / 0 / 3 |
| E_appendages_basal_foot | 19 | 4 / 5 / 4 | 4 / 6 / 3 | 4 / 5 / 3 | 4 / 5 / 4 | 6 / 3 / 4 | 4 / 6 / 3 | 2 / 5 / 2 | 6 / 0 / 7 |
| F_transition_zone | 24 | 7 / 8 / 4 | 8 / 3 / 8 | 8 / 5 / 5 | 10 / 3 / 5 | 9 / 5 / 4 | 9 / 8 / 4 | 8 / 7 / 4 | 15 / 0 / 4 |
| G_IFT_BBSome_motors | 48 | 27 / 13 / 2 | 30 / 8 / 4 | 25 / 14 / 4 | 27 / 11 / 4 | 30 / 10 / 2 | 29 / 11 / 3 | 29 / 9 / 4 | 34 / 6 / 4 |
| H_motile_axoneme | 51 | 41 / 4 / 3 | 44 / 2 / 2 | 38 / 7 / 4 | 39 / 7 / 1 | 41 / 5 / 3 | 41 / 5 / 1 | 42 / 3 / 4 | 38 / 4 / 9 |
| I_multiciliogenesis | 14 | 2 / 5 / 2 | 3 / 3 / 2 | 2 / 4 / 1 | 2 / 5 / 0 | 2 / 3 / 3 | 2 / 4 / 1 | 1 / 5 / 2 | 4 / 5 / 2 |
| J_ciliary_signaling | 21 | 0 / 6 / 7 | 2 / 5 / 6 | 1 / 5 / 7 | 2 / 4 / 7 | 1 / 5 / 7 | 1 / 6 / 5 | 2 / 5 / 6 | 7 / 3 / 7 |
| K_linkers_satellites_IDPs | 8 | 2 / 1 / 0 | 2 / 0 / 1 | 1 / 2 / 0 | 2 / 0 / 1 | 2 / 0 / 1 | 2 / 1 / 0 | 1 / 1 / 0 | 2 / 0 / 4 |
| L_tubulin_PTM | 23 | 7 / 2 / 6 | 10 / 2 / 3 | 7 / 4 / 4 | 8 / 2 / 5 | 8 / 5 / 3 | 9 / 2 / 4 | 6 / 5 / 5 | 9 / 0 / 7 |
| M_controls | 4 | 4 / 0 / 0 | 4 / 0 / 0 | 4 / 0 / 0 | 4 / 0 / 0 | 4 / 0 / 0 | 4 / 0 / 0 | 4 / 0 / 0 | 4 / 0 / 0 |

Totals: genes with HMM/RBH evidence at thresholds (+ HMM-partial) / 272: Odio_OKI 168 (+39), Odio_OKI_refseq 163 (+39), Odio_OSA 164 (+39), Odio_BCN 163 (+37), Odio_KUM 170 (+36), Odio_AOM 170 (+31), Odio_OdB3 166 (+36), Crob 172 (+57).

## A_tubulins: per-gene detail

● HMM+RBH  ◐ HMM only  ◑ RBH only  ◔ HMM-partial (strong E, low profile coverage)  ○ no proteome hit (genome search pending); G = miniprot genome hit at thresholds

| gene | Odio_OKI | Odio_OKI_refseq | Odio_OSA | Odio_BCN | Odio_KUM | Odio_AOM | Odio_OdB3 | Crob |
|---|---|---|---|---|---|---|---|---|
| TUBA1A | ◐ 0 | ◐ 0 | ◐ 9.88e-171 | ◐ 0 | ◐ 0 | ◐ 1.33e-132 | ◐ 0 | ● 0 |
| TUBB | ◐ 1.11e-318 | ◐ 1.06e-318 | ◐ 1.32e-318 | ◐ 9.17e-319 | ◐ 9.99e-319 | ◐ 1.28e-318 | ◐ 1.1e-318 | ◐ 0 |
| TUBG1 | ◐ 2.08e-215 | ◐ 1.99e-215 | ◐ 4.56e-216 | ◐ 3.03e-216 | ◐ 1.87e-215 | ◐ 4.41e-216 | ● 2.2e-216 | ● 1.06e-272 |
| TUBG2 | ◑ 1.503E-198 | ◑ 1.384E-198 | ◑ 5.439E-199 | ◑ 2.385E-199 | ◑ 1.349E-198 | ◑ 5.539E-199 | ◑ 1.182E-198 | ○ (6.04e-264, cov 0.996) |
| TUBD1 | ● 3.2e-103 | ● 3.05e-103 | ● 2.2e-87 | ○ G | ● 8.57e-103 | ● 1.11e-87 | ● 2.92e-102 | ● 1.33e-212 |
| TUBE1 | ● 6.89e-109 | ● 4.43e-183 | ◔ E=4.6e-45 profile cov=0.258 | ◔ E=5.63e-99 profile cov=0.492 | ● 8.53e-134 | ◔ E=1.97e-45 profile cov=0.286 | ● 6.38e-178 | ● 3.44e-247 |
| TUBZ1 | ○ (2.04e-57, cov 0.979) | ○ (1.95e-57, cov 0.979) | ○ (1.4e-55, cov 0.988) | ○ (7.1e-57, cov 0.988) | ○ (1.01e-57, cov 0.976) | ○ (7.51e-57, cov 0.988) | ○ (1.38e-58, cov 0.986) | ● 2.67e-175 |

Notes:
- TUBG2 / Odio_OKI: cross-hit only: CAG5105093.1 is best explained by TUBG1 (E=5.82e-201); RBH target CAG5105093.1 is best explained by profile TUBG1
- TUBG2 / Odio_OKI_refseq: cross-hit only: XP_081645013.1 is best explained by TUBG1 (E=5.55e-201); RBH target XP_081645013.1 is best explained by profile TUBG1
- TUBG2 / Odio_OSA: cross-hit only: OSKA2016v1.9.g2032.t1 is best explained by TUBG1 (E=2.21e-201); RBH target OSKA2016v1.9.g2032.t1 is best explained by profile TUBG1
- TUBG2 / Odio_BCN: cross-hit only: Bar2_p4_Flye.g2051.t1 is best explained by TUBG1 (E=7.67e-202); RBH target Bar2_p4_Flye.g2051.t1 is best explained by profile TUBG1
- TUBG2 / Odio_KUM: cross-hit only: KUM-M3-7f.g16436.t1 is best explained by TUBG1 (E=5.23e-201); RBH target KUM-M3-7f.g16436.t1 is best explained by profile TUBG1
- TUBG2 / Odio_AOM: cross-hit only: AOM-5-5f.g5910.t1 is best explained by TUBG1 (E=2.14e-201); RBH target AOM-5-5f.g5910.t1 is best explained by profile TUBG1
- TUBG2 / Odio_OdB3: cross-hit only: CBY21563.1 is best explained by TUBG1 (E=6.32e-201); RBH target CBY21563.1 is best explained by profile TUBG1
- TUBG2 / Crob: cross-hit only: tr|H2XQX6|H2XQX6_CIOIN is best explained by TUBG1 (E=6.04e-264); near-miss E=6.04e-264 cov=0.996
- TUBD1 / Odio_BCN: cross-hit only: Bar2_p4_Flye.g12553.t1 is best explained by TUBB (E=1.62e-63); near-miss E=1.62e-63 cov=0.936; genome hit (miniprot)
- TUBE1 / Odio_OSA: partial: E=4.6e-45 profile cov=0.258; cross-hit only: OSKA2016v1.9.g422.t1 is best explained by TUBB (E=5.73e-97)
- TUBE1 / Odio_BCN: partial: E=5.63e-99 profile cov=0.492; cross-hit only: Bar2_p4_Flye.g11071.t1 is best explained by TUBB (E=3.6e-97)
- TUBE1 / Odio_AOM: partial: E=1.97e-45 profile cov=0.286; cross-hit only: AOM-5-5f.g9098.t1 is best explained by TUBB (E=3.7e-97)
- TUBZ1 / Odio_OKI: cross-hit only: CAG5091327.1 is best explained by TUBD1 (E=2.04e-57); near-miss E=2.04e-57 cov=0.979
- TUBZ1 / Odio_OKI_refseq: cross-hit only: XP_081636851.1 is best explained by TUBD1 (E=1.95e-57); near-miss E=1.95e-57 cov=0.979
- TUBZ1 / Odio_OSA: cross-hit only: OSKA2016v1.9.g10010.t1 is best explained by TUBB (E=1.4e-55); near-miss E=1.4e-55 cov=0.988
- TUBZ1 / Odio_BCN: cross-hit only: Bar2_p4_Flye.g11071.t1 is best explained by TUBB (E=7.1e-57); near-miss E=7.1e-57 cov=0.988
- TUBZ1 / Odio_KUM: cross-hit only: KUM-M3-7f.g6434.t1 is best explained by TUBD1 (E=1.01e-57); near-miss E=1.01e-57 cov=0.976
- TUBZ1 / Odio_AOM: cross-hit only: AOM-5-5f.g9098.t1 is best explained by TUBB (E=7.51e-57); near-miss E=7.51e-57 cov=0.988
- TUBZ1 / Odio_OdB3: cross-hit only: CBY19037.1 is best explained by TUBD1 (E=1.38e-58); near-miss E=1.38e-58 cov=0.986

## I_multiciliogenesis: per-gene detail

● HMM+RBH  ◐ HMM only  ◑ RBH only  ◔ HMM-partial (strong E, low profile coverage)  ○ no proteome hit (genome search pending); G = miniprot genome hit at thresholds

| gene | Odio_OKI | Odio_OKI_refseq | Odio_OSA | Odio_BCN | Odio_KUM | Odio_AOM | Odio_OdB3 | Crob |
|---|---|---|---|---|---|---|---|---|
| MCIDAS | ○ (1.15e-05, cov 0.159) | ○ (0.000638, cov 0.149) | ○ (4.11e-05, cov 0.191) | ○ (1.42e-05, cov 0.196) | ○ (1.03e-05, cov 0.159) | ○ (0.000304, cov 0.202) | ○ (1.7e-05, cov 0.196) | ○ (1.96e-08, cov 0.260) |
| GMNC | ○ | ○ | ○ | ○ | ○ | ○ | ○ | ○ |
| GMNN | ◐ 8.13e-11 | ● 7.75e-11 | ◑ 1.499E-07 | ● 7.27e-11 | ● 7.31e-11 | ◑ 1.526E-07 | ◐ 6.02e-11 | ● 3.17e-80 |
| E2F4 | ◐ 1.76e-41 | ◐ 2.63e-43 | ◐ 3.99e-43 | ◐ 1.19e-43 | ◐ 2.48e-43 | ◐ 1.1e-40 | ◐ 1.32e-43 | ◐ 1.31e-61 |
| E2F5 | ◔ E=1.92e-13 profile cov=0.281 | ◔ E=1.83e-13 profile cov=0.281 | ○ G | ○ G | ◔ E=3.34e-13 profile cov=0.281 | ○ G | ◔ E=4.12e-12 profile cov=0.406 | ● 2.83e-130 |
| CCNO | ● 1.14e-36 | ● 1.13e-36 | ● 1.21e-34 | ● 6.31e-35 | ◐ 1.06e-36 | ● 1.17e-34 | ● 6.79e-35 | ◐ 1.18e-43 |
| CDC20B | ◐ 2.19e-77 | ◐ 2.09e-77 | ◐ 1.81e-77 | ◐ 1.64e-76 | ◐ 1.97e-77 | ◐ 1.75e-77 | ◐ 1.97e-77 | ◐ 1.48e-91 |
| DEUP1 | ○ | ○ | ○ | ○ | ○ | ○ | ○ | ○ |
| CCDC78 | ○ | ○ | ○ | ○ | ○ | ○ | ○ | ● 3.59e-229 |
| FOXJ1 | ◔ E=2.48e-38 profile cov=0.333 | ◔ E=1.9e-38 profile cov=0.333 | ◔ E=4.17e-38 profile cov=0.396 | ◐ 3.98e-36 | ◔ E=2.23e-38 profile cov=0.333 | ◔ E=3.61e-38 profile cov=0.401 | ◐ 5.07e-36 | ◐ 5.48e-37 |
| RFX2 | ◑ 3.832E-100 | ○ (9.14e-124, cov 0.644) | ○ G | ○ G | ○ G | ○ G | ○ G | ◔ E=3.31e-23 profile cov=0.241 |
| RFX3 | ◐ 9.48e-128 | ◐ 4.37e-136 | ◐ 2.06e-134 | ◐ 3.87e-136 | ◔ E=2.98e-89 profile cov=0.450 | ◐ 1.02e-136 | ◐ 4.32e-136 | ● 0 |
| MYB | ● 1.37e-89 | ● 1.31e-89 | ● 6e-91 | ◑ 2.572E-73 | ● 6.1e-90 | ● 4.04e-66 | ◔ E=3.65e-29 profile cov=0.209 | ◔ E=1.71e-116 profile cov=0.464 |
| TP73 | ○ (0.000561, cov 0.150) | ○ (0.000579, cov 0.150) | ○ | ○ | ○ (0.000528, cov 0.150) | ○ | ○ | ◐ 4.25e-64 |

Notes:
- MCIDAS / Odio_OKI: near-miss E=1.15e-05 cov=0.159
- MCIDAS / Odio_OKI_refseq: near-miss E=0.000638 cov=0.149
- MCIDAS / Odio_OSA: near-miss E=4.11e-05 cov=0.191
- MCIDAS / Odio_BCN: near-miss E=1.42e-05 cov=0.196
- MCIDAS / Odio_KUM: near-miss E=1.03e-05 cov=0.159
- MCIDAS / Odio_AOM: near-miss E=0.000304 cov=0.202
- MCIDAS / Odio_OdB3: near-miss E=1.7e-05 cov=0.196
- MCIDAS / Crob: near-miss E=1.96e-08 cov=0.260
- GMNN / Odio_OKI_refseq: HMM best XP_081638445.1 != RBH XP_081644364.1
- GMNN / Odio_BCN: HMM best Bar2_p4_Flye.g9547.t1 != RBH Bar2_p4_Flye.g10490.t1
- GMNN / Odio_KUM: HMM best KUM-M3-7f.g5770.t1 != RBH KUM-M3-7f.g11708.t1
- E2F5 / Odio_OKI: partial: E=1.92e-13 profile cov=0.281; cross-hit only: CAG5113812.1 is best explained by E2F4 (E=3.58e-37)
- E2F5 / Odio_OKI_refseq: partial: E=1.83e-13 profile cov=0.281; cross-hit only: XP_081630492.1 is best explained by E2F4 (E=4.99e-41)
- E2F5 / Odio_OSA: cross-hit only: OSKA2016v1.9.g6854.t1 is best explained by E2F4 (E=4.31e-41); near-miss E=4.31e-41 cov=0.777; genome hit (miniprot)
- E2F5 / Odio_BCN: cross-hit only: Bar2_p4_Flye.g5509.t1 is best explained by E2F4 (E=1.21e-40); near-miss E=1.21e-40 cov=0.733; genome hit (miniprot)
- E2F5 / Odio_KUM: partial: E=3.34e-13 profile cov=0.281; cross-hit only: KUM-M3-7f.g14781.t1 is best explained by E2F4 (E=4.71e-41)
- E2F5 / Odio_AOM: cross-hit only: AOM-5-5f.g4831.t1 is best explained by E2F4 (E=1.99e-39); near-miss E=1.99e-39 cov=0.716; genome hit (miniprot)
- E2F5 / Odio_OdB3: partial: E=4.12e-12 profile cov=0.406; cross-hit only: CBY20531.1 is best explained by E2F4 (E=1.64e-40)
- CCNO / Odio_OKI: HMM best CAG5112868.1 != RBH CAG5105116.1
- CCNO / Odio_OKI_refseq: HMM best XP_081633323.1 != RBH XP_081645036.1
- CCNO / Odio_OSA: HMM best OSKA2016v1.9.g4637.t2 != RBH OSKA2016v1.9.g2043.t1
- CCNO / Odio_BCN: HMM best Bar2_p4_Flye.g4888.t2 != RBH Bar2_p4_Flye.g2039.t1
- CCNO / Odio_AOM: HMM best AOM-5-5f.g13025.t3 != RBH AOM-5-5f.g5922.t1
- CCNO / Odio_OdB3: HMM best CBY22811.1 != RBH CBY21542.1
- FOXJ1 / Odio_OKI: partial: E=2.48e-38 profile cov=0.333
- FOXJ1 / Odio_OKI_refseq: partial: E=1.9e-38 profile cov=0.333
- FOXJ1 / Odio_OSA: partial: E=4.17e-38 profile cov=0.396
- FOXJ1 / Odio_KUM: partial: E=2.23e-38 profile cov=0.333
- FOXJ1 / Odio_AOM: partial: E=3.61e-38 profile cov=0.401
- RFX2 / Odio_OKI: cross-hit only: CAG5087566.1 is best explained by RFX3 (E=3.21e-120); RBH target CAG5087566.1 is best explained by profile RFX3
- RFX2 / Odio_OKI_refseq: cross-hit only: XP_081635138.1 is best explained by RFX3 (E=9.14e-124); near-miss E=9.14e-124 cov=0.644
- RFX2 / Odio_OSA: cross-hit only: OSKA2016v1.9.g9596.t2 is best explained by RFX3 (E=5.23e-122); near-miss E=5.23e-122 cov=0.719; genome hit (miniprot)
- RFX2 / Odio_BCN: cross-hit only: Bar2_p4_Flye.g8288.t2 is best explained by RFX3 (E=1.73e-123); near-miss E=1.73e-123 cov=0.638; genome hit (miniprot)
- RFX2 / Odio_KUM: cross-hit only: KUM-M3-7f.g5137.t1 is best explained by RFX3 (E=8.47e-78); near-miss E=8.47e-78 cov=0.525; genome hit (miniprot)
- RFX2 / Odio_AOM: cross-hit only: AOM-5-5f.g2172.t2 is best explained by RFX3 (E=3.84e-122); near-miss E=3.84e-122 cov=0.637; genome hit (miniprot)
- RFX2 / Odio_OdB3: cross-hit only: CBY15668.1 is best explained by RFX3 (E=2.69e-123); near-miss E=2.69e-123 cov=0.640; genome hit (miniprot)
- RFX2 / Crob: partial: E=3.31e-23 profile cov=0.241; cross-hit only: tr|F6YLZ5|F6YLZ5_CIOIN is best explained by RFX3 (E=4e-212)
- RFX3 / Odio_KUM: partial: E=2.98e-89 profile cov=0.450
- MYB / Odio_OdB3: partial: E=3.65e-29 profile cov=0.209
- MYB / Crob: partial: E=1.71e-116 profile cov=0.464
- TP73 / Odio_OKI: near-miss E=0.000561 cov=0.150
- TP73 / Odio_OKI_refseq: near-miss E=0.000579 cov=0.150
- TP73 / Odio_KUM: near-miss E=0.000528 cov=0.150

## G_IFT_BBSome_motors: per-gene detail

● HMM+RBH  ◐ HMM only  ◑ RBH only  ◔ HMM-partial (strong E, low profile coverage)  ○ no proteome hit (genome search pending); G = miniprot genome hit at thresholds

| gene | Odio_OKI | Odio_OKI_refseq | Odio_OSA | Odio_BCN | Odio_KUM | Odio_AOM | Odio_OdB3 | Crob |
|---|---|---|---|---|---|---|---|---|
| IFT20 | ● 2.05e-24 | ● 1.95e-24 | ● 3.08e-23 | ◐ 2.68e-22 | ● 5.32e-24 | ● 6.02e-21 | ● 2.82e-24 | ● 1.18e-57 |
| IFT22 | ● 5.9e-26 | ● 8.44e-27 | ◐ 7.1e-26 | ◐ 1.09e-25 | ● 7.96e-27 | ◐ 1.38e-25 | ● 4.22e-27 | ● 2.75e-70 |
| IFT25 | ○ | ○ | ○ | ○ | ○ | ○ | ○ | ○ |
| IFT27 | ● 5.86e-50 | ● 1.95e-50 | ◐ 2.44e-37 | ◐ 2.21e-37 | ● 5.26e-50 | ◐ 2.36e-37 | ◐ 2.66e-37 | ● 2.28e-64 |
| IFT46 | ● 1.07e-70 | ● 1.02e-70 | ● 3.94e-73 | ● 1.99e-71 | ● 1.1e-70 | ● 3.82e-73 | ● 2.39e-71 | ● 2.07e-151 |
| IFT52 | ● 1.07e-116 | ● 1.02e-116 | ● 6.45e-116 | ◐ 5.41e-110 | ● 3.51e-84 | ◐ 2.37e-109 | ● 1.75e-115 | ● 3.84e-221 |
| IFT54 | ◐ 1.63e-99 | ● 4.4e-100 | ● 3.34e-94 | ● 8.32e-97 | ● 7.33e-101 | ● 7.15e-93 | ● 1.07e-96 | ● 3.25e-200 |
| IFT56 | ● 2.46e-225 | ● 2.35e-225 | ● 1.13e-221 | ○ (2.2e-10, cov 0.403) | ● 8.06e-226 | ● 1.5e-221 | ● 2.51e-224 | ● 9.48e-318 |
| IFT57 | ● 1.3e-65 | ● 1.31e-94 | ◐ 7.32e-58 | ● 9.01e-93 | ● 6.66e-66 | ◐ 4.19e-57 | ◔ E=1.67e-46 profile cov=0.420 | ● 1.35e-188 |
| IFT70A | ◐ 7.08e-245 | ◐ 5.16e-260 | ◐ 3.46e-260 | ◐ 5e-261 | ◐ 9.82e-257 | ◐ 4.94e-260 | ◐ 7.2e-261 | ◐ 7.25e-305 |
| IFT74 | ● 1.25e-105 | ● 1.19e-105 | ● 2.35e-89 | ● 1.48e-101 | ● 1.83e-105 | ● 2.65e-103 | ● 1.76e-100 | ● 1.94e-237 |
| IFT80 | ● 5.11e-119 | ● 4.87e-119 | ● 5.9e-161 | ● 1.22e-161 | ● 2.44e-151 | ● 1.68e-161 | ● 2.18e-159 | ● 0 |
| IFT81 | ● 2.08e-117 | ● 1.94e-91 | ● 1.02e-51 | ● 3.07e-116 | ● 2.63e-117 | ● 4.01e-117 | ● 6.9e-70 | ● 5.34e-294 |
| IFT88 | ● 1.84e-258 | ● 1.81e-253 | ● 1.74e-248 | ● 1.59e-258 | ● 5.02e-259 | ● 5.53e-259 | ● 1.98e-258 | ● 0 |
| IFT172 | ◐ 1.63e-272 | ● 0 | ● 0 | ● 0 | ● 0 | ● 0 | ◐ 1.31e-298 | ◔ E=6.61e-275 profile cov=0.291 |
| CLUAP1 | ● 1.5e-80 | ● 1.43e-80 | ● 8.93e-70 | ● 8.93e-74 | ● 1.95e-80 | ● 3.8e-75 | ● 6.77e-69 | ● 6.63e-148 |
| IFT43 | ○ (3.47e-06, cov 0.282) | ○ (3.31e-06, cov 0.282) | ○ (1.63e-05, cov 0.248) | ○ | ○ (4.57e-06, cov 0.262) | ○ (2.91e-05, cov 0.248) | ○ (0.000101, cov 0.257) | ● 3.07e-67 |
| IFT121 | ● 3.56e-317 | ● 3.39e-317 | ● 2.13e-309 | ● 4.18e-317 | ● 9.08e-318 | ● 1.22e-317 | ● 2.14e-316 | ◔ E=3.88e-188 profile cov=0.342 |
| IFT122 | ● 5.9e-292 | ● 2.03e-291 | ◐ 4.12e-168 | ● 8.05e-297 | ● 1.75e-296 | ● 1.27e-297 | ◐ 9.42e-169 | ● 0 |
| IFT139 | ● 7.61e-74 | ● 1.52e-137 | ◐ 7.32e-81 | ● 1.31e-143 | ● 1.33e-139 | ◔ E=5.61e-61 profile cov=0.274 | ● 4.56e-137 | ● 0 |
| IFT140 | ● 1.89e-262 | ● 9.74e-314 | ● 3.31e-303 | ● 1.07e-237 | ● 3.53e-315 | ● 2.14e-286 | ● 4.66e-312 | ● 0 |
| IFT144 | ● 6.18e-229 | ◔ E=1.6e-121 profile cov=0.360 | ● 0 | ● 0 | ● 0 | ● 0 | ● 0 | ● 0 |
| BBS1 | ◐ 1.17e-108 | ◔ E=2.28e-53 profile cov=0.394 | ◔ E=4.68e-48 profile cov=0.406 | ● 7.64e-136 | ◐ 1.1e-108 | ● 1.42e-135 | ● 3.49e-102 | ● 1.85e-262 |
| BBS2 | ● 9.16e-129 | ● 8.72e-129 | ◐ 1.21e-87 | ◐ 5.45e-119 | ● 5.12e-129 | ● 7.98e-127 | ● 4.92e-127 | ● 0 |
| BBS4 | ◐ 2.64e-38 | ◐ 2.89e-38 | ● 8.1e-53 | ● 7.22e-43 | ◐ 3.49e-53 | ● 7.85e-53 | ● 7.33e-43 | ● 4.28e-214 |
| BBS5 | ● 1.11e-99 | ● 1.06e-99 | ● 7.2e-100 | ● 2.8e-98 | ● 9.99e-100 | ● 1.11e-99 | ● 2.22e-98 | ● 1.49e-171 |
| BBS7 | ● 8.17e-41 | ◔ E=1.2e-17 profile cov=0.237 | ● 2.1e-48 | ● 1.21e-48 | ● 6.15e-47 | ● 5.61e-41 | ● 2.48e-48 | ● 0 |
| BBS8 | ◐ 8.29e-74 | ● 9.58e-89 | ● 8.08e-93 | ● 2.14e-94 | ◐ 6.68e-91 | ● 2.37e-89 | ○ (5.6e-32, cov 0.482) | ● 1.88e-268 |
| BBS9 | ○ (1.85e-08, cov 0.135) | ○ (6.13e-06, cov 0.099) | ◔ E=2.31e-17 profile cov=0.480 | ◔ E=6.81e-18 profile cov=0.463 | ○ (4.76e-08, cov 0.135) | ◔ E=1.51e-17 profile cov=0.487 | ◔ E=4.65e-13 profile cov=0.327 | ◐ 1.55e-149 |
| BBIP1 | ○ | ○ | ○ | ○ | ○ | ○ | ○ | ○ |
| ARL6 | ● 3.85e-68 | ● 3.67e-68 | ● 2.25e-68 | ● 2.18e-68 | ● 2.51e-68 | ● 1.75e-68 | ● 2.62e-68 | ● 4.26e-95 |
| LZTFL1 | ● 4.92e-47 | ● 1.75e-55 | ● 1.38e-54 | ● 9.41e-54 | ● 8.43e-56 | ● 7.82e-55 | ● 4.96e-53 | ● 3.75e-131 |
| BBS10 | ◐ 3.34e-25 | ◐ 4.91e-25 | ◐ 1.47e-23 | ◐ 4.24e-24 | ◐ 4.23e-25 | ◐ 1.13e-23 | ◐ 4.1e-24 | ◐ 5.05e-29 |
| BBS12 | ○ (7.19e-11, cov 0.345) | ○ (6.85e-11, cov 0.345) | ○ (1.18e-10, cov 0.352) | ○ (2.15e-10, cov 0.421) | ○ (6.46e-11, cov 0.345) | ○ (9.06e-11, cov 0.352) | ○ (1.81e-10, cov 0.428) | ○ (1.01e-11, cov 0.582) |
| MKKS | ◐ 4.79e-21 | ◐ 4.57e-21 | ◐ 5.96e-20 | ◐ 1.33e-19 | ◐ 4.31e-21 | ◐ 5.78e-20 | ◐ 1.3e-19 | ◔ E=3.46e-14 profile cov=0.202 |
| DYNC2H1 | ● 0 | ● 0 | ● 0 | ● 0 | ◐ 0 | ● 0 | ● 0 | ● 0 |
| DYNC2I1 | ◔ E=1.16e-17 profile cov=0.384 | ◔ E=1.1e-17 profile cov=0.384 | ◔ E=9.74e-17 profile cov=0.380 | ◔ E=4.13e-19 profile cov=0.382 | ◔ E=5.51e-17 profile cov=0.331 | ◔ E=1.54e-17 profile cov=0.376 | ◔ E=1.01e-17 profile cov=0.341 | ◔ E=1.46e-60 profile cov=0.402 |
| DYNC2I2 | ◐ 1.24e-14 | ◐ 1.18e-14 | ◐ 5.98e-15 | ◔ E=2.32e-13 profile cov=0.253 | ◐ 8.46e-15 | ◐ 1.67e-14 | ◐ 1.99e-28 | ● 3.61e-182 |
| DYNC2LI1 | ◐ 5.17e-18 | ◐ 4.93e-18 | ◐ 5.12e-18 | ◐ 3.51e-18 | ● 2.38e-22 | ◐ 4.24e-18 | ◐ 2.21e-18 | ● 1.04e-144 |
| DYNLT2B | ◐ 4.29e-17 | ◐ 9.04e-18 | ● 3.03e-33 | ● 1.28e-32 | ◐ 5.44e-18 | ● 2.93e-33 | ● 1.37e-33 | ● 2.2e-62 |
| KIF3A | ● 3.47e-273 | ● 1.91e-273 | ◐ 5.28e-166 | ● 7.51e-272 | ● 9.15e-228 | ● 3.75e-271 | ● 4.19e-264 | ● 0 |
| KIF3B | ◐ 4.26e-278 | ● 1.77e-294 | ● 1.16e-259 | ● 1.85e-254 | ● 4.28e-289 | ● 2.59e-289 | ● 8.17e-290 | ● 0 |
| KIF3C | ○ G | ○ (5.03e-239, cov 0.825) | ○ G | ○ G | ○ G | ○ G | ○ G | ○ (1.66e-289, cov 0.900) |
| KIFAP3 | ● 1.67e-124 | ● 4.75e-221 | ◔ E=3.88e-51 profile cov=0.261 | ● 4.11e-117 | ● 3.84e-105 | ● 2.18e-198 | ● 3.1e-224 | ● 0 |
| KIF17 | ● 2.15e-266 | ● 2.05e-266 | ● 8.52e-230 | ◐ 2.3e-253 | ● 1.1e-266 | ● 1.78e-264 | ● 1.87e-264 | ◐ 5.34e-237 |
| TULP3 | ◐ 1.82e-120 | ◐ 1.96e-120 | ◐ 4.41e-120 | ◐ 3.24e-122 | ◐ 3.04e-120 | ◐ 9.22e-120 | ◐ 3.6e-121 | ◐ 1.35e-150 |
| ARL3 | ● 2.62e-88 | ● 2.5e-88 | ● 4.95e-89 | ● 2.49e-88 | ● 9.08e-88 | ● 4.79e-89 | ● 7.5e-89 | ● 7.15e-102 |
| UNC119 | ◔ E=2.1e-11 profile cov=0.262 | ● 9.2e-57 | ◐ 1.97e-52 | ◔ E=1.38e-11 profile cov=0.266 | ◔ E=2.04e-11 profile cov=0.262 | ◐ 3.16e-52 | ◔ E=2.12e-11 profile cov=0.266 | ◐ 4.79e-96 |

Notes:
- IFT56 / Odio_BCN: near-miss E=2.2e-10 cov=0.403
- IFT57 / Odio_OdB3: partial: E=1.67e-46 profile cov=0.420
- IFT81 / Odio_OKI: HMM best CAG5095129.1 != RBH CAG5095131.1
- IFT172 / Crob: partial: E=6.61e-275 profile cov=0.291
- IFT43 / Odio_OKI: near-miss E=3.47e-06 cov=0.282
- IFT43 / Odio_OKI_refseq: near-miss E=3.31e-06 cov=0.282
- IFT43 / Odio_OSA: near-miss E=1.63e-05 cov=0.248
- IFT43 / Odio_KUM: near-miss E=4.57e-06 cov=0.262
- IFT43 / Odio_AOM: near-miss E=2.91e-05 cov=0.248
- IFT43 / Odio_OdB3: near-miss E=0.000101 cov=0.257
- IFT121 / Crob: partial: E=3.88e-188 profile cov=0.342
- IFT139 / Odio_AOM: partial: E=5.61e-61 profile cov=0.274
- IFT144 / Odio_OKI_refseq: partial: E=1.6e-121 profile cov=0.360; cross-hit only: XP_081631388.1 is best explained by IFT140 (E=7.57e-25)
- BBS1 / Odio_OKI_refseq: partial: E=2.28e-53 profile cov=0.394
- BBS1 / Odio_OSA: partial: E=4.68e-48 profile cov=0.406
- BBS7 / Odio_OKI_refseq: partial: E=1.2e-17 profile cov=0.237
- BBS8 / Odio_OdB3: cross-hit only: CBY09563.1 is best explained by BBS4 (E=2.74e-20); near-miss E=5.6e-32 cov=0.482
- BBS9 / Odio_OKI: near-miss E=1.85e-08 cov=0.135
- BBS9 / Odio_OKI_refseq: near-miss E=6.13e-06 cov=0.099
- BBS9 / Odio_OSA: partial: E=2.31e-17 profile cov=0.480
- BBS9 / Odio_BCN: partial: E=6.81e-18 profile cov=0.463
- BBS9 / Odio_KUM: near-miss E=4.76e-08 cov=0.135
- BBS9 / Odio_AOM: partial: E=1.51e-17 profile cov=0.487
- BBS9 / Odio_OdB3: partial: E=4.65e-13 profile cov=0.327
- BBS12 / Odio_OKI: cross-hit only: CAG5090001.1 is best explained by BBS10 (E=1.14e-09); near-miss E=7.19e-11 cov=0.345
- BBS12 / Odio_OKI_refseq: cross-hit only: XP_081635890.1 is best explained by BBS10 (E=1.09e-09); near-miss E=6.85e-11 cov=0.345
- BBS12 / Odio_OSA: cross-hit only: OSKA2016v1.9.g11480.t1 is best explained by BBS10 (E=1e-09); near-miss E=1.18e-10 cov=0.352
- BBS12 / Odio_BCN: cross-hit only: Bar2_p4_Flye.g1738.t1 is best explained by BBS10 (E=9.63e-07); near-miss E=2.15e-10 cov=0.421
- BBS12 / Odio_KUM: cross-hit only: KUM-M3-7f.g5956.t1 is best explained by BBS10 (E=1.11e-09); near-miss E=6.46e-11 cov=0.345
- BBS12 / Odio_AOM: cross-hit only: AOM-5-5f.g7280.t1 is best explained by BBS10 (E=5.57e-10); near-miss E=9.06e-11 cov=0.352
- BBS12 / Odio_OdB3: near-miss E=1.81e-10 cov=0.428
- BBS12 / Crob: cross-hit only: tr|F6W7K5|F6W7K5_CIOIN is best explained by BBS10 (E=1.01e-11); near-miss E=1.01e-11 cov=0.582
- MKKS / Crob: partial: E=3.46e-14 profile cov=0.202; cross-hit only: tr|F6W7K5|F6W7K5_CIOIN is best explained by BBS10 (E=1.81e-23)
- DYNC2I1 / Odio_OKI: partial: E=1.16e-17 profile cov=0.384; cross-hit only: CAG5112518.1 is best explained by DNAI2 (E=1.16e-14)
- DYNC2I1 / Odio_OKI_refseq: partial: E=1.1e-17 profile cov=0.384; cross-hit only: XP_081634063.1 is best explained by DNAI2 (E=1.1e-14)
- DYNC2I1 / Odio_OSA: partial: E=9.74e-17 profile cov=0.380; cross-hit only: OSKA2016v1.9.g4620.t1 is best explained by DNAI2 (E=2.32e-13)
- DYNC2I1 / Odio_BCN: partial: E=4.13e-19 profile cov=0.382
- DYNC2I1 / Odio_KUM: partial: E=5.51e-17 profile cov=0.331
- DYNC2I1 / Odio_AOM: partial: E=1.54e-17 profile cov=0.376; cross-hit only: AOM-5-5f.g5129.t1 is best explained by DNAI2 (E=7.05e-16)
- DYNC2I1 / Odio_OdB3: partial: E=1.01e-17 profile cov=0.341
- DYNC2I1 / Crob: partial: E=1.46e-60 profile cov=0.402
- DYNC2I2 / Odio_BCN: partial: E=2.32e-13 profile cov=0.253; cross-hit only: Bar2_p4_Flye.g7609.t1 is best explained by DNAI1 (E=1.95e-50)
- KIF3C / Odio_OKI: cross-hit only: CAG5086379.1 is best explained by KIF3B (E=3.21e-233); near-miss E=3.21e-233 cov=0.888; genome hit (miniprot)
- KIF3C / Odio_OKI_refseq: cross-hit only: XP_081637647.1 is best explained by KIF3B (E=5.03e-239); near-miss E=5.03e-239 cov=0.825
- KIF3C / Odio_OSA: cross-hit only: OSKA2016v1.9.g10532.t1 is best explained by KIF3B (E=4.01e-212); near-miss E=4.01e-212 cov=0.788; genome hit (miniprot)
- KIF3C / Odio_BCN: cross-hit only: Bar2_p4_Flye.g9350.t1 is best explained by KIF3B (E=7.54e-210); near-miss E=7.54e-210 cov=0.825; genome hit (miniprot)
- KIF3C / Odio_KUM: cross-hit only: KUM-M3-7f.g15696.t1 is best explained by KIF3B (E=3.09e-239); near-miss E=3.09e-239 cov=0.818; genome hit (miniprot)
- KIF3C / Odio_AOM: cross-hit only: AOM-5-5f.g3126.t1 is best explained by KIF3B (E=2.44e-239); near-miss E=2.44e-239 cov=0.820; genome hit (miniprot)
- KIF3C / Odio_OdB3: cross-hit only: CBY09744.1 is best explained by KIF3B (E=1.6e-241); near-miss E=1.6e-241 cov=0.821; genome hit (miniprot)
- KIF3C / Crob: cross-hit only: tr|F7B875|F7B875_CIOIN is best explained by KIF3B (E=1.66e-289); near-miss E=1.66e-289 cov=0.900
- KIFAP3 / Odio_OSA: partial: E=3.88e-51 profile cov=0.261
- KIF17 / Odio_OKI: HMM best CAG5111451.1 != RBH CAG5111450.1
- UNC119 / Odio_OKI: partial: E=2.1e-11 profile cov=0.262
- UNC119 / Odio_BCN: partial: E=1.38e-11 profile cov=0.266
- UNC119 / Odio_KUM: partial: E=2.04e-11 profile cov=0.262
- UNC119 / Odio_OdB3: partial: E=2.12e-11 profile cov=0.266

## Genes with no proteome hit at thresholds in any of these 8 species (54; HMM-partial allowed)

Candidates for the genome search; listed by module.

- B_duplication_cartwheel: STIL, CEP152, CEP192, CEP57L1, CEP295, CEP44, CNTROB, SPICE1, WDR62
- C_distal_end_length: KIF24, TTBK2, KIAA0586
- D_PCM_gTuRC: TUBGCP6, MZT2B, PCNT, CDK5RAP2, CEP170
- E_appendages_basal_foot: NIN, CEP128, CCDC120, CCDC68, CEP89, KIAA0753
- F_transition_zone: TCTN2, CEP162, NPHP3
- G_IFT_BBSome_motors: IFT25, BBIP1, BBS12, DYNC2I1, KIF3C
- H_motile_axoneme: CFAP54, DNAAF2
- I_multiciliogenesis: MCIDAS, GMNC, DEUP1
- J_ciliary_signaling: DHH, PTCH2, GLI1, GLI3, KIF7, ARL13B, PKD1L1, EVC, EVC2
- K_linkers_satellites_IDPs: CEP250, CEP68, PCM1, ALMS1
- L_tubulin_PTM: TTLL4, TTLL6, TTLL10, AGBL3, AGBL1
