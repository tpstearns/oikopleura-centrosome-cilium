# Gene anchor resolution (genes_resolve), 2026-10-20

272 of 272 genes in `config/gene_modules.tsv` resolved to a UniProt anchor sequence; 0 unresolved.
Resolution: 271 by exact gene symbol against reviewed human (Swiss-Prot) entries; 1 by explicit accession (TUBZ1 → *Xenopus tropicalis* tubz1, A0A6I8RNG6, unreviewed; see docs/decisions.md).

Table: `results/anchors.tsv` (module, gene, anchor species, taxid, accession, entry, UniProt primary gene, protein name, length, how resolved).

## Genes per module

| module | n |
|---|---|
| A_tubulins | 7 |
| B_duplication_cartwheel | 24 |
| C_distal_end_length | 11 |
| D_PCM_gTuRC | 18 |
| E_appendages_basal_foot | 19 |
| F_transition_zone | 24 |
| G_IFT_BBSome_motors | 48 |
| H_motile_axoneme | 51 |
| I_multiciliogenesis | 14 |
| J_ciliary_signaling | 21 |
| K_linkers_satellites_IDPs | 8 |
| L_tubulin_PTM | 23 |
| M_controls | 4 |

## Symbols whose UniProt primary gene name differs (official renames; the resolved protein is the intended one)

| list symbol | UniProt primary | accession | protein |
|---|---|---|---|
| CENPJ | CPAP | Q9HC77 | Centrosomal P4.1-associated protein |
| FGFR1OP | CEP43 | O95684 | Centrosomal protein 43 |
| FOPNL | CEP20 | Q96NB1 | Centrosomal protein 20 |
| CLUAP1 | IFT38 | Q96AJ1 | Intraflagellar transport protein 38 |
| IFT121 | WDR35 | Q9P2L0 | WD repeat-containing protein 35 |
| IFT139 | TTC21B | Q7Z4L5 | Tetratricopeptide repeat protein 21B |
| IFT144 | WDR19 | Q8NEZ3 | WD repeat-containing protein 19 |
| BBS8 | TTC8 | Q8TAM2 | Tetratricopeptide repeat protein 8 |
| CCDC65 | DRC2 | Q8IXS2 | Dynein regulatory complex subunit 2 |
| GAS8 | DRC4 | O95995 | Dynein regulatory complex subunit 4 |
| LRRC6 | DNAAF11 | Q86X45 | Dynein axonemal assembly factor 11 |
