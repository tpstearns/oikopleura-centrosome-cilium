# Species panel after `panel_manifest` + `fetch_all` (2026-10-20)

39 entries; 33 proteomes and 12 genomes downloaded (checksums in `data/manifest/files.tsv`). Full table: `results/panel_summary.tsv`. Proteome counts are sequences in the file: UniProt reference proteomes are one canonical sequence per gene; NCBI RefSeq/GenBank protein sets include all isoforms; Augustus sets (Plessy 2024 Zenodo supplement) are longest-isoform ab initio models and are uncurated.

| key | species / lineage | role | proteome (source) | n proteins | genome (source) | seqs / bp | status |
|---|---|---|---|---|---|---|---|
| Odio_OKI | Oikopleura dioica — Okinawa, OKI2018_I69_1.0 | core | GCA_907165135.1 (ncbi_datasets) | 18751 | GCA_907165135.1 (ncbi_datasets) | 19 / 64,281,565 | verified |
| Odio_OKI_refseq | Oikopleura dioica — Okinawa, OKI2018 — same assembly as Odio_OKI, NCBI RefSeq annotation RS_2026_03 | core | GCF_907165135.1 (ncbi_datasets) | 17869 | — (—) | — | verified |
| Odio_OKI_hap | Oikopleura dioica — Okinawa, haplotype-resolved assembly | core | — (—) | — | PRJEB91365 (ena_bioproject) | — | to-verify |
| Odio_OSA | Oikopleura dioica — Osaka, OSKA2016v1.9 | core | OSKA2016v1.9.prot.longest.fa (zenodo_archive) | 15720 | OSKA2016v1.9.fa (zenodo_archive) | 483 / 56,625,162 | verified |
| Odio_BCN | Oikopleura dioica — Barcelona, Bar2_p4 | core | Bar2_p4.Flye.prot.longest.fa (zenodo_archive) | 14272 | Bar2_p4.fa (zenodo_archive) | 68 / 55,793,437 | verified |
| Odio_KUM | Oikopleura dioica — Kume island, KUM-M3-7f | core | KUM-M3-7f.prot.longest.fa (zenodo_archive) | 16852 | KUM-M3-7f.fa (zenodo_archive) | 55 / 64,653,574 | verified |
| Odio_AOM | Oikopleura dioica — Aomori, AOM-5-5f | core | AOM-5-5f.prot.longest.fa (zenodo_archive) | 15224 | AOM-5-5f.fa (zenodo_archive) | 33 / 56,753,784 | verified |
| Odio_OdB3 | Oikopleura dioica — Bergen | core | GCA_000209535.1 (ncbi_datasets) | 17152 | GCA_000209535.1 (ncbi_datasets) | 1,260 / 70,471,451 | verified |
| Oalb | Oikopleura albicans | larvacean | O_alb.prot.longest.fa (zenodo_archive) | 24943 | GCA_004367875.1 (ncbi_datasets) | 92,906 / 365,958,051 | verified |
| Olon | Oikopleura longicauda | larvacean | O_lon.prot.longest.fa (zenodo_archive) | 27566 | GCA_004367895.1 (ncbi_datasets) | 178,440 / 308,740,099 | verified |
| Ovan | Oikopleura vanhoeffeni | larvacean | O_van.prot.longest.fa (zenodo_archive) | 18686 | GCA_004367855.1 (ncbi_datasets) | 79,827 / 643,662,591 | verified |
| Bsty | Bathochordaeus stygius | larvacean | B_sty.prot.longest.fa (zenodo_archive) | 25023 | GCA_004367955.1 (ncbi_datasets) | 467,514 / 396,466,979 | verified |
| Mery | Mesochordaeus erythrocephalus | larvacean | M_ery.prot.longest.fa (zenodo_archive) | 48793 | GCA_004367975.1 (ncbi_datasets) | 745,784 / 874,018,965 | verified |
| Fbor | Fritillaria borealis | larvacean | F_bor.prot.longest.fa (zenodo_archive) | 17861 | GCA_004368075.1 (ncbi_datasets) | 142,328 / 143,085,447 | verified |
| Crob | Ciona intestinalis | ascidian | UP000008144 (uniprot_reference) | 16680 | GCF_000224145.3 (ncbi_datasets) | — | verified |
| Cint | Ciona intestinalis | ascidian | GCF_018327825.1 (ncbi_datasets) | 19614 | GCF_018327825.1 (ncbi_datasets) | — | verified |
| Csav | Ciona savignyi | ascidian | UP000007875 (uniprot_reference) | 11592 | GCA_054131345.1 (ncbi_datasets) | — | verified |
| Pmam | Phallusia mammillata | ascidian | GCF_965637545.1 (ncbi_datasets) | 27517 | GCF_965637545.1 (ncbi_datasets) | — | verified |
| Hror | Halocynthia roretzi | ascidian | — (to-resolve) | — | GCA_013436055.1 (ncbi_datasets) | — | to-verify |
| Bsch | Botryllus schlosseri | ascidian | GCF_051294905.1 (ncbi_datasets) | 29009 | GCF_051294905.1 (ncbi_datasets) | — | verified |
| Blea | Botrylloides leachii | ascidian | — (to-resolve) | — | — (—) | — | to-verify |
| Mocc | Molgula occulta | ascidian | — (to-resolve) | — | — (—) | — | to-verify |
| Mocu | Molgula oculata | ascidian | — (to-resolve) | — | — (—) | — | to-verify |
| Scla | Styela clava | ascidian | GCF_964204865.1 (ncbi_datasets) | 26946 | GCF_964204865.1 (ncbi_datasets) | — | verified |
| Clep | Clavelina lepadiformis | ascidian | UP001642483 (uniprot_reference) | 22422 | GCF_947623445.1 (ncbi_datasets) | — | verified |
| Stho | Salpa thompsoni | thaliacean | — (to-resolve) | — | GCA_047495815.1 (ncbi_datasets) | — | to-verify |
| Tdem | Thalia democratica | thaliacean | GCF_965202585.1 (ncbi_datasets) | 19084 | GCF_965202585.1 (ncbi_datasets) | — | verified |
| Bflo | Branchiostoma floridae | backbone | UP000001554 (uniprot_reference) | 26616 | — (—) | — | verified |
| Blan | Branchiostoma lanceolatum | backbone | UP000838412 (uniprot_reference) | 26322 | — (—) | — | verified |
| Drer | Danio rerio | backbone | UP000000437 (uniprot_reference) | 26567 | — (—) | — | verified |
| Xtro | Xenopus tropicalis | backbone | UP000008143 (uniprot_reference) | 21732 | — (—) | — | verified |
| Hsap | Homo sapiens | backbone | UP000005640 (uniprot_reference) | 20652 | — (—) | — | verified |
| Mmus | Mus musculus | backbone | UP000000589 (uniprot_reference) | 21860 | — (—) | — | verified |
| Crei | Chlamydomonas reinhardtii | control | UP000006906 (uniprot_reference) | 17625 | — (—) | — | verified |
| Tthe | Tetrahymena thermophila | control | UP000009168 (uniprot_reference) | 26975 | — (—) | — | verified |
| Dmel | Drosophila melanogaster | control | UP000000803 (uniprot_reference) | 13814 | — (—) | — | verified |
| Cele | Caenorhabditis elegans | control | UP000001940 (uniprot_reference) | 19792 | — (—) | — | verified |
| Atha | Arabidopsis thaliana | control | UP000006548 (uniprot_reference) | 27496 | — (—) | — | verified |
| Scer | Saccharomyces cerevisiae | control | UP000002311 (uniprot_reference) | 6066 | — (—) | — | verified |

## Not fetched

- Odio_OKI_hap: PRJEB91365 has reads only (no assembly record).
- Hror, Blea, Mocc, Mocu: ANISEED-only proteomes; ANISEED download links redirect for scripted requests.
- Stho: no annotated assembly (genome GCA_047495815.1 is contig-level, unannotated; not in genome_search_species).

## Also in the Plessy 2024 Augustus archive (not used yet)

- A second OIST annotation of OKI2018_I69 (`OKI2018_I69_annotv2`), Augustus runs for both *Ciona intestinalis* types, and `GSE64421_transcripts.trinity.fasta` (a Trinity transcriptome used as hints; a candidate transcript source for the OdB3/Bergen lineage once its provenance is confirmed).
