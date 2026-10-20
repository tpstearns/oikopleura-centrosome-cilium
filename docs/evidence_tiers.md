# Evidence tiers

Every cell of the gene × species matrix carries one tier. Thresholds are in `config/params.yaml`; the logic is in `scripts/tiers.py`. A tier is a statement about evidence, not about biology: E4 means "not found by these searches at these thresholds in this assembly," and it becomes a loss claim only after the sensitivity controls for that species pass and the genome search has been run.

| Tier | Meaning | Requirements |
|---|---|---|
| E1 | Confident ortholog with intact architecture | Reciprocal best hit to the anchor (E ≤ 1e-5, ≥50% anchor coverage) **and** per-gene profile hit (E ≤ 1e-5, best domain i-E ≤ 1e-3, ≥50% profile coverage) **and**, when Pfam is enabled, ≥50% of the anchor's Pfam families present |
| E2 | Homolog present, identity or architecture not fully established | Profile hit, RBH, or Foldseek match passing thresholds, but not all E1 conditions (paralog ambiguity, partial architecture, short divergent protein) |
| E3 | Present in genome or transcripts, missing from annotation | No proteome hit; miniprot spliced alignment (identity ≥ 25%, ≥50% query coverage) or transcript hit passing thresholds |
| E4 | No evidence | Nothing above threshold in proteome, genome, or transcripts |

## Paralog families
A profile hit counts toward a gene only when that gene's profile is the best-scoring profile for the target protein. Tubulin, cyclin, E2F and kinesin profiles all hit their whole family; without this rule a β-tubulin would register as ε-tubulin wherever the true ε protein is missing from the annotation. Cross-hits are recorded in the notes column. Reciprocal best hits are reciprocal against the full anchor-species proteome, so a paralog outside the gene list (CDC20 for CDC20B) cannot be returned as the listed gene.

## Why absence is hard here
- *Oikopleura* proteins are highly divergent, so sequence similarity to human orthologs is often marginal. Per-gene profiles built from diverse seeds (vertebrates, amphioxus, *Ciona*, *Chlamydomonas*, *Tetrahymena*, fly, worm) are used instead of pairwise search alone.
- Annotations are incomplete and differ among the cryptic-lineage assemblies; the genome-level search (miniprot) is mandatory before any E4 is reported for a core species.
- Short proteins and coiled-coil-rich proteins (many centriolar proteins) give weak profile statistics; profile coverage thresholds are therefore moderate and E2 is reserved for these cases.

## Sensitivity controls
For each core species, a set of genes known to be present (ubiquitous controls in module M plus genes previously documented in *O. dioica*) is run blind through the same pipeline. The fraction recovered at E1/E2 is reported alongside the matrix; absence calls are not reported for a species whose control recovery is below the threshold set in `config/params.yaml` (to be added once controls are curated).

## Lineage-specific differences
A gene at E1/E2 in one *O. dioica* lineage and E4 in another is flagged for manual review (assembly gap vs. true loss), using the genome search and, where possible, the haplotype-resolved Okinawa assembly as the arbiter.
