# Gene List Matching — CellFM Vocabulary Compatibility

## Context

As part of replicating HINGE (Fang et al., CVPR 2026) — which adapts the pre-trained
single-cell foundation model **CellFM** to spatial gene expression generation from
histology — the ST dataset's gene list must be intersected with CellFM's fixed gene
vocabulary before HMHVG (Highly Mean–High Variance Gene) selection can be performed.

This intersection step was not matching cleanly, so the gene lists were audited and
partially corrected. This README documents the files in this folder and the steps
taken so far.

## Files

| File | Description |
|---|---|
| `used_selected_gene_list_original.txt` | Original gene list (18,085 genes) before any correction. |
| `used_selected_gene_list.txt` | Working list after a first attempt at alias correction (17,662 genes). |
| `genemap.csv` | Gene → index lookup table intended to represent the CellFM gene vocabulary (27,934 unique genes, columns: gene symbol, `idx`). |
| `added_genes.txt` | 444 genes present in `used_selected_gene_list.txt` but absent from the original list — i.e., genes introduced during the correction pass (mostly renamed/alias symbols). |
| `unmatched_genes.txt` | 159 genes from `used_selected_gene_list.txt` that still do **not** match any entry in `genemap.csv`. |

## What was found

1. **Symbol mismatches (mostly resolved).** A large share of genes in the original
   list used older/alternate HGNC symbols not present in `genemap.csv`, while
   `genemap.csv` uses updated symbols. Examples:

   | Old symbol (original list) | New symbol (in genemap.csv) |
   |---|---|
   | AARS | AARS1 |
   | ARNTL | BMAL1 |
   | ADSS | ADSS1, ADSS2 |
   | ATP5MD | ATP5MK |
   | ARSE | ARSL |

   These were corrected by replacing the old symbol with the genemap-matching one,
   accounting for most of the 444 genes in `added_genes.txt`.

2. **546 genes were unnecessarily dropped.** During the first correction pass, 546
   genes that were already present as-is in `genemap.csv` (e.g., ACAT1, ACAT2, ACP1,
   ADA2, ADCY3, AGT, AIP, AK6, ALB, APC2, APPL1, AR, ARC) were removed from the list
   for reasons unrelated to symbol matching. These have since been identified and can
   be added back.

3. **159 genes remain unmatched (`unmatched_genes.txt`).** These do not correspond to
   any entry in `genemap.csv` under any symbol variant checked so far. Almost all of
   them are Ensembl clone/accession-style identifiers rather than standard HGNC gene
   symbols (e.g., `AC004593.3`, `AL109810.2`, `AC010325.1`), plus a small number of
   symbols such as `ADGRF2P` with no corresponding genemap entry. These are likely
   novel transcripts / non-standard identifiers that were never part of the CellFM
   training vocabulary, but this has not been independently confirmed against
   CellFM's original vocabulary file.

## What has been tried

- Diffed the original and working gene lists to identify all additions/removals
  (867 removed, 444 added between the two lists).
- Cross-checked every gene in both lists directly against `genemap.csv` by exact
  string match.
- Manually verified several rename cases (AARS/AARS1, ARNTL/BMAL1, ADSS/ADSS1+ADSS2,
  ATP5MD/ATP5MK, ARSE/ARSL) by confirming the new symbol's presence in `genemap.csv`.
- Checked for case-sensitivity issues (none found — no gene matches only after
  case-folding).
- Computed a corrected candidate list (not yet finalized/uploaded) that combines
  valid original-list genes with valid renamed genes, achieving a 100% match rate
  against `genemap.csv` (18,049 genes) — this still excludes the 159 unmatched genes,
  which have no resolvable counterpart in `genemap.csv`.

## Open question / discrepancy to flag

The HINGE paper states that CellFM was pre-trained on a **fixed 24,078-gene
vocabulary**, but `genemap.csv` contains **27,934 unique genes**. This discrepancy
hasn't been resolved yet — `genemap.csv` may be a different or newer version of the
vocabulary than the one used in the paper, which would need to be verified against
the original CellFM release before the 159 unmatched genes can be conclusively
labeled as "not in vocabulary" versus "not in this particular genemap file."

## Next steps

- Confirm the correct/original CellFM gene vocabulary file (24,078 genes) and re-run
  the matching against that, rather than `genemap.csv`, to resolve the count
  discrepancy.
- Decide whether to add back the 546 unnecessarily dropped genes.
- Decide how to handle the 159 unmatched genes (drop from HMHVG candidate pool, since
  they cannot be represented in CellFM's embedding space).
