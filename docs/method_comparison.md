# What the paper did vs what the guide specifies vs what this pipeline does

This table exists so that no sentence in the report claims to have reproduced a
method that was not actually run. Column C is the truth about this repository.

| Step | A. Bise et al. 2023 (Seurat, R) | B. Project guide (Scanpy) | C. This implementation | Equivalent? |
|---|---|---|---|---|
| Alignment | Cell Ranger 3.0.2, GRCz11 (NCBI) **manually modified to add the EGFP sequence** | not performed by the student | not performed; we start from the supplied filtered matrices | n/a — inherited |
| Cell QC | `scater`: low library size, low gene count, high mitochondrial fraction; ~12 ± 2% of cells removed per condition as doublets/dying | `min_genes=200`, `min_cells=3`, `pct_mt<15` | same as B, applied after per-sample diagnostics | **No.** Different criteria and no explicit doublet step |
| Doublet removal | explicit, via scater-based filtering | none | none | **No.** A stated limitation |
| Normalisation | **SCTransform** — regularised negative-binomial Pearson residuals, with mitochondrial and ribosomal percentage regressed out during normalisation | `normalize_total(1e4)` + `log1p` | same as B | **No.** Different estimator, not a reimplementation |
| Confounder removal | inside SCTransform (mito + ribo %) | `regress_out(["total_counts","pct_counts_mt"])` | same covariates as B, but applied after HVG subsetting for tractability (documented in `config.py`) | Partially |
| Feature selection | implicit in SCTransform variable-feature selection | 2,000 HVGs, `batch_key="sample"` | same as B | **No** |
| Integration | Seurat anchor-based (CCA), `dims = 30` | not specified | **not applied by default**; decided from PC-variance diagnostics in Step 04 and recorded either way | **No.** Deliberate, see below |
| Graph / clustering | KNN + SNN modularity optimisation, **resolution 0.2** | Leiden, **resolution 0.6** | Leiden 0.6 as baseline, with a 0.2–1.0 stability sweep reported | **No.** Different algorithm and resolution |
| Embedding | UMAP | UMAP | UMAP | Yes |
| Cluster markers | `FindAllMarkers` | `rank_genes_groups` | `rank_genes_groups`, `method="wilcoxon"`, on log-normalised values | Yes, same test |
| DE vs control | `FindMarkers`, Wilcoxon rank-sum, per cell type per timepoint | not specified | same test, per cell type per timepoint | Yes, same test |
| Enrichment | `topGO`, Fisher's exact test | `gseapy` available | `gseapy` (optional) | **No.** Different gene-set backend and background handling |
| EGFP⁺ cells | Seurat `subset` on EGFP transcripts | violin/UMAP of EGFP | raw counts ≥ 1, threshold justified from the observed count distribution | Comparable in spirit; the paper does not state its exact cut |
| Cluster count | 17, then bipolar and GABAergic subclusters merged | not specified | whatever Leiden 0.6 yields; annotated by marker, merged only with a written rationale | **No** |

## Why integration is off by default

The paper integrated across samples. This project does not, unless the Step 04
diagnostics justify it, for a design reason: each timepoint is its own pair of 10x
runs, so `sample` and `condition` are fully confounded. Correcting on `sample`
therefore removes the injury response together with the batch effect.

`scripts/clustering.batch_diagnostics()` decomposes PC variance by `condition`,
`replicate` and `sample`. `replicate` is the informative one — replicates within a
condition differ only technically. If replicate structure exceeds condition
structure, integration is warranted; set `cfg.USE_INTEGRATION = True`, install
`harmonypy`, and show the UMAP before and after.

Either way the confound remains, and the report should say so rather than claim
the batch effect was removed.

## Numbers from the paper, for comparison only

Do not copy these into the Results section. They belong in the Discussion, as
comparison points.

| Quantity | Bise et al. |
|---|---|
| Cells removed as low-quality / doublets | ~12 ± 2% per condition |
| Clusters | 17 (before merging bipolar and GABAergic subclusters) |
| Retinal ganglion cells | essentially not captured; no cluster with unique `rbpms2b` |
| Rod cluster A (immature-like) | 3.5% of control cells → 22 / 28 / 23% at 3 / 7 / 10 dpMNU |
| Rod A vs B DEGs | ~400 |
| EGFP⁺ cells | 0.54% (ctrl) → 5.40% (3dp) → 2.64% (7dp) → 2.66% (10dp) |
| EGFP⁺ within Müller glia | ~10% of MG cells |
| Genes up in EGFP⁺ vs EGFP⁻ MG | 42 |
| Rod capture bias | authors attribute the higher rod fraction after MNU to easier release of rods from a damaged outer nuclear layer, not to more rods |

## Figure numbering in the paper

The project guide's reading list points to "Figures 1 & 2" for the cell atlas and
"Figures 3 & 4" for EGFP kinetics. In the published article those are
immunofluorescence figures. The scRNA-seq content is:

- **Figure 5** — experimental design, UMAP atlas, canonical marker dot plot, per-timepoint cluster percentages
- **Figure 6** — rod clusters A/B, paralog inversion, GO terms, volcano plots
- **Figure 7** — UV vs non-UV cones
- **Figure 8** — EGFP⁺ cell distribution, counts per timepoint, EGFP⁺ vs EGFP⁻ Müller glia heatmap
- **Figure 9–10** — TOR signalling and rapamycin (no scRNA-seq counterpart in this project)

Cite the real figure numbers.

## Things this dataset cannot show

- **Recovery.** The series stops at 10 dpMNU. Structural restoration in the paper
  is at 30 dpMNU (rods, outer plexiform layer) and 40 dpMNU (UV cone outer
  segments), by immunofluorescence. No transcriptomic statement about recovery is
  available here.
- **Protein-level reporter dynamics.** EGFP transcript ≠ EGFP protein. The
  immunofluorescence time course runs to 90 dpMNU and reflects protein, which
  persists differently from mRNA.
- **Population-level significance.** With n = 2 replicates per condition, a
  cell-level Wilcoxon test treats correlated cells as independent observations. A
  pseudobulk design over the eight samples would be the correct test. Report the
  DE rankings as descriptive.
- **Lineage.** "Immature" and "activated" are transcriptional descriptions of a
  single snapshot, not demonstrated trajectories.
