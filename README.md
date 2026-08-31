# Zebrafish retina regeneration single-cell RNA-seq analysis

Reanalysis of *careg:EGFP* adult zebrafish retina scRNA-seq across an MNU
photoreceptor-injury time course (control, 3, 7 and 10 days post-injury; two
replicates each).

**Primary reference:** Bise T, Pfefferli C, Bonvin M, Taylor L, Lischer HEL,
Bruggmann R, Jaźwińska A (2023). The regeneration-responsive element *careg*
monitors activation of Müller glia after MNU-induced damage of photoreceptors in
the zebrafish retina. *Front. Mol. Neurosci.* 16:1160707.
doi:10.3389/fnmol.2023.1160707. Data: GEO **GSE202212**.

This repository implements the analysis in Python/Scanpy. The paper used Seurat
with SCTransform and anchor-based integration. The two are not equivalent, and
`docs/method_comparison.md` sets out exactly where they diverge.

## Headline results

| | Published | This reanalysis |
|---|---|---|
| EGFP⁺ cells, control | 0.54% | **0.56%** |
| EGFP⁺ cells, 3 dpMNU | 5.40% | **5.40%** |
| EGFP⁺ cells, 7 dpMNU | 2.64% | **2.45%** |
| EGFP⁺ cells, 10 dpMNU | 2.66% | **2.80%** |
| EGFP⁺ within Müller glia | ~10% | **11.34%** |
| Published EGFP⁺ MG signature genes | 42 up | **23 of 23 tested recovered** |
| Cone *arr3a*/*arr3b* inversion | present | **reproduced** |
| Rod immature/mature paralog inversion | present | **not reproduced** |

14 of 18 reproducible panels reproduced. The reporter time course agrees within
0.2 percentage points at every timepoint despite a completely independent
pipeline. The rod result does not reproduce under any tested treatment, and two
sensitivity analyses (steps 08 and 09) eliminate batch structure and clustering
granularity as explanations.

## Quick start

```bash
git clone https://github.com/mr-mahabubur-rahman/sc-project1-zebrafish-retina.git
cd sc-project1-zebrafish-retina

uv venv --python 3.11
source .venv/bin/activate                 # Windows: .venv\Scripts\Activate.ps1
uv pip install -r requirements.txt
python -m ipykernel install --user --name retina --display-name "Python 3 (retina)"

# place the eight 10x sample folders under data/ (see data/README.md)

uv run jupyter lab
```

Run the notebooks in order. Each loads its input from `results/` and writes its
output there, so there is no hidden state between them and any notebook can be
re-run on its own.

| Notebook | Does | Writes |
|---|---|---|
| `step_0.0_Environment_Management_with_uv.ipynb` | environment, kernel, version record | `docs/session_info.json` |
| `step_01_Loading_Libraries_and_10x_Data.ipynb` | load 8 matrices, metadata, EGFP presence check | `results/01_loaded.h5ad` |
| `step_02_Quality_Control_and_Cell_Filtering.ipynb` | QC metrics, threshold diagnostics, filtering | `results/02_qc_filtered.h5ad`, Figure 1 |
| `step_03_Normalization_HVG_and_Scaling.ipynb` | normalise, HVG, regress, scale | `results/03_preprocessed.h5ad` |
| `step_04_PCA_UMAP_and_Clustering.ipynb` | PCA, batch diagnostics, UMAP, Leiden | `results/04_clustered.h5ad`, Figure 2 |
| `step_05_Cell_Type_Annotation_and_Marker_Discovery.ipynb` | markers, annotation, proportions | `results/05_annotated.h5ad`, Figures 3 and 5 |
| `step_06_EGFP_and_Muller_Glia_Dynamics.ipynb` | EGFP, MG sub-states, rods, cones, depth stratification | `results/06_mg_egfp.h5ad`, Figure 4 |
| `step_07_Paper_Figure_Reproduction.ipynb` | Bise et al. Fig 5, 6, 7, 8 panels | `figures/paper_figure_reproduction/` |
| `step_08_Integration_Sensitivity.ipynb` | Harmony correction, compared against the primary result | `tables/integration_check/`, `figures/integration_check/` |
| `step_09_Resolution_Sensitivity.ipynb` | rod paralog test at the published resolution | `tables/resolution_check/`, `figures/resolution_check/` |

Steps 08 and 09 are **sensitivity analyses**, not part of the primary pipeline.
They test two analytical decisions rather than replacing them, and they require
steps 01–05 to have been run. Step 08 additionally requires `harmonypy`
(`uv pip install harmonypy`), which the primary pipeline does not use.

Notebooks are generated from `tools/build_notebooks.py`, so their content is
reviewable as plain text. `python tools/build_notebooks.py --force` regenerates
them — this discards any manual edits, so put lasting changes in the builder or in
`scripts/`.

## Repository structure

```
project1_zebrafish_retina/
  data/                     8 x 10x filtered matrices        [git-ignored]
  notebooks/                the ten analysis notebooks
  scripts/
    config.py               paths, parameters, colours, marker panels
    io_utils.py             discovery, loading, metadata, EGFP detection, checkpoints
    qc.py                   QC metrics, threshold diagnostics, filtering
    preprocessing.py        layers, normalisation, HVG, scaling
    clustering.py           PCA, batch diagnostics, integration, UMAP, Leiden
    annotation.py           marker discovery, panel scoring, annotation record
    egfp_analysis.py        EGFP positivity, MG sub-states, rod/cone analysis
    plotting.py             style, PNG+PDF saving, the five project figures
    paper_figures.py        Bise et al. panel equivalents, GO enrichment
    integration_check.py    Harmony sensitivity analysis (step 08)
    resolution_check.py     clustering resolution sensitivity (step 09)
  tools/build_notebooks.py
  tables/                   CSV results                      [committed]
    integration_check/      step 08 comparison tables
    resolution_check/       step 09 comparison tables
  figures/                  PNG + PDF                        [committed]
    figure_01_qc/ ... figure_05_cell_proportions/
    paper_figure_reproduction/
    integration_check/
    resolution_check/
  results/                  .h5ad checkpoints                [git-ignored]
  results_integrated/       step 08 corrected object         [git-ignored]
  docs/
    method_comparison.md            paper vs guide vs this pipeline
    ai_usage_disclosure.md          appendix template
    INTEGRATION_CHECK_README.md     how to run and read step 08
    session_info.json               package versions         [generated]
  pyproject.toml  requirements.txt  .gitignore  LICENSE
```

**Committed:** notebooks, scripts, tables, figures, docs, environment files.
**Not committed:** `data/` (large, redistributable from GEO) and the `.h5ad`
checkpoints in `results/` and `results_integrated/` (regenerable, and large enough
to breach GitHub's 100 MB file limit).

## Research questions and where they are answered

| # | Question | Notebook | Output |
|---|---|---|---|
| 1 | Major retinal cell types | 05 | Figure 3, `cluster_annotations.csv` |
| 2 | Composition across ctrl → 3dp → 7dp → 10dp | 05 | Figure 5, `cell_type_proportions.csv` |
| 3 | EGFP populations and dynamics | 06 | Figure 4, `egfp_summary.csv` |
| 4 | Müller glia sub-states | 06 | Figure 4F, `mg_substate_evidence.csv` |
| 5 | Rod and cone injury response | 06 | `rod_subcluster_evidence.csv`, `*_injury_response_de.csv` |

Thirteen cell types were recovered, including horizontal cells (712 cells,
identified by *cx52.6*) and a small retinal ganglion cell population (97 cells,
*rbpms2b*⁺) that the published analysis reported as essentially absent.

## Figures

| Figure | Panels | Built by | Directory |
|---|---|---|---|
| 1 | QC violins by sample, counts-vs-genes, filter waterfall | `plotting.figure_01_qc` (nb 02) | `figures/figure_01_qc/` |
| 2 | UMAP by cluster / condition / sample, plus split-by-condition | `plotting.figure_02_*` (nb 04) | `figures/figure_02_umap/` |
| 3 | marker dotplot, matrixplot, annotated UMAP | `plotting.figure_03_*` (nb 05) | `figures/figure_03_markers/` |
| 4 | EGFP UMAP, by timepoint, positive fraction, by cell type, activation markers, MG sub-states | `plotting.figure_04_*` (nb 06) | `figures/figure_04_egfp_mg/` |
| 5 | composition per sample, proportion trajectories with replicate range | `plotting.figure_05_proportions` (nb 05) | `figures/figure_05_cell_proportions/` |
| S1–S2 | HVG diagnostic, PCA elbow | `plotting.figure_hvg`, `figure_pca_variance` | `figures/figure_01_qc/`, `figure_02_umap/` |
| Paper panels | 18 equivalents of Bise et al. Fig 5, 6, 7, 8 | `paper_figures.*` (nb 07) | `figures/paper_figure_reproduction/` |
| Integration | embeddings, variance, confusion heatmap, key results | `integration_check.figure_*` (nb 08) | `figures/integration_check/` |
| Resolution | paralog ratios by clustering resolution | `resolution_check.figure_*` (nb 09) | `figures/resolution_check/` |

Every figure is written as PNG (300 dpi) and vector PDF, with editable text
(`pdf.fonttype = 42`).

## Tables

| File | Contents | Notebook |
|---|---|---|
| `qc_summary.csv` | per-sample cells/genes/counts/mito/ribo, before and after filtering | 02 |
| `cell_counts_by_sample.csv` | cells by sample, condition, replicate | 02 |
| `batch_diagnostics.csv` | PC variance explained by condition / replicate / sample | 04 |
| `leiden_resolution_sweep.csv` | cluster counts across resolutions 0.2–1.0 | 04 |
| `cluster_composition.csv` | sample composition of each cluster | 04 |
| `marker_genes.csv` | ranked markers per cluster (Wilcoxon, full gene set) | 05 |
| `cluster_annotations.csv` | cluster, cell type, markers, confidence, rationale | 05 |
| `cell_type_proportions.csv` | proportions per sample and replicate | 05 |
| `egfp_summary.csv` | EGFP counts and fractions by sample, condition, cell type | 06 |
| `egfp_enrichment_by_cell_type.csv` | EGFP enrichment per cell type | 06 |
| `egfp_depth_quartiles.csv`, `egfp_depth_by_condition.csv` | depth-stratified positivity | 06 |
| `mg_substate_evidence.csv` | MG subcluster marker and condition profile | 06 |
| `rod_subcluster_evidence.csv`, `cone_subcluster_evidence.csv` | photoreceptor subcluster profiles | 06 |
| `*_injury_response_de.csv` | per-cell-type DE, each timepoint vs control | 06 |
| `biological_validation.csv` | final validation answers | 06 |
| `integration_check/*.csv` | variance before/after, cluster agreement, purity, EGFP and paralog comparisons | 08 |
| `resolution_check/*.csv` | rod marker means, paralogs at resolution 0.2, summary | 09 |

## Parameters

Baseline values come from the project guide and are the defaults in
`scripts/config.py`: `min_genes=200`, `min_cells=3`, `pct_mt<15`,
`target_sum=1e4`, 2,000 HVGs, `n_pcs=20`, `n_neighbors=15`, Leiden
`resolution=0.6`, Wilcoxon for DE.

To change one: edit `config.py`, add an entry to `cfg.PARAM_OVERRIDES` recording
the old value, the new value and the evidence, and justify it in the notebook
markdown and the Methods section. Nothing is changed silently.

Two deviations are documented in the code:

- `regress_out` runs after HVG subsetting rather than on all genes — same
  covariates, same result, hours faster (`cfg.REGRESS_ON_HVG_SUBSET`).
- Batch integration is off by default and decided from diagnostics, because
  condition and sequencing run are confounded in this design
  (`cfg.USE_INTEGRATION`). **Step 08 tests this decision empirically:** applying
  Harmony removes 79.6% of the condition variance against only 30.4% of the
  technical variance, confirming that correcting on `sample` removes the injury
  response preferentially.

**One further point on parameters.** Clustering at resolution 0.2 — the value used
in the published analysis — produces 12 clusters here against their 17 (step 09).
Leiden and Seurat's SNN modularity optimisation solve related but distinct
problems, so matching a published resolution value does not match a published
granularity.

## Reproducibility

- One seed (`cfg.RANDOM_SEED = 0`) passed to PCA, UMAP, Leiden and gene scoring.
- Package versions captured to `docs/session_info.json` in Step 0.0.
- No absolute paths anywhere; everything resolves from the repository root.
- Every notebook loads its input from a checkpoint rather than relying on
  variables left over from another notebook.
- **All differential expression and marker discovery run on the full 22,813-gene
  matrix in `.raw`**, not on the 2,000-gene highly variable subset used for
  clustering. This matters: *cx52.6* is not a highly variable gene, and horizontal
  cells could not have been identified from the reduced set.

Not fully deterministic: UMAP coordinates can differ slightly across platforms and
BLAS builds even with a fixed seed. Cluster membership is stable; exact embedding
coordinates may not be. Harmony (step 08) introduces further run-to-run variation
even with a fixed seed.

## Known limitations

Stated here as well as in the report, because they bound what can be concluded:

- **No doublet detection.** The paper filtered ~12 ± 2% of cells as doublets or
  dying. Clusters with mixed identity markers may be doublets; cluster 1 (15
  cells) and rod subcluster 3 (23 cells) are evident multiplets caught by
  clustering rather than by filtering.
- **Not SCTransform.** Log-normalisation is a different estimator; cluster
  boundaries will not match the paper's exactly. This is the leading remaining
  explanation for the rod result that did not reproduce.
- **Condition and 10x run are confounded by design.** No integration method
  resolves this — step 08 demonstrates the consequence directly.
- **Cell proportions are capture-biased.** The paper attributes the post-MNU rod
  increase to easier release of rods from a damaged outer nuclear layer.
  Compositional data also mean one type's over-capture depresses all others, so
  declines are not interpretable; only increases against that dilution are.
- **Cell-level DE overstates significance** with n = 2 replicates per condition;
  pseudobulk would be the correct design. Rankings are descriptive. The practical
  consequence is visible: directional skew in DE tracks sequencing-depth imbalance
  between the compared groups almost exactly.
- **No recovery claim is possible** — the series ends at day 10, while the paper's
  structural restoration is at 30–40 dpMNU by immunofluorescence.

## AI usage

Parts of this pipeline were AI-assisted. `docs/ai_usage_disclosure.md` is the
appendix template required by the course academic-integrity policy: tool, model,
verbatim prompts, purpose, what was modified, and what was personally verified.

## Licence

MIT (see `LICENSE`). The primary paper is CC BY; the sequencing data are from GEO
GSE202212 and are subject to their own terms.
