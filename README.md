# Zebrafish retina regeneration â€” single-cell RNA-seq analysis

Reanalysis of `careg:EGFP` adult zebrafish retina scRNA-seq across an MNU
photoreceptor-injury time course (control, 3, 7 and 10 days post-injury; two
replicates each).

Primary reference: Bise T, Pfefferli C, Bonvin M, Taylor L, Lischer HEL,
Bruggmann R, JaÅºwiÅ„ska A (2023). *The regeneration-responsive element careg
monitors activation of MÃ¼ller glia after MNU-induced damage of photoreceptors in
the zebrafish retina.* Front. Mol. Neurosci. 16:1160707.
doi:10.3389/fnmol.2023.1160707. Data: GEO **GSE202212**.

This repository implements the analysis in **Python/Scanpy**. The paper used
Seurat with SCTransform and anchor-based integration. The two are **not**
equivalent, and `docs/method_comparison.md` sets out exactly where they diverge.

---

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

Run the notebooks in order. Each one loads its input from `results/` and writes
its output there, so there is no hidden state between them and any notebook can be
re-run on its own.

| Notebook | Does | Writes |
|---|---|---|
| `step_0.0_Environment_Management_with_uv.ipynb` | environment, kernel, version record | `docs/session_info.json` |
| `step_01_Loading_Libraries_and_10x_Data.ipynb` | load 8 matrices, metadata, EGFP presence check | `results/01_loaded.h5ad` |
| `step_02_Quality_Control_and_Cell_Filtering.ipynb` | QC metrics, threshold diagnostics, filtering | `results/02_qc_filtered.h5ad`, Figure 1 |
| `step_03_Normalization_HVG_and_Scaling.ipynb` | normalise, HVG, regress, scale | `results/03_preprocessed.h5ad` |
| `step_04_PCA_UMAP_and_Clustering.ipynb` | PCA, batch diagnostics, UMAP, Leiden | `results/04_clustered.h5ad`, Figure 2 |
| `step_05_Cell_Type_Annotation_and_Marker_Discovery.ipynb` | markers, panel scores, annotation, proportions | `results/05_annotated.h5ad`, Figures 3 and 5 |
| `step_06_EGFP_and_Muller_Glia_Dynamics.ipynb` | EGFP, MG sub-states, rods, cones, validation | `results/06_mg_egfp.h5ad`, Figure 4 |

Notebooks are generated from `tools/build_notebooks.py`, so their content is
reviewable as plain text. `python tools/build_notebooks.py --force` regenerates
them (this discards any manual edits â€” put lasting changes in the builder or in
`scripts/`).

---

## Repository structure

```
project1_zebrafish_retina/
â”œâ”€â”€ data/                 8 x 10x filtered matrices          [git-ignored]
â”œâ”€â”€ notebooks/            the seven analysis notebooks
â”œâ”€â”€ scripts/              all reusable logic
â”‚   â”œâ”€â”€ config.py         paths, parameters, colours, marker panels
â”‚   â”œâ”€â”€ io_utils.py       discovery, loading, metadata, EGFP detection, checkpoints
â”‚   â”œâ”€â”€ qc.py             QC metrics, threshold diagnostics, filtering
â”‚   â”œâ”€â”€ preprocessing.py  layers, normalisation, HVG, scaling
â”‚   â”œâ”€â”€ clustering.py     PCA, batch diagnostics, integration, UMAP, Leiden
â”‚   â”œâ”€â”€ annotation.py     marker discovery, panel scoring, annotation record
â”‚   â”œâ”€â”€ egfp_analysis.py  EGFP positivity, MG sub-states, rod/cone analysis
â”‚   â””â”€â”€ plotting.py       style, PNG+PDF saving, the five figures
â”œâ”€â”€ tools/build_notebooks.py
â”œâ”€â”€ tables/               CSV results                        [committed]
â”œâ”€â”€ figures/              PNG + PDF                          [committed]
â”œâ”€â”€ results/              .h5ad checkpoints                  [git-ignored]
â”œâ”€â”€ docs/
â”‚   â”œâ”€â”€ method_comparison.md      paper vs guide vs this pipeline
â”‚   â”œâ”€â”€ ai_usage_disclosure.md    appendix template
â”‚   â””â”€â”€ session_info.json         package versions            [generated]
â”œâ”€â”€ pyproject.toml Â· requirements.txt Â· .gitignore Â· LICENSE
```

**Committed:** notebooks, scripts, tables, figures, docs, environment files.
**Not committed:** `data/` (large, redistributable from GEO) and `results/*.h5ad`
(regenerable, and large enough to breach GitHub's 100 MB file limit).

---

## Research questions and where they are answered

| # | Question | Notebook | Output |
|---|---|---|---|
| 1 | Major retinal cell types | 05 | Figure 3, `cluster_annotations.csv` |
| 2 | Composition across ctrl â†’ 3dp â†’ 7dp â†’ 10dp | 05 | Figure 5, `cell_type_proportions.csv` |
| 3 | EGFPâº populations and dynamics | 06 | Figure 4, `egfp_summary.csv` |
| 4 | MÃ¼ller glia sub-states | 06 | Figure 4F, `mg_substate_evidence.csv` |
| 5 | Rod and cone injury response | 06 | `rod_subcluster_evidence.csv`, `*_injury_response_de.csv` |

## Figures

| Figure | Panels | Built by | Directory |
|---|---|---|---|
| 1 | QC violins by sample, counts-vs-genes, filter waterfall | `plotting.figure_01_qc` (nb 02) | `figures/figure_01_qc/` |
| 2 | UMAP by cluster / condition / sample, plus split-by-condition | `plotting.figure_02_umap`, `figure_02_split_by_condition` (nb 04) | `figures/figure_02_umap/` |
| 3 | marker dotplot, matrixplot, annotated UMAP | `plotting.figure_03_*` (nb 05) | `figures/figure_03_markers/` |
| 4 | EGFP UMAP, by timepoint, positive fraction, by cell type, activation markers, MG sub-states | `plotting.figure_04_*` (nb 06) | `figures/figure_04_egfp_mg/` |
| 5 | composition per sample, proportion trajectories with replicate range | `plotting.figure_05_proportions` (nb 05) | `figures/figure_05_cell_proportions/` |
| S1â€“S2 | HVG diagnostic, PCA elbow | `plotting.figure_hvg`, `figure_pca_variance` | `figures/figure_01_qc/`, `figure_02_umap/` |

Every figure is written as **PNG (300 dpi) and vector PDF**, with editable text
(`pdf.fonttype = 42`).

## Tables

| File | Contents | Notebook |
|---|---|---|
| `qc_summary.csv` | per-sample cells/genes/counts/mito/ribo, before and after filtering | 02 |
| `cell_counts_by_sample.csv` | cells by sample, condition, replicate | 02 |
| `batch_diagnostics.csv` | PC variance explained by condition / replicate / sample | 04 |
| `leiden_resolution_sweep.csv` | cluster counts across resolutions 0.2â€“1.0 | 04 |
| `cluster_composition.csv` | sample composition of each cluster | 04 |
| `marker_genes.csv` | ranked markers per cluster (Wilcoxon) | 05 |
| `cluster_annotations.csv` | cluster, cell type, markers, confidence, rationale | 05 |
| `cell_type_proportions.csv` | proportions per sample and replicate | 05 |
| `egfp_summary.csv` | EGFPâº counts and fractions by sample Ã— condition Ã— cell type | 06 |
| `egfp_enrichment_by_cell_type.csv` | EGFPâº enrichment per cell type | 06 |
| `mg_substate_evidence.csv` | MG subcluster marker and condition profile | 06 |
| `rod_subcluster_evidence.csv` / `cone_subcluster_evidence.csv` | photoreceptor subcluster profiles | 06 |
| `*_injury_response_de.csv` | per-cell-type DE, each timepoint vs control | 06 |
| `biological_validation.csv` | final validation answers | 06 |

---

## Parameters

Baseline values come from the project guide and are the defaults in
`scripts/config.py`: `min_genes=200`, `min_cells=3`, `pct_mt<15`,
`target_sum=1e4`, 2,000 HVGs, `n_pcs=20`, `n_neighbors=15`, Leiden `resolution=0.6`,
Wilcoxon for DE.

To change one: edit `config.py`, add an entry to `cfg.PARAM_OVERRIDES` recording
the old value, the new value and the evidence, and justify it in the notebook
markdown and the Methods section. Nothing is changed silently.

Two deviations are already documented in the code:

1. `regress_out` runs after HVG subsetting rather than on all genes â€” same
   covariates, same result, hours faster (`cfg.REGRESS_ON_HVG_SUBSET`).
2. Batch integration is **off by default** and decided from diagnostics, because
   condition and sequencing run are confounded in this design
   (`cfg.USE_INTEGRATION`).

## Reproducibility

- One seed (`cfg.RANDOM_SEED = 0`) passed to PCA, UMAP, Leiden and gene scoring.
- Package versions captured to `docs/session_info.json` in Step 0.0.
- No absolute paths anywhere; everything resolves from the repository root.
- Every notebook loads its input from a checkpoint rather than relying on
  variables left over from another notebook.

**Not fully deterministic:** UMAP coordinates can differ slightly across
platforms and BLAS builds even with a fixed seed. Cluster membership is stable;
exact embedding coordinates may not be.

## Known limitations

Stated here as well as in the report, because they bound what can be concluded:

1. **No doublet detection** (the paper filtered ~12 Â± 2% of cells as doublets or
   dying). Clusters with mixed identity markers may be doublets.
2. **Not SCTransform.** Log-normalisation is a different estimator; cluster
   boundaries will not match the paper's exactly.
3. **Condition and 10x run are confounded** by design; no integration method
   resolves this.
4. **Cell proportions are capture-biased.** The paper attributes the post-MNU rod
   increase to easier release of rods from a damaged outer nuclear layer.
   Compositional data also mean one type's over-capture depresses all others.
5. **Cell-level DE overstates significance** with n = 2 replicates per condition;
   pseudobulk would be the correct design. Rankings are descriptive.
6. **No recovery claim is possible** â€” the series ends at day 10, while the
   paper's structural restoration is at 30â€“40 dpMNU by immunofluorescence.

## AI usage

Parts of this pipeline were AI-assisted. `docs/ai_usage_disclosure.md` is the
appendix template required by the course academic-integrity policy: tool, model,
verbatim prompts, purpose, what was modified, and what was personally verified.

## Licence

MIT (see `LICENSE`). The primary paper is CC BY; the sequencing data are from GEO
GSE202212 and are subject to their own terms.

