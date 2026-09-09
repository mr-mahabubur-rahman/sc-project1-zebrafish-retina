# Independent Scanpy reanalysis of *careg:EGFP* reporter dynamics and Müller glia activation in the regenerating zebrafish retina after MNU-induced photoreceptor ablation

Single-cell RNA-seq reanalysis of adult zebrafish retina across an MNU
photoreceptor-injury time course — control, 3, 7 and 10 days post-injury, two
replicates each.

**Primary reference:** Bise T, Pfefferli C, Bonvin M, Taylor L, Lischer HEL,
Bruggmann R, Jaźwińska A (2023). The regeneration-responsive element *careg*
monitors activation of Müller glia after MNU-induced damage of photoreceptors in
the zebrafish retina. *Front. Mol. Neurosci.* 16:1160707.
doi:10.3389/fnmol.2023.1160707 · **Data:** GEO
[GSE202212](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202212)

## The question, and why reanalyse

Zebrafish regenerate lost photoreceptors by reprogramming Müller glia into
proliferative progenitors — a capacity mammals lack. Bise *et al.* showed that the
regeneration-responsive element *careg* acts as a reporter of that activation. This
repository asks a different question of the same data: **which of their conclusions
survive a change of analytical method?**

Re-executing published code proves only that a pipeline is deterministic.
Recovering the same conclusion through independent choices — Scanpy instead of
Seurat, log-normalisation instead of SCTransform, Leiden instead of SNN, no
integration instead of anchor-based integration — tests whether a finding depends
on the decisions every analyst must make. Where results agree, the agreement is
informative. Where they diverge, the divergence identifies which conclusions are
method-dependent.

## Headline results

| | Published | This reanalysis |
|---|---|---|
| EGFP⁺ cells, control | 0.54% | **0.56%** |
| EGFP⁺ cells, 3 dpMNU | 5.40% | **5.40%** |
| EGFP⁺ cells, 7 dpMNU | 2.64% | **2.45%** |
| EGFP⁺ cells, 10 dpMNU | 2.66% | **2.80%** |
| EGFP⁺ within Müller glia | ~10% | **11.34%** (3.75× enriched) |
| Published EGFP⁺ MG signature | 42 genes up | **23 of 23 tested recovered** |
| Cone *arr3a*/*arr3b* inversion | present | **reproduced** |
| Rod immature/mature paralog inversion | present | **not reproduced** |

**14 of 18 reproducible panels reproduced.** The reporter time course agrees within
0.2 percentage points at every timepoint despite a completely independent pipeline.

**The rod result does not reproduce under any treatment tested.** Two sensitivity
analyses eliminate the most likely explanations: applying Harmony batch correction
leaves the verdict unchanged (step 08), as does re-clustering at the published
resolution of 0.2 (step 09). Normalisation remains the leading explanation.

Thirteen cell types were recovered, including horizontal cells (712 cells,
identified by *cx52.6*) and a small retinal ganglion cell population (97 cells,
*rbpms2b*⁺) that the published analysis reported as essentially absent.

## Group members

| | Name | Affiliation | Email |
|---|---|---|---|
| **Group Leader** | Mahabubur Rahman | North South University & National University | mr.mahabub@gmail.com |
| Member | Chandrani Dey | University of Dhaka | cds01307@gmail.com |
| Member | Md Abir Hasan Siddique Rahi | BRAC University | abirhasan2525@gmail.com |
| Member | Md Ariful Amin | University of Dhaka | ariful.amin12@gmail.com |
| Member | Anik Kumar Saha | University of Dhaka | anik.ksaha21@gmail.com |
| Member | Md Sohel Rahman | Gazipur Agricultural University | sohel5926@stu.gau.edu.bd |
| Member | Anamika Jahan Tuli | East West University | jahananamika001@gmail.com |
| Member | Khadija Al Ferdous | Hajee Mohammad Danesh Science and Technology University | azfat2017@gmail.com |
| Member | Md Tangimul Islam | — | rifat.tangimul@gmail.com |

Correspondence to the group leader.

## Quickstart

```bash
git clone https://github.com/mr-mahabubur-rahman/sc-project1-zebrafish-retina.git
cd sc-project1-zebrafish-retina

uv venv --python 3.11
source .venv/bin/activate                 # Windows: .venv\Scripts\Activate.ps1
uv pip install -r requirements.txt
python -m ipykernel install --user --name retina --display-name "Python 3 (retina)"

# download the eight 10x matrices into data/ -- see data/README.md
uv run jupyter lab
```

Run the notebooks in order. Each loads its input from a checkpoint in `results/`
and writes its output there, so there is no hidden state between them and any
notebook can be re-run on its own.

| Notebook | Does | Writes |
|---|---|---|
| `step_0.0_Environment_Management_with_uv.ipynb` | environment, kernel, version record | `docs/session_info.json` |
| `step_01_Loading_Libraries_and_10x_Data.ipynb` | load 8 matrices, metadata, EGFP presence check | `results/01_loaded.h5ad` |
| `step_02_Quality_Control_and_Cell_Filtering.ipynb` | QC metrics, threshold diagnostics, filtering | `results/02_qc_filtered.h5ad`, Figure 1 |
| `step_03_Normalization_HVG_and_Scaling.ipynb` | normalise, HVG, regress, scale | `results/03_preprocessed.h5ad` |
| `step_04_PCA_UMAP_and_Clustering.ipynb` | PCA, batch diagnostics, UMAP, Leiden | `results/04_clustered.h5ad`, Figure 2 |
| `step_05_Cell_Type_Annotation_and_Marker_Discovery.ipynb` | markers, annotation, proportions | `results/05_annotated.h5ad`, Figures 3 and 5 |
| `step_06_EGFP_and_Muller_Glia_Dynamics.ipynb` | EGFP, MG sub-states, rods, cones, depth stratification | `results/06_mg_egfp.h5ad`, Figure 4 |
| `step_07_Paper_Figure_Reproduction.ipynb` | Bise et al. Fig 5, 6, 7, 8 panels | `results/figures/paper_figure_reproduction/` |
| `step_08_Integration_Sensitivity.ipynb` | Harmony correction, compared against the primary result | `results/tables/integration_check/` |
| `step_09_Resolution_Sensitivity.ipynb` | rod paralog test at the published resolution | `results/tables/resolution_check/` |

Steps 08 and 09 are **sensitivity analyses**, not part of the primary pipeline.
They test two analytical decisions rather than replacing them, and require steps
01–05 to have been run. Step 08 additionally needs `harmonypy`.

Notebooks are generated from `tools/build_notebooks.py`, so their content is
reviewable as plain text. `python tools/build_notebooks.py --force` regenerates
them, discarding manual edits — put lasting changes in the builder or in
`scripts/`.

## Repository structure

```
sc-project1-zebrafish-retina/
  README.md
  pyproject.toml  requirements.txt  uv.lock  .gitignore  LICENSE
  data/README.md            download instructions       [matrices git-ignored]
  notebooks/                the ten analysis notebooks
  scripts/
    config.py               paths, parameters, colours, marker panels
    io_utils.py             discovery, loading, metadata, EGFP detection
    qc.py                   QC metrics, threshold diagnostics, filtering
    preprocessing.py        layers, normalisation, HVG, scaling
    clustering.py           PCA, batch diagnostics, UMAP, Leiden
    annotation.py           marker discovery, annotation record
    egfp_analysis.py        EGFP positivity, MG sub-states, rod/cone analysis
    plotting.py             style, PNG+PDF saving, the five project figures
    paper_figures.py        Bise et al. panel equivalents, GO enrichment
    integration_check.py    Harmony sensitivity analysis (step 08)
    resolution_check.py     clustering resolution sensitivity (step 09)
    colour_check.py         colour-vision accessibility verification
  tools/build_notebooks.py
  results/
    tables/                 CSV outputs                 [committed]
    figures/                PNG + PDF                   [committed]
    *.h5ad                  checkpoints                 [git-ignored]
    integrated/             step 08 corrected object    [git-ignored]
  report/
    Rahman_Project1_Report.pdf
    Rahman_Project1_Report.docx
  peer_review/
    review_received.md      peer review from Group 4
    rebuttal_response.md    point-by-point response
  docs/
    method_comparison.md            paper vs guide vs this pipeline
    ai_usage_disclosure.md          appendix template
    INTEGRATION_CHECK_README.md     how to run and read step 08
    validation_checklist.md
    git_workflow.md
    session_info.json               package versions    [generated]
    prompts/                        AI prompt records
```

**Not committed:** the 10x matrices (large, available from GEO) and the `.h5ad`
checkpoints (regenerable, and large enough to breach GitHub's 100 MB file limit).

## Research questions and where they are answered

| # | Question | Notebook | Output |
|---|---|---|---|
| 1 | Major retinal cell types | 05 | Figure 3, `cluster_annotations.csv` |
| 2 | Composition across ctrl → 3dp → 7dp → 10dp | 05 | Figure 5, `cell_type_proportions.csv` |
| 3 | EGFP populations and dynamics | 06 | Figure 4, `egfp_summary.csv` |
| 4 | Müller glia sub-states | 06 | Figure 4F, `mg_substate_evidence.csv` |
| 5 | Rod and cone injury response | 06 | `rod_subcluster_evidence.csv`, `*_injury_response_de.csv` |

## Parameters

Baseline values come from the project guide and are the defaults in
`scripts/config.py`: `min_genes=200`, `min_cells=3`, `pct_mt<15`,
`target_sum=1e4`, 2,000 HVGs, `n_pcs=20`, `n_neighbors=15`, Leiden
`resolution=0.6`, Wilcoxon for DE at `padj<0.05` and `|log2FC|≥0.25`.

To change one: edit `config.py`, add an entry to `cfg.PARAM_OVERRIDES` recording
the old value, the new value and the evidence, and justify it in the notebook
markdown and the Methods section. Nothing is changed silently.

Three deviations are documented in the code:

- **`regress_out` runs after HVG subsetting** rather than on all genes. Same two
  covariates (`total_counts`, `pct_counts_mt`), same result on the retained genes,
  hours faster (`cfg.REGRESS_ON_HVG_SUBSET`).
- **Batch integration is off by default**, decided from diagnostics rather than
  copied from the paper, because condition and sequencing run are confounded in
  this design. **Step 08 tests this empirically:** Harmony removes 79.6% of the
  condition variance against only 30.4% of the technical variance, confirming that
  correcting on `sample` removes the injury response preferentially.
- **Clustering resolution 0.2 gives 12 clusters here against the published 17**
  (step 09). Leiden and Seurat's SNN optimise different objectives, so matching a
  published resolution value does not match a published granularity.

## Figure accessibility

The categorical palette is drawn from colour-vision-deficiency-safe schemes
(Okabe-Ito; Tol, 2021). The assignment of colour to cell type was optimised to
maximise the smallest perceptual distance between any two categories in CAM02-UCS,
evaluated under normal vision, deuteranopia, protanopia and tritanopia
simultaneously.

This was not assumed. Testing the previous palette found a real defect: rods and
retinal ganglion cells sat 2.3 apart under protanopia — effectively the same
colour, and both appear on the same UMAP. The minimum is now 8.8, and no pair falls
below 5. Reproduce the verification with:

```bash
uv run python scripts/colour_check.py
```

Thirteen categories cannot be separated by hue alone under any palette, so every
figure carrying categorical colour also identifies populations by direct label or
by a legend ordered to match the plot. Continuous scales use perceptually uniform
colour maps.

## Reproducibility

- One seed (`cfg.RANDOM_SEED = 0`) passed to PCA, UMAP, Leiden and gene scoring.
- Package versions captured to `docs/session_info.json`.
- No absolute paths; everything resolves from the repository root.
- Every notebook loads its input from a checkpoint rather than relying on
  variables left over from another notebook.
- **All differential expression and marker discovery run on the full 22,813-gene
  matrix in `.raw`**, not the 2,000-gene highly variable subset used for
  clustering. This matters: *cx52.6* is not a highly variable gene, and horizontal
  cells could not have been identified from the reduced set.

Not fully deterministic: UMAP coordinates can differ across platforms and BLAS
builds even with a fixed seed. Cluster membership is stable; exact coordinates may
not be. Harmony (step 08) introduces further run-to-run variation.

## Known limitations

Stated here as well as in the report, because they bound what can be concluded:

- **No doublet detection.** The paper filtered ~12 ± 2% of cells as doublets or
  dying. Cluster 1 (15 cells) and rod subcluster 3 (23 cells) are evident
  multiplets caught by clustering rather than by filtering.
- **Not SCTransform.** Log-normalisation is a different estimator; cluster
  boundaries will not match the paper's exactly. This is the leading remaining
  explanation for the rod result that did not reproduce.
- **Condition and 10x run are confounded by design.** No integration method
  resolves this — step 08 demonstrates the consequence directly.
- **Cell proportions are capture-biased.** Over-capture of rods mechanically
  depresses every other population, so declines are not interpretable; only
  increases against that dilution are.
- **Cell-level DE overstates significance** with n = 2 replicates per condition;
  pseudobulk would be the correct design. Directional skew in DE tracks
  sequencing-depth imbalance between compared groups almost exactly.
- **No recovery claim is possible** — the series ends at day 10, while the paper's
  structural restoration is at 30–40 dpMNU by immunofluorescence.

## AI usage

Parts of this pipeline were AI-assisted. Appendix A of the report discloses the
tools and models used, the prompts submitted, the division of authorship, and the
errors identified and corrected during the work. Prompt records are in
`docs/prompts/`.

## Licence

MIT (see `LICENSE`). The primary paper is CC BY; the sequencing data are from GEO
GSE202212 and are subject to their own terms.
