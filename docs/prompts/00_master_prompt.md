# MASTER PROMPT — Reproducible Scanpy Workflow for Zebrafish Retina Regeneration scRNA-seq Project

## Role

Act as a **senior computational biologist, bioinformatician, single-cell RNA-seq analyst, Scanpy/AnnData expert, zebrafish retina regeneration researcher, scientific programmer, and reproducible-research engineer**.

I am completing a student single-cell RNA-seq project on **adult zebrafish (Danio rerio) retinal regeneration following MNU-induced injury**, using careg:EGFP reporter data.

Your task is to design and write a **complete, reproducible, scientifically defensible Python/Scanpy workflow for JupyterLab** that I can execute from start to finish and commit to GitHub.

Do not give me a generic tutorial. Build the workflow specifically around the project structure, research questions, expected figures, and primary paper described below.

---

# 1. PROJECT OBJECTIVE

Analyze 8 10x Genomics scRNA-seq samples:

* ctrl1 = uninjured control replicate 1
* ctrl2 = uninjured control replicate 2
* 3dp1 = 3 days post-MNU injury replicate 1
* 3dp2 = 3 days post-MNU injury replicate 2
* 7dp1 = 7 days post-MNU injury replicate 1
* 7dp2 = 7 days post-MNU injury replicate 2
* 10dp1 = 10 days post-MNU injury replicate 1
* 10dp2 = 10 days post-MNU injury replicate 2

The primary biological objectives are:

1. Identify the major retinal cell types.
2. Determine how cellular composition changes across Ctrl → 3dp → 7dp → 10dp.
3. Determine which cell populations express EGFP and characterize EGFP dynamics.
4. Investigate Müller glia activation/reprogramming and distinguish resting MG, activated MG, and proliferative/progenitor states where supported by the data.
5. Investigate changes in rod and cone photoreceptors after MNU injury.
6. Reproduce the instructor-required figures in a clean, publication-quality, reproducible manner.

---

# 2. PROJECT DIRECTORY STRUCTURE

Assume the repository is:

project1_zebrafish_retina/

Use this structure:

project1_zebrafish_retina/
│
├── data/
│   ├── ctrl1/
│   │   ├── filtered_feature_bc_matrix/
│   │   └── web_summary.html
│   ├── ctrl2/
│   │   ├── filtered_feature_bc_matrix/
│   │   └── web_summary.html
│   ├── 3dp1/
│   │   ├── filtered_feature_bc_matrix/
│   │   └── web_summary.html
│   ├── 3dp2/
│   ├── 7dp1/
│   ├── 7dp2/
│   ├── 10dp1/
│   └── 10dp2/
│
├── notebooks/
│   ├── step_0.0_Environment_Management_with_uv.ipynb
│   ├── step_01_Loading_Libraries_and_10x_Data.ipynb
│   ├── step_02_Quality_Control_and_Cell_Filtering.ipynb
│   ├── step_03_Normalization_HVG_and_Scaling.ipynb
│   ├── step_04_PCA_UMAP_and_Clustering.ipynb
│   ├── step_05_Cell_Type_Annotation_and_Marker_Discovery.ipynb
│   └── step_06_EGFP_and_Muller_Glia_Dynamics.ipynb
│
├── scripts/
│   ├── config.py
│   ├── io_utils.py
│   ├── qc.py
│   ├── preprocessing.py
│   ├── clustering.py
│   ├── annotation.py
│   ├── egfp_analysis.py
│   └── plotting.py
│
├── tables/
│   ├── qc_summary.csv
│   ├── cell_counts_by_sample.csv
│   ├── marker_genes.csv
│   ├── cluster_annotations.csv
│   ├── egfp_summary.csv
│   └── cell_type_proportions.csv
│
├── figures/
│   ├── figure_01_qc/
│   ├── figure_02_umap/
│   ├── figure_03_markers/
│   ├── figure_04_egfp_mg/
│   └── figure_05_cell_proportions/
│
├── README.md
├── requirements.txt
├── pyproject.toml
├── .gitignore
└── LICENSE

Do not change this structure unless there is a strong technical reason. If a change is necessary, explain why before implementing it.

---

# 3. INPUT DATA REQUIREMENTS

The data are standard 10x Genomics filtered feature-barcode matrices:

* barcodes.tsv.gz
* features.tsv.gz
* matrix.mtx.gz

The feature table may contain the EGFP reporter.

Automatically discover the 8 sample directories rather than hard-coding absolute paths.

Use pathlib.Path and make the workflow runnable from the repository root.

Never use absolute local paths such as:

C:/Users/...
/home/user/...
/mnt/data/...

Use relative paths such as:

Path("data")

---

# 4. ENVIRONMENT

Use **uv** for reproducible Python environment management.

Target:

* Python 3.11
* Scanpy
* AnnData
* NumPy
* pandas
* SciPy
* matplotlib
* seaborn
* Jupyter
* ipykernel
* gseapy where useful

Do not unnecessarily introduce additional packages.

If a package is needed, explain its purpose.

Provide:

1. pyproject.toml
2. requirements.txt
3. commands for creating the environment
4. commands for installing dependencies
5. commands for registering the Jupyter kernel
6. commands for running the workflow

Make the environment reproducible.

---

# 5. IMPORTANT SCIENTIFIC CONSTRAINT

The project guide provides a baseline Scanpy workflow, including:

* minimum genes = 200
* minimum cells per gene = 3
* mitochondrial percentage < 15%
* normalization target = 10,000 counts/cell
* 2,000 HVGs
* PCA
* 15 nearest neighbors
* approximately 20 PCs for neighbors
* Leiden resolution = 0.6

Use these as the **initial baseline parameters**.

Do not silently change these values.

If you believe a parameter should be changed because the actual data quality or Scanpy version requires it:

1. explain why,
2. show the diagnostic evidence,
3. make the parameter configurable,
4. preserve the baseline as the default unless evidence supports changing it.

---

# 6. IMPORTANT PAPER-REPRODUCTION PRINCIPLE

The primary paper used Seurat and describes:

* filtering based on low library size, low expressed genes and mitochondrial reads,
* SCTransform normalization,
* removal of mitochondrial/ribosomal confounding,
* integration,
* k-nearest-neighbor/SNN clustering,
* UMAP,
* marker discovery,
* EGFP-positive-cell subsetting,
* differential expression.

My project, however, must be implemented in **Python/Scanpy**.

Therefore:

DO NOT pretend that Scanpy's standard log-normalization pipeline is identical to Seurat SCTransform.

Explicitly distinguish:

A. What the published paper did.

B. What the student project guide specifies.

C. What our Scanpy implementation actually does.

If an exact methodological reproduction of the paper is impossible in Scanpy without additional tools, state that clearly rather than claiming exact reproduction.

---

# 7. SAMPLE METADATA

Create a clean metadata table containing at minimum:

* sample
* condition
* replicate
* timepoint
* batch

Use:

ctrl → Control
3dp → 3 days
7dp → 7 days
10dp → 10 days

Preserve sample identity throughout the analysis.

Make categorical ordering explicit:

condition_order = ["ctrl", "3dp", "7dp", "10dp"]

Never allow alphabetical sorting to change the biological time order.

---

# 8. STEP 01 — LOAD 10X DATA

Create a notebook that:

1. imports libraries,
2. configures Scanpy,
3. discovers sample directories,
4. reads each 10x matrix,
5. validates the expected files,
6. detects gene names,
7. checks whether EGFP exists,
8. adds metadata,
9. makes gene names unique,
10. concatenates the eight datasets into one AnnData object.

Use informative assertions and error messages.

Report for every sample:

* number of cells
* number of genes
* total counts
* presence/absence of EGFP

Save an initial metadata/QC table.

---

# 9. STEP 02 — QUALITY CONTROL

Implement QC appropriate for zebrafish retina.

Identify:

* mitochondrial genes: mt-
* ribosomal genes: rps-, rpl-
* hemoglobin genes: hba, hbb, hb-

Calculate:

* total_counts
* n_genes_by_counts
* pct_counts_mt
* pct_counts_ribo
* pct_counts_hb

Generate:

1. violin plots by condition,
2. violin plots by sample,
3. scatter plot of total counts vs detected genes,
4. mitochondrial percentage diagnostic,
5. pre/post-filter cell counts.

Use the baseline thresholds:

* min_genes = 200
* min_cells = 3
* pct_counts_mt < 15%

Make thresholds configurable.

Before filtering, calculate and save a QC summary table.

After filtering, calculate another summary table.

Do not perform filtering blindly. Include a diagnostic interpretation explaining whether the thresholds appear reasonable for these data.

---

# 10. STEP 03 — NORMALIZATION, HVG, AND SCALING

Preserve raw counts:

adata.layers["counts"]

Use the project guide's baseline workflow:

1. normalize_total(target_sum=1e4)
2. log1p
3. preserve normalized data appropriately
4. identify 2,000 HVGs using sample as batch key
5. regress technical covariates only if justified
6. scale

Be careful with AnnData/raw/layers so that downstream differential expression can still access appropriate expression values.

Explain the difference between:

* raw counts
* normalized counts
* log-transformed expression
* scaled expression

Do not overwrite raw counts.

Generate the HVG diagnostic figure.

---

# 11. STEP 04 — PCA, NEIGHBORS, UMAP, LEIDEN

Perform:

* PCA
* PCA variance explained diagnostic
* neighbor graph
* UMAP
* Leiden clustering

Baseline:

n_neighbors = 15
n_pcs = 20
resolution = 0.6

However, make these configurable.

Generate:

### Figure 2A

UMAP colored by Leiden cluster.

### Figure 2B

UMAP colored by condition.

### Figure 2C

UMAP colored by sample.

The UMAP must preserve sample/condition information sufficiently to allow assessment of batch effects.

Also generate cluster-size summaries.

Explain whether clustering appears biologically plausible.

---

# 12. BATCH / REPLICATE AWARENESS

Do not confuse biological timepoint with technical replicate.

Replicates are:

ctrl1/ctrl2
3dp1/3dp2
7dp1/7dp2
10dp1/10dp2

Always retain sample identity.

If there is evidence of strong batch effects, diagnose them before applying correction.

Do not automatically apply Harmony, BBKNN, scVI, or another integration method merely because integration is mentioned in the paper.

If integration is used, make it an explicit, reproducible analysis choice and explain:

* why it is needed,
* what method is used,
* which representation is integrated,
* how the result differs from the uncorrected analysis.

The biological signal of interest is the injury time course, so avoid overcorrection that removes real biological differences.

---

# 13. STEP 05 — CELL TYPE ANNOTATION

Use marker-based annotation.

Start with the marker genes provided by the project guide.

### Müller Glia

* rlbp1a
* glula
* slc1a3b

### Activated Müller Glia / progenitor-like

* EGFP
* pcna
* mki67
* ascl1a
* her4.1

### Rod photoreceptors

* rho
* gnat1
* nr2e3

### Cone photoreceptors

* opn1sw1
* opn1sw2
* opn1mw1
* opn1lw1

### Bipolar cells

* vsx1
* vsx2
* islet1

### Amacrine cells

* gad1b
* gad2
* tfap2a

### Retinal ganglion cells

* rbpms2b
* isl2b
* pou4f1

### Microglia / immune cells

* mpeg1.1
* coro1a
* cxcr4b

Before plotting, verify which genes actually exist in adata.var_names.

Never invent replacement genes.

If a marker is absent, report it as unavailable and continue with available markers.

---

# 14. MARKER DISCOVERY

Use:

sc.tl.rank_genes_groups(..., method="wilcoxon")

Generate:

1. ranked marker plot,
2. top marker table per cluster,
3. dotplot of canonical markers,
4. heatmap of selected markers.

Export marker results to:

tables/marker_genes.csv

Do not define cell types solely from one marker.

Require multiple concordant markers where possible.

For every cluster, create an annotation table:

* cluster
* proposed cell type
* positive markers
* conflicting markers
* confidence
* biological rationale

Use labels such as:

high / medium / low confidence

where appropriate.

---

# 15. IMPORTANT: DO NOT FORCE ANNOTATIONS

If the data do not support a cell-type assignment, say:

"Unresolved / candidate X"

rather than forcing a label.

Distinguish:

* cell identity markers
* activation markers
* proliferation markers
* stress markers
* injury-response genes

EGFP alone must not be treated as proof that a cell is Müller glia.

---

# 16. FIGURE 3 — MARKER-BASED CELL TYPE ANNOTATION

Create a publication-quality figure containing:

A. Marker dotplot.

B. Marker heatmap.

C. UMAP colored by manually assigned cell type.

Use consistent cell-type colors throughout the entire project.

Do not assign arbitrary colors independently in different notebooks.

Save the color dictionary in:

scripts/config.py

---

# 17. STEP 06 — EGFP DYNAMICS

This is a major biological objective.

First verify whether the gene is named exactly:

EGFP

If not, inspect the feature table and determine the actual feature name.

Do not assume capitalization.

Report:

* EGFP detection status,
* number of EGFP-positive cells,
* fraction EGFP-positive by sample,
* fraction EGFP-positive by condition,
* expression distribution by condition,
* expression by cell type/cluster.

Define EGFP positivity transparently.

Do not use an arbitrary threshold without justification.

Prefer a clearly documented threshold such as expression > 0 if appropriate for sparse UMI counts, while also reporting continuous expression.

---

# 18. EGFP FIGURES

Create:

### Figure 4A

EGFP expression on UMAP.

### Figure 4B

EGFP expression distribution across Ctrl, 3dp, 7dp, 10dp.

### Figure 4C

Percentage of EGFP-positive cells across timepoints.

### Figure 4D

EGFP expression by annotated cell type.

### Figure 4E

EGFP vs proliferation markers:

* pcna
* mki67
* ascl1a
* her4.1

Only plot genes that actually exist.

Use violin, dotplot, feature plot, or equivalent appropriate visualization.

---

# 19. MÜLLER GLIA DEEP DIVE

Subset Müller glia based on defensible marker evidence.

Then investigate:

1. resting/quiescent MG,
2. activated MG,
3. proliferative/progenitor-like MG.

Compare:

* EGFP
* pcna
* mki67
* ascl1a
* her4.1
* rlbp1a
* glula
* slc1a3b

If the data do not support three distinct populations, do not manufacture them.

Perform a focused MG re-analysis if justified.

Generate:

* MG-only UMAP,
* MG marker dotplot,
* EGFP expression,
* proliferation-marker expression,
* condition distribution.

---

# 20. ROD ANALYSIS

The paper describes distinct immature and mature rod populations and identifies markers including:

* rho/rhol
* gnat1
* gngt1
* rom1b
* meig1
* phototransduction genes

Use the data to investigate rod heterogeneity.

Do not automatically reproduce the paper's cluster labels.

Instead:

1. identify rod clusters in the dataset,
2. determine their marker genes,
3. identify immature/maturing vs mature-like transcriptional programs,
4. compare their abundance across Ctrl, 3dp, 7dp, 10dp,
5. report whether the data support the paper's immature/mature rod interpretation.

Clearly separate observed results from interpretation.

---

# 21. CONE ANALYSIS

Investigate cone populations using available opsin markers:

* opn1sw1
* opn1sw2
* opn1mw1
* opn1lw1

Determine:

* whether cone subtypes are detected,
* their abundance across timepoints,
* whether their expression changes following MNU injury.

Do not claim recovery if the dataset only covers 10 days and does not directly measure later recovery.

---

# 22. CELL-TYPE PROPORTION ANALYSIS

Calculate cell-type proportions by:

* sample
* condition
* replicate

Do not only calculate percentages from pooled cells without showing replicate-level information.

Generate:

### Figure 5

Cell-type proportions across regeneration timepoints.

Prefer a visualization that clearly shows:

Ctrl → 3dp → 7dp → 10dp

If possible, show replicate-level observations in addition to pooled proportions.

Explicitly warn that scRNA-seq cell proportions can be affected by:

* cell survival,
* tissue dissociation,
* capture efficiency,
* sequencing depth,
* changes in cell size/adhesion,
* biological differences.

This is particularly important because the primary paper itself discusses increased rod representation after MNU injury as potentially influenced by survival/capture bias.

---

# 23. REQUIRED FINAL FIGURES

The final project must produce at least:

## Figure 1 — QC

* n_genes
* total_counts
* mitochondrial percentage
* ribosomal percentage
* filtering summary

## Figure 2 — Global Cell Atlas

* UMAP by cluster
* UMAP by condition
* UMAP by sample

## Figure 3 — Cell Type Annotation

* marker dotplot
* marker heatmap
* annotated UMAP

## Figure 4 — EGFP / Müller Glia Dynamics

* EGFP UMAP
* EGFP across timepoints
* EGFP-positive fraction
* MG markers
* proliferation markers

## Figure 5 — Cell-Type Dynamics

* cell-type proportions across Ctrl → 3dp → 7dp → 10dp
* replicate-aware representation where possible

Each figure must be saved as both:

* PNG
* PDF

Use high resolution.

---

# 24. PUBLICATION-QUALITY PLOTTING

Create a central plotting configuration.

Use:

* consistent figure dimensions,
* consistent font sizes,
* readable legends,
* descriptive axis labels,
* biological timepoint ordering,
* consistent cell-type colors,
* vector PDF output where appropriate.

Every plotting function should have:

* title,
* axis labels,
* legend,
* deterministic output filename.

Do not rely on interactive-only plots.

---

# 25. OUTPUT TABLES

Automatically generate:

### qc_summary.csv

Sample-level:

* cells_before
* cells_after
* median_genes
* median_counts
* median_mt_percent
* median_ribo_percent

### marker_genes.csv

* cluster
* gene
* score
* logfoldchange
* pval
* pval_adj

### cluster_annotations.csv

* cluster
* cell_type
* confidence
* rationale

### egfp_summary.csv

* sample
* condition
* cell_type
* n_cells
* egfp_positive
* egfp_positive_fraction
* mean_egfp
* median_egfp

### cell_type_proportions.csv

* sample
* condition
* replicate
* cell_type
* cell_count
* proportion

---

# 26. REPRODUCIBILITY

Set random seeds wherever supported.

Record:

* Python version
* Scanpy version
* AnnData version
* NumPy version
* pandas version

Store these in the analysis output.

Avoid hidden state between notebooks.

Each notebook should either:

A. load a clearly defined intermediate AnnData object, or

B. explicitly execute the necessary previous steps.

Do not assume that variables from another notebook already exist.

---

# 27. INTERMEDIATE DATA

Save appropriate intermediate AnnData objects, for example:

results/
├── 01_loaded.h5ad
├── 02_qc_filtered.h5ad
├── 03_preprocessed.h5ad
├── 04_clustered.h5ad
├── 05_annotated.h5ad
└── 06_mg_egfp.h5ad

However, DO NOT commit huge generated data files to GitHub unless explicitly required.

Configure .gitignore appropriately.

Explain which files should remain local and which should be committed.

---

# 28. CODE QUALITY REQUIREMENTS

Every notebook should contain:

1. Markdown explaining the biological purpose.
2. Code.
3. Expected output.
4. Diagnostic interpretation.
5. A short "What to check before continuing" section.

Use functions instead of repeatedly copying code.

Use type hints where useful.

Use descriptive variable names.

Avoid unnecessary one-line magic code.

Add assertions for important assumptions.

Examples:

assert adata.n_obs > 0
assert "condition" in adata.obs
assert "sample" in adata.obs

---

# 29. ERROR HANDLING

Anticipate common Scanpy problems:

* EGFP absent or differently named
* duplicate gene symbols
* sparse matrix handling
* Leiden package missing
* incompatible Scanpy/NumPy versions
* insufficient PCs
* missing marker genes
* zero-expression genes
* empty subsets
* plotting failures
* categorical ordering problems

When writing code, include informative error messages.

Do not silently catch errors.

---

# 30. BIOLOGICAL VALIDATION

At the end of the workflow, create a concise validation report answering:

1. Are major retinal cell types detected?
2. Are Müller glia identifiable?
3. Are EGFP-positive cells concentrated in the expected regenerative population?
4. Does EGFP change across the injury time course?
5. Do pcna/mki67 support proliferative activity?
6. Are rods altered after MNU injury?
7. Are cone populations detectable?
8. Does the cell-type composition change across time?
9. Are replicate effects larger or smaller than biological timepoint effects?
10. Which observations agree with the primary paper?
11. Which observations cannot be directly compared?
12. Which conclusions are uncertain?

Do not manufacture biological conclusions if the data do not support them.

---

# 31. PAPER COMPARISON

Use the supplied paper as the primary biological reference.

The paper reports that careg:EGFP is induced in activated Müller glia following MNU injury and describes its temporal activation during regeneration.

It also reports that EGFP-positive cells overlap with proliferative markers during early stages, while later proliferating progenitor populations can be EGFP-negative.

The paper's scRNA-seq Figure 5 specifically presents integrated UMAP clusters, canonical marker expression, and cluster percentages across timepoints.

Use these observations as comparison points, NOT as results that should be copied into the student's analysis.

---

# 32. ACADEMIC INTEGRITY

The project guide explicitly requires disclosure of AI assistance, including the AI tool/model, exact prompts, purpose of each prompt, and validation notes.

Therefore:

* Do not hide AI involvement.
* Do not fabricate validation.
* Do not claim that code was independently written if it was AI-assisted.
* Include an AI Usage Disclosure section in README/report documentation.
* Preserve the exact prompt used for major code-generation steps.
* Explain what I must personally verify.

The final workflow must be understandable enough that I can explain every major line during an oral defense.

---

# 33. GIT / GITHUB REQUIREMENTS

The final project must be ready to push to GitHub.

The GitHub repository/account target is:

mr-mahabubur-rahman

Do not assume the exact repository name unless I provide it.

Prepare:

* README.md
* .gitignore
* requirements.txt
* pyproject.toml
* notebooks
* scripts
* tables
* figures
* reproducibility instructions

Provide exact commands for:

git init
git add
git commit
git branch
git remote
git push

But DO NOT execute destructive Git commands.

Do not use:

git push --force

unless I explicitly request it.

Before recommending git add, identify files/directories that should NOT be committed, particularly raw/large sequencing data and generated intermediate objects.

---

# 34. NOTEBOOK DESIGN

For each notebook, provide complete Jupyter-ready code cells.

Use this format:

### Markdown Cell

Purpose and scientific explanation.

### Code Cell

Complete executable Python code.

### Markdown Cell

What the output means.

### Code Cell

Next analysis step.

Do not provide pseudocode where executable code is requested.

Do not omit imports.

Do not assume variables exist from previous notebooks unless the notebook explicitly loads them.

---

# 35. EXPECTED RESPONSE FORMAT

First provide:

## A. Final project architecture

Show the directory tree.

Then:

## B. Environment setup

Provide the uv commands.

Then:

## C. Notebook-by-notebook workflow

For:

1. step_0.0
2. step_01
3. step_02
4. step_03
5. step_04
6. step_05
7. step_06

provide the complete code.

Then:

## D. Reusable scripts

Provide the contents of the necessary scripts.

Then:

## E. Expected figures

Map each figure to the notebook/code that creates it.

Then:

## F. Expected tables

Map each table to the code that creates it.

Then:

## G. Validation checklist

Provide a checklist I can use before considering the analysis complete.

Then:

## H. GitHub preparation

Provide .gitignore, README structure, and safe Git commands.

Then:

## I. AI Usage Disclosure

Provide a template that I can include in my final report.

---

# 36. MOST IMPORTANT RULE

Do not optimize for producing the largest amount of code.

Optimize for:

**correctness → reproducibility → biological validity → interpretability → figure reproduction → GitHub readiness.**

Whenever a scientific decision is uncertain, stop and identify the uncertainty rather than silently guessing.

Whenever a gene is missing, report that it is missing rather than inventing a substitute.

Whenever a result differs from the paper, investigate whether the difference is due to:

* preprocessing,
* filtering,
* normalization,
* integration,
* clustering,
* annotation,
* biological sampling,
* or technical/capture effects.

Clearly label:

**Observed result**

vs.

**Biological interpretation**

vs.

**Comparison with published paper**

vs.

**Speculation / limitation**.

The final workflow should be suitable for a student research project, reproducible in JupyterLab, scientifically defensible during oral examination, and ready for GitHub submission.
