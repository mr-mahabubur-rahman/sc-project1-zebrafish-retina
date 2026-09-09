"""Central configuration for the zebrafish retina regeneration scRNA-seq project.

Every notebook and script imports from here so that paths, parameters, colours and
category ordering are defined exactly once. Nothing in this file reads the data;
it only declares constants.

Parameter provenance is labelled explicitly:

  [GUIDE]    value prescribed by the student project guide -> baseline default
  [PROJECT]  our own choice, made because the guide is silent on the point
  [PAPER]    value used by Bise et al. 2023 (Seurat), recorded for comparison only

Baseline [GUIDE] values are never changed silently. To deviate, override the
constant in the notebook, record the override in `PARAM_OVERRIDES`, and justify it
in the notebook markdown.
"""

from __future__ import annotations

from pathlib import Path

# --------------------------------------------------------------------------- #
# Paths (all relative to the repository root, so the workflow is portable)
# --------------------------------------------------------------------------- #

REPO_ROOT: Path = Path(__file__).resolve().parents[1]

DATA_DIR: Path = REPO_ROOT / "data"
RESULTS_DIR: Path = REPO_ROOT / "results"
# Tables and figures live under results/ so the repository follows the layout in
# the peer-review guidelines. Both are committed; only the .h5ad checkpoints and
# results/integrated/ are git-ignored.
TABLES_DIR: Path = RESULTS_DIR / "tables"
FIGURES_DIR: Path = RESULTS_DIR / "figures"
DOCS_DIR: Path = REPO_ROOT / "docs"
REPORT_DIR: Path = REPO_ROOT / "report"
PEER_REVIEW_DIR: Path = REPO_ROOT / "peer_review"

FIGURE_DIRS: dict[str, Path] = {
    "qc": FIGURES_DIR / "figure_01_qc",
    "umap": FIGURES_DIR / "figure_02_umap",
    "markers": FIGURES_DIR / "figure_03_markers",
    "egfp_mg": FIGURES_DIR / "figure_04_egfp_mg",
    "proportions": FIGURES_DIR / "figure_05_cell_proportions",
}

# Intermediate AnnData checkpoints, one per notebook. These are regenerable and
# are therefore git-ignored (see .gitignore).
H5AD = {
    "loaded": RESULTS_DIR / "01_loaded.h5ad",
    "qc_filtered": RESULTS_DIR / "02_qc_filtered.h5ad",
    "preprocessed": RESULTS_DIR / "03_preprocessed.h5ad",
    "clustered": RESULTS_DIR / "04_clustered.h5ad",
    "annotated": RESULTS_DIR / "05_annotated.h5ad",
    "mg_egfp": RESULTS_DIR / "06_mg_egfp.h5ad",
}

MATRIX_SUBDIR: str = "filtered_feature_bc_matrix"
REQUIRED_10X_FILES: tuple[str, ...] = (
    "barcodes.tsv.gz",
    "features.tsv.gz",
    "matrix.mtx.gz",
)

# --------------------------------------------------------------------------- #
# Experimental design
# --------------------------------------------------------------------------- #

EXPECTED_SAMPLES: tuple[str, ...] = (
    "ctrl1", "ctrl2", "3dp1", "3dp2", "7dp1", "7dp2", "10dp1", "10dp2",
)

# Biological time order. Never rely on alphabetical sorting: "10dp" sorts before
# "3dp" as a string, which would silently reverse the time course in every plot.
CONDITION_ORDER: list[str] = ["ctrl", "3dp", "7dp", "10dp"]

CONDITION_LABELS: dict[str, str] = {
    "ctrl": "Control",
    "3dp": "3 days",
    "7dp": "7 days",
    "10dp": "10 days",
}

# Days post-MNU as a numeric covariate (control = 0). Useful for ordering and for
# any trend test; it does not imply the control is "day 0 post-injury".
CONDITION_DAYS: dict[str, int] = {"ctrl": 0, "3dp": 3, "7dp": 7, "10dp": 10}

# --------------------------------------------------------------------------- #
# Baseline analysis parameters
# --------------------------------------------------------------------------- #

RANDOM_SEED: int = 0

# --- QC --------------------------------------------------------------------- #
MIN_GENES_PER_CELL: int = 200        # [GUIDE]
MIN_CELLS_PER_GENE: int = 3          # [GUIDE]
MAX_PCT_MT: float = 15.0             # [GUIDE]
# Deliberately not applied by default: an upper bound on counts/genes is a common
# doublet heuristic, but the guide does not ask for it and a hard cap would remove
# real high-RNA cells. Set to None to disable; a number to enable and document.
MAX_COUNTS_PER_CELL: int | None = None   # [PROJECT]
MAX_GENES_PER_CELL: int | None = None    # [PROJECT]

# --- Normalisation / feature selection -------------------------------------- #
TARGET_SUM: float = 1e4              # [GUIDE]
N_TOP_HVG: int = 2000                # [GUIDE]
HVG_BATCH_KEY: str = "sample"        # [GUIDE]

# Technical covariate regression. Two covariates are regressed out with
# sc.pp.regress_out: total_counts and pct_counts_mt.
#
# The guide's walkthrough calls regress_out on the full gene set. We keep the same
# covariates and keep the step on by default (baseline preserved) but run it after
# subsetting to HVGs, because regressing ~30k genes x ~16k cells is hours of
# compute for genes that are then discarded before PCA. The covariates and the
# result on the retained genes are identical. This is a documented deviation, not
# a silent one.
REGRESS_OUT: bool = True                                          # [GUIDE]
REGRESS_COVARIATES: list[str] = ["total_counts", "pct_counts_mt"]  # [GUIDE]
REGRESS_ON_HVG_SUBSET: bool = True                                # [PROJECT]
SCALE_MAX_VALUE: float = 10.0        # [GUIDE]

# --- Dimensionality reduction / clustering ---------------------------------- #
N_PCS_NEIGHBORS: int = 20            # [GUIDE] "approximately 20", confirm on elbow plot
N_NEIGHBORS: int = 15                # [GUIDE]
LEIDEN_RESOLUTION: float = 0.6       # [GUIDE]
LEIDEN_KEY: str = "leiden_res_0.6"   # [GUIDE]
# Additional resolutions computed for a stability check only. The reported
# clustering stays at LEIDEN_RESOLUTION unless a change is justified in writing.
LEIDEN_RESOLUTION_SWEEP: list[float] = [0.2, 0.4, 0.6, 0.8, 1.0]

# --- Batch handling --------------------------------------------------------- #
# OFF by default, on purpose. The paper integrated (Seurat, dims = 30); we only
# integrate if the Step 04 diagnostics show replicate-driven rather than
# condition-driven structure. The injury time course IS the biological signal, so
# over-correction is a real risk here. Step 08 tests this decision empirically:
# Harmony removes 79.6% of the condition variance against 30.4% of the technical
# variance, confirming the risk is real for this design.
USE_INTEGRATION: bool = False        # [PROJECT]
INTEGRATION_METHOD: str = "harmony"  # only consulted when USE_INTEGRATION is True
INTEGRATION_BATCH_KEY: str = "sample"

# --- Differential expression ------------------------------------------------ #
DE_METHOD: str = "wilcoxon"          # [GUIDE] and [PAPER] use the same test
DE_TOP_N: int = 50                   # markers exported per cluster
# Significance thresholds applied when a gene is called differentially expressed.
# Benjamini-Hochberg adjusted p-value, and an absolute log2 fold-change floor so
# that no result rests on a statistically significant but trivially small effect.
DE_PADJ_CUTOFF: float = 0.05         # [PROJECT]
DE_LOG2FC_CUTOFF: float = 0.25       # [PROJECT]

# --- EGFP ------------------------------------------------------------------- #
# The feature name is NOT assumed. io_utils.find_egfp_feature() searches var_names
# case-insensitively for these candidates and reports what it found.
EGFP_NAME_CANDIDATES: tuple[str, ...] = (
    "EGFP", "eGFP", "egfp", "GFP", "gfp", "EGFP-transgene", "careg-EGFP", "careg_EGFP",
)
# Positivity rule. With sparse 3' UMI data and a single-copy transgene, >0 counts
# is the defensible default; anything higher is an arbitrary cut unless justified
# from the observed count distribution. Continuous expression is always reported
# alongside the binary call.
EGFP_POSITIVE_MIN_COUNTS: int = 1    # [PROJECT] raw counts, not normalised

# Record here any parameter you override in a notebook, e.g.
# PARAM_OVERRIDES["MAX_PCT_MT"] = (15.0, 10.0, "long tail of high-mito rods, see Fig 1C")
PARAM_OVERRIDES: dict[str, tuple] = {}

# --------------------------------------------------------------------------- #
# Marker panels
# --------------------------------------------------------------------------- #

# Panel exactly as supplied by the project guide. Genes are checked against
# adata.var_names before use; missing genes are reported, never substituted.
MARKERS_GUIDE: dict[str, list[str]] = {
    "Muller glia (resting)": ["rlbp1a", "glula", "slc1a3b"],
    "Activated MG / progenitor": ["pcna", "mki67", "ascl1a", "her4.1"],
    "Rod photoreceptors": ["rho", "gnat1", "nr2e3"],
    "Cone photoreceptors": ["opn1sw1", "opn1sw2", "opn1mw1", "opn1lw1"],
    "Bipolar cells": ["vsx1", "vsx2", "islet1"],
    "Amacrine cells": ["gad1b", "gad2", "tfap2a"],
    "Retinal ganglion cells": ["rbpms2b", "isl2b", "pou4f1"],
    "Microglia / immune": ["mpeg1.1", "coro1a", "cxcr4b"],
}

# Supplementary panel drawn from Bise et al. 2023 (Fig. 5D and the rod/cone
# sections) plus the cell types the guide's panel does not cover. Used as
# secondary evidence for clusters the guide panel leaves unresolved.
MARKERS_PAPER: dict[str, list[str]] = {
    "Muller glia": ["glulb", "six3b", "gfap", "apoeb", "atp1a1b"],
    "Rods (mature-like)": ["rhol", "pde6gb", "guca1b", "sagb", "gnb1a", "rcvrna"],
    "Rods (immature-like)": ["rho", "pde6ga", "guca1a", "meig1", "ppdpfa", "rom1b", "gngt1"],
    "Cones (UV)": ["opn1sw1", "arr3b", "tbx2a", "cngb3.2", "guca1e"],
    "Cones (non-UV)": ["arr3a", "opn1mw1", "opn1lw1", "opn1sw2"],
    "Horizontal cells": ["rem1", "cx52.6"],
    "Bipolar cells": ["vsx1", "cabp5a", "prkca"],
    "Oligodendrocytes": ["mbpa", "plp1b"],
    "Erythrocytes": ["hbba1", "hbaa1"],
    "Pericytes": ["pdgfrb", "myl9a"],
    "Microglia": ["mpeg1.1", "apoc1"],
    "RPE": ["rpe65a", "bhlhe40"],
}

# Genes plotted together in the EGFP / regenerative-activation figure.
ACTIVATION_MARKERS: list[str] = ["pcna", "mki67", "ascl1a", "her4.1", "stat3", "mdka"]
MG_IDENTITY_MARKERS: list[str] = ["rlbp1a", "glula", "glulb", "slc1a3b", "gfap", "apoeb"]

# --------------------------------------------------------------------------- #
# Colours — defined once, used everywhere, so a cell type keeps one colour in
# every figure of the report.
# --------------------------------------------------------------------------- #

# Colours are drawn from colour-vision-deficiency-safe qualitative schemes (Paul
# Tol's muted, light and vibrant sets, and Okabe-Ito). The assignment of colour to
# cell type was not chosen by eye: it was optimised to maximise the smallest
# perceptual distance (CAM02-UCS) between any two categories, evaluated
# simultaneously under normal vision, deuteranopia, protanopia and tritanopia.
#
# The previous palette placed Rods and Retinal ganglion cells 2.3 apart under
# protanopia — effectively the same colour, and both appear on the same UMAP. The
# minimum across all categories and all three deficiency forms is now 8.8.
#
# Entries match exactly the cell types the annotation produces. Reproduce the
# verification with:  uv run python scripts/colour_check.py
CELL_TYPE_COLORS: dict[str, str] = {
    "Muller glia":            "#999933",
    "Rods":                   "#009988",
    "Cones":                  "#0077BB",
    "Bipolar cells":          "#CC3311",
    "Amacrine cells":         "#332288",
    "Horizontal cells":       "#882255",
    "Retinal ganglion cells": "#AA4499",
    "Microglia":              "#44BB99",
    "Oligodendrocytes":       "#EEDD88",
    "Erythrocytes":           "#EE8866",
    "Pericytes":              "#117733",
    "RPE":                    "#CC6677",
    "Rods (low quality)":     "#99DDFF",
    "Unresolved":             "#BBBBBB",
}

CONDITION_COLORS: dict[str, str] = {
    "ctrl": "#4D4D4D",
    "3dp": "#F4A582",
    "7dp": "#D6604D",
    "10dp": "#B2182B",
}

# --------------------------------------------------------------------------- #
# Figure defaults
# --------------------------------------------------------------------------- #

FIGURE_DPI: int = 300
FIGURE_FORMATS: tuple[str, ...] = ("png", "pdf")
BASE_FONT_SIZE: int = 9


def ensure_directories() -> None:
    """Create every output directory the workflow writes to."""
    for path in (RESULTS_DIR, TABLES_DIR, FIGURES_DIR, DOCS_DIR,
                 *FIGURE_DIRS.values()):
        path.mkdir(parents=True, exist_ok=True)
