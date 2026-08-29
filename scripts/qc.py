"""Quality control for zebrafish retina 10x data.

Zebrafish gene symbols are lower case, so the prefixes differ from the human
convention: mitochondrial genes start with 'mt-', ribosomal protein genes with
'rps'/'rpl', haemoglobins with 'hba'/'hbb'/'hb-'.

Filtering is deliberately split from diagnostics: `annotate_qc` and
`qc_summary_table` never modify the data, so thresholds can be judged on the
observed distributions before anything is discarded.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from . import config as cfg


def annotate_qc(adata: ad.AnnData) -> ad.AnnData:
    """Flag mito/ribo/haemoglobin gene families and compute per-cell QC metrics."""
    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    adata.var["ribo"] = adata.var_names.str.startswith(("rps", "rpl"))
    adata.var["hb"] = adata.var_names.str.startswith(("hba", "hbb", "hb-"))

    n_mt, n_ribo, n_hb = (int(adata.var[c].sum()) for c in ("mt", "ribo", "hb"))
    print(f"Gene families flagged: {n_mt} mitochondrial, {n_ribo} ribosomal, {n_hb} haemoglobin")
    if n_mt == 0:
        print(
            "WARNING: no gene matched the 'mt-' prefix. pct_counts_mt will be 0 for "
            "every cell and the mitochondrial filter will do nothing. Inspect "
            "adata.var_names for the actual mitochondrial naming before filtering."
        )

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], percent_top=None, log1p=False, inplace=True
    )
    return adata


def qc_summary_table(adata: ad.AnnData, stage: str) -> pd.DataFrame:
    """Per-sample QC summary at a named stage ('before_filter' / 'after_filter')."""
    grouped = adata.obs.groupby("sample", observed=True)
    df = pd.DataFrame({
        "n_cells": grouped.size(),
        "median_genes": grouped["n_genes_by_counts"].median(),
        "median_counts": grouped["total_counts"].median(),
        "median_mt_percent": grouped["pct_counts_mt"].median(),
        "median_ribo_percent": grouped["pct_counts_ribo"].median(),
        "median_hb_percent": grouped["pct_counts_hb"].median(),
    }).reset_index()
    df.insert(1, "condition", df["sample"].map(
        adata.obs.drop_duplicates("sample").set_index("sample")["condition"]
    ))
    df.insert(0, "stage", stage)
    return df.round(2)


def threshold_diagnostics(adata: ad.AnnData) -> pd.DataFrame:
    """How many cells each threshold would remove, per sample, before filtering.

    Run this BEFORE `filter_cells_and_genes` so the thresholds are a decision, not
    an accident. The percentages are computed independently per criterion, so they
    can overlap and need not sum to the total loss.
    """
    obs = adata.obs
    fails = pd.DataFrame({
        "sample": obs["sample"],
        "fail_min_genes": obs["n_genes_by_counts"] < cfg.MIN_GENES_PER_CELL,
        "fail_max_mt": obs["pct_counts_mt"] >= cfg.MAX_PCT_MT,
    })
    if cfg.MAX_COUNTS_PER_CELL is not None:
        fails["fail_max_counts"] = obs["total_counts"] > cfg.MAX_COUNTS_PER_CELL
    if cfg.MAX_GENES_PER_CELL is not None:
        fails["fail_max_genes"] = obs["n_genes_by_counts"] > cfg.MAX_GENES_PER_CELL

    fail_cols = [c for c in fails.columns if c.startswith("fail_")]
    grouped = fails.groupby("sample", observed=True)
    out = grouped[fail_cols].sum()
    out["n_cells"] = grouped.size()
    for col in fail_cols:
        out[col.replace("fail_", "pct_")] = (100 * out[col] / out["n_cells"]).round(2)
    out["fail_any"] = grouped.apply(
        lambda g: g[fail_cols].any(axis=1).sum(), include_groups=False
    )
    out["pct_any"] = (100 * out["fail_any"] / out["n_cells"]).round(2)
    return out.reset_index()


def interpret_thresholds(diagnostics: pd.DataFrame) -> list[str]:
    """Plain-language flags on whether the baseline thresholds look reasonable.

    These are prompts for the analyst to look at the distributions, not automatic
    decisions. Nothing here changes a parameter.
    """
    notes: list[str] = []
    worst = diagnostics.loc[diagnostics["pct_any"].idxmax()]
    notes.append(
        f"Largest loss: {worst['sample']} at {worst['pct_any']:.1f}% of cells failing "
        "at least one criterion."
    )
    if (diagnostics["pct_any"] > 30).any():
        heavy = diagnostics.loc[diagnostics["pct_any"] > 30, "sample"].tolist()
        notes.append(
            f"Over 30% of cells would be removed in {heavy}. That is high enough to "
            "check whether the threshold is wrong for these samples rather than the "
            "cells being bad -- inspect the violin plots before accepting it."
        )
    if "pct_max_mt" in diagnostics and (diagnostics["pct_max_mt"] > 20).any():
        notes.append(
            "The mitochondrial cut alone removes >20% of cells somewhere. Injured "
            "retina dissociates harder than control, so a condition-dependent loss "
            "here is itself a result worth reporting, not just a technical nuisance."
        )
    spread = diagnostics["pct_any"].max() - diagnostics["pct_any"].min()
    if spread > 20:
        notes.append(
            f"Loss varies by {spread:.1f} percentage points across samples; check "
            "whether the variation tracks condition (biological) or replicate "
            "(technical) before interpreting cell-type proportions."
        )
    return notes


def filter_cells_and_genes(
    adata: ad.AnnData,
    min_genes: int | None = None,
    min_cells: int | None = None,
    max_pct_mt: float | None = None,
) -> tuple[ad.AnnData, pd.DataFrame]:
    """Apply the baseline QC filters and return the filtered object plus a log.

    Order matters: cells first, then genes, so that `min_cells` is evaluated on the
    surviving cells.
    """
    min_genes = cfg.MIN_GENES_PER_CELL if min_genes is None else min_genes
    min_cells = cfg.MIN_CELLS_PER_GENE if min_cells is None else min_cells
    max_pct_mt = cfg.MAX_PCT_MT if max_pct_mt is None else max_pct_mt

    log: list[dict] = []
    n_obs, n_vars = adata.n_obs, adata.n_vars
    log.append({"step": "input", "cells": n_obs, "genes": n_vars})

    sc.pp.filter_cells(adata, min_genes=min_genes)
    log.append({"step": f"min_genes>={min_genes}", "cells": adata.n_obs, "genes": adata.n_vars})

    adata = adata[adata.obs["pct_counts_mt"] < max_pct_mt].copy()
    log.append({"step": f"pct_mt<{max_pct_mt}", "cells": adata.n_obs, "genes": adata.n_vars})

    if cfg.MAX_COUNTS_PER_CELL is not None:
        adata = adata[adata.obs["total_counts"] <= cfg.MAX_COUNTS_PER_CELL].copy()
        log.append({"step": f"counts<={cfg.MAX_COUNTS_PER_CELL}",
                    "cells": adata.n_obs, "genes": adata.n_vars})
    if cfg.MAX_GENES_PER_CELL is not None:
        adata = adata[adata.obs["n_genes_by_counts"] <= cfg.MAX_GENES_PER_CELL].copy()
        log.append({"step": f"genes<={cfg.MAX_GENES_PER_CELL}",
                    "cells": adata.n_obs, "genes": adata.n_vars})

    sc.pp.filter_genes(adata, min_cells=min_cells)
    log.append({"step": f"min_cells>={min_cells}", "cells": adata.n_obs, "genes": adata.n_vars})

    assert adata.n_obs > 0, "QC removed every cell. Re-examine the thresholds."

    log_df = pd.DataFrame(log)
    log_df["cells_lost"] = log_df["cells"].diff().fillna(0).astype(int).abs()
    print(f"Cells: {n_obs} -> {adata.n_obs} ({100 * adata.n_obs / n_obs:.1f}% kept)")
    print(f"Genes: {n_vars} -> {adata.n_vars}")
    return adata, log_df


def per_sample_cell_counts(adata: ad.AnnData) -> pd.DataFrame:
    """Cell counts by sample / condition / replicate, for cell_counts_by_sample.csv."""
    df = (
        adata.obs.groupby(["sample", "condition", "replicate"], observed=True)
        .size().reset_index(name="n_cells")
    )
    total = df["n_cells"].sum()
    df["pct_of_dataset"] = (100 * df["n_cells"] / total).round(2)
    return df.sort_values("sample")


def egfp_qc_note(adata: ad.AnnData, egfp_name: str | None) -> str:
    """One-line statement of transgene detectability, for the QC section of the report."""
    if egfp_name is None:
        return "EGFP feature absent from the count matrix; reporter analysis not possible."
    counts = np.asarray(adata[:, egfp_name].X.todense()).ravel() \
        if hasattr(adata[:, egfp_name].X, "todense") \
        else np.asarray(adata[:, egfp_name].X).ravel()
    n_pos = int((counts >= cfg.EGFP_POSITIVE_MIN_COUNTS).sum())
    return (
        f"EGFP feature '{egfp_name}': {n_pos} of {adata.n_obs} cells "
        f"({100 * n_pos / adata.n_obs:.2f}%) have >= {cfg.EGFP_POSITIVE_MIN_COUNTS} count(s); "
        f"max {counts.max():.0f} counts in a single cell."
    )
