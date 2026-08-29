"""Normalisation, feature selection and scaling.

Four representations of the same data exist after this step, and confusing them is
the most common source of wrong differential expression results:

  adata.layers["counts"]  raw integer UMI counts. Never overwritten. Input to any
                          count-based model and to the EGFP positivity call.
  adata.layers["lognorm"] library-size normalised (target_sum) then log1p. This is
                          the correct matrix for rank_genes_groups, dotplots,
                          violins and any expression number quoted in the report.
  adata.raw               a frozen copy of the lognorm matrix over ALL genes, so
                          plots can show genes that were not selected as HVGs.
  adata.X (after scaling) z-scored HVGs, clipped. Only for PCA/neighbours/UMAP.
                          Values are unitless and can be negative; never quote them
                          as expression and never use them for DE.

What the paper did instead: SCTransform (regularised negative binomial residuals)
with mitochondrial and ribosomal percentage regressed out during normalisation.
That is a different estimator, not a cosmetic difference, and this pipeline does
not reproduce it. See docs/method_comparison.md.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from . import config as cfg


def layer_names(adata: ad.AnnData) -> list[str]:
    """Sorted layer names, tolerant of the anndata quirk that puts a None in keys().

    anndata 0.13.x returns a spurious ``None`` entry from ``adata.layers.keys()``,
    which makes a bare ``sorted(...)`` raise a TypeError. Membership tests
    (``"counts" in adata.layers``) are unaffected, so only listing needs this.
    """
    return sorted(str(k) for k in adata.layers.keys() if k is not None)


def stash_counts(adata: ad.AnnData) -> ad.AnnData:
    """Store raw counts in a layer before any transformation touches .X."""
    if "counts" in adata.layers:
        print("layers['counts'] already present; leaving it untouched.")
        return adata
    X = adata.X
    values = X.data if hasattr(X, "data") else np.asarray(X).ravel()
    sample = values[:1000]
    if sample.size and not np.allclose(sample, np.round(sample)):
        raise ValueError(
            "adata.X does not look like integer counts, so it has already been "
            "transformed. Re-run from the QC checkpoint rather than normalising twice."
        )
    adata.layers["counts"] = X.copy()
    return adata


def normalize_and_log(adata: ad.AnnData, target_sum: float | None = None) -> ad.AnnData:
    """Library-size normalise to `target_sum` counts per cell, then log1p."""
    target_sum = cfg.TARGET_SUM if target_sum is None else target_sum
    stash_counts(adata)
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    adata.layers["lognorm"] = adata.X.copy()
    adata.raw = adata          # all genes, lognorm, for plotting and DE
    print(f"Normalised to {target_sum:.0f} counts/cell and log1p-transformed; "
          "stored in layers['lognorm'] and .raw")
    return adata


def select_hvgs(
    adata: ad.AnnData,
    n_top_genes: int | None = None,
    batch_key: str | None = None,
) -> ad.AnnData:
    """Select highly variable genes, computed per batch and combined.

    `batch_key` makes the selection robust to genes that are variable in only one
    sample -- important here, because a gene that is variable only in an injured
    sample would otherwise dominate feature selection.
    """
    n_top_genes = cfg.N_TOP_HVG if n_top_genes is None else n_top_genes
    batch_key = cfg.HVG_BATCH_KEY if batch_key is None else batch_key
    if batch_key not in adata.obs:
        raise KeyError(f"HVG batch_key '{batch_key}' is not a column of adata.obs.")

    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, batch_key=batch_key)
    n_hvg = int(adata.var["highly_variable"].sum())
    print(f"Selected {n_hvg} highly variable genes (requested {n_top_genes}).")
    if n_hvg < 0.8 * n_top_genes:
        print("NOTE: fewer HVGs than requested, which happens when many genes are "
              "variable in too few batches. Check the HVG diagnostic plot.")
    return adata


def preserve_gene_counts_in_obs(adata: ad.AnnData, gene: str,
                                obs_key: str | None = None) -> bool:
    """Copy a gene's raw counts into .obs so it survives HVG subsetting.

    Subsetting to highly variable genes removes every other gene from .X and from
    the layers. EGFP is a single-copy transgene and is not guaranteed to be
    selected as an HVG, so its counts are copied into .obs BEFORE the subset. This
    changes no analysis -- it only keeps a value reachable.
    """
    if gene not in adata.var_names:
        return False
    obs_key = obs_key or f"{gene.lower()}_counts"
    layer = "counts" if "counts" in adata.layers else "X"
    adata.obs[obs_key] = get_expression(adata, gene, layer=layer)
    print(f"Preserved raw counts of '{gene}' in obs['{obs_key}'] before HVG subsetting.")
    return True


def regress_and_scale(
    adata: ad.AnnData,
    regress: bool | None = None,
    covariates: list[str] | None = None,
    on_hvg_subset: bool | None = None,
    max_value: float | None = None,
) -> ad.AnnData:
    """Optionally regress out technical covariates, then z-score and clip.

    Deviation from the guide's snippet, applied deliberately and reported: the
    guide regresses over the full gene set; we subset to HVGs first. The genes
    dropped are the ones that never enter PCA, so the analysis is unchanged while
    the runtime falls from hours to minutes. Set `on_hvg_subset=False` to reproduce
    the guide's exact order.
    """
    regress = cfg.REGRESS_OUT if regress is None else regress
    covariates = cfg.REGRESS_COVARIATES if covariates is None else covariates
    on_hvg_subset = cfg.REGRESS_ON_HVG_SUBSET if on_hvg_subset is None else on_hvg_subset
    max_value = cfg.SCALE_MAX_VALUE if max_value is None else max_value

    if on_hvg_subset:
        if "highly_variable" not in adata.var:
            raise KeyError("Run select_hvgs() before regress_and_scale(on_hvg_subset=True).")
        adata = adata[:, adata.var["highly_variable"]].copy()
        print(f"Subset to {adata.n_vars} HVGs before regression/scaling.")

    if regress:
        missing = [c for c in covariates if c not in adata.obs]
        if missing:
            raise KeyError(f"Regression covariates not in adata.obs: {missing}")
        print(f"Regressing out {covariates} (this is the slow step) ...")
        sc.pp.regress_out(adata, covariates)

    sc.pp.scale(adata, max_value=max_value, zero_center=True)
    print(f"Scaled to unit variance, clipped at +/-{max_value}. "
          "adata.X is now z-scores -- use layers['lognorm'] or .raw for expression.")
    return adata


def preprocessing_summary(adata: ad.AnnData) -> pd.DataFrame:
    """One-row record of what was actually applied, for the Methods section."""
    return pd.DataFrame([{
        "n_cells": adata.n_obs,
        "n_genes_in_X": adata.n_vars,
        "n_genes_in_raw": adata.raw.n_vars if adata.raw is not None else np.nan,
        "target_sum": cfg.TARGET_SUM,
        "n_top_hvg_requested": cfg.N_TOP_HVG,
        "hvg_batch_key": cfg.HVG_BATCH_KEY,
        "regress_out": cfg.REGRESS_OUT,
        "regress_covariates": ",".join(cfg.REGRESS_COVARIATES) if cfg.REGRESS_OUT else "",
        "regress_on_hvg_subset": cfg.REGRESS_ON_HVG_SUBSET,
        "scale_max_value": cfg.SCALE_MAX_VALUE,
        "layers_present": ",".join(layer_names(adata)),
    }])


def get_expression(adata: ad.AnnData, gene: str, layer: str = "lognorm") -> np.ndarray:
    """Return a dense 1-D vector of a gene's expression from a named layer.

    Use layer='counts' for the EGFP positivity call and layer='lognorm' for any
    plotted or quoted expression value.
    """
    if gene not in adata.var_names:
        raise KeyError(f"Gene '{gene}' is not in adata.var_names.")
    if layer == "X":
        X = adata[:, gene].X
    else:
        if layer not in adata.layers:
            raise KeyError(
                f"Layer '{layer}' not found. Available: {layer_names(adata)}. "
                "If the object was subset to HVGs, use .raw or reload the "
                "preprocessed checkpoint."
            )
        X = adata[:, gene].layers[layer]
    return np.asarray(X.todense()).ravel() if hasattr(X, "todense") else np.asarray(X).ravel()
