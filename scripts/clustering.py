"""PCA, batch diagnostics, neighbour graph, UMAP and Leiden clustering.

Integration is not applied by default. The paper integrated across samples, but
the biological signal we are after -- the injury time course -- is confounded with
sample identity by design (each timepoint is its own pair of 10x runs). Correcting
on `sample` can therefore erase the effect the project is meant to measure. The
policy here: quantify the batch structure first, then decide, then document.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from . import config as cfg


def run_pca(adata: ad.AnnData, n_comps: int = 50, seed: int | None = None) -> ad.AnnData:
    """PCA on the scaled HVG matrix."""
    seed = cfg.RANDOM_SEED if seed is None else seed
    n_comps = int(min(n_comps, adata.n_vars - 1, adata.n_obs - 1))
    sc.tl.pca(adata, n_comps=n_comps, svd_solver="arpack", random_state=seed)
    var_ratio = adata.uns["pca"]["variance_ratio"]
    print(f"PCA: {n_comps} components, {100 * var_ratio.sum():.1f}% of variance retained.")
    return adata


def pca_variance_table(adata: ad.AnnData) -> pd.DataFrame:
    """Per-PC variance explained plus cumulative, to justify the n_pcs choice."""
    var_ratio = adata.uns["pca"]["variance_ratio"]
    return pd.DataFrame({
        "PC": np.arange(1, len(var_ratio) + 1),
        "variance_ratio": var_ratio,
        "cumulative_variance_ratio": np.cumsum(var_ratio),
    })


def suggest_n_pcs(adata: ad.AnnData, drop_threshold: float = 0.001) -> int:
    """Heuristic elbow: first PC whose marginal variance ratio falls below a floor.

    Advisory only. The baseline stays at cfg.N_PCS_NEIGHBORS unless the analyst
    changes it explicitly after looking at the plot.
    """
    var_ratio = adata.uns["pca"]["variance_ratio"]
    below = np.where(var_ratio < drop_threshold)[0]
    suggestion = int(below[0]) if below.size else len(var_ratio)
    print(f"Elbow heuristic suggests ~{suggestion} PCs; baseline in use is "
          f"{cfg.N_PCS_NEIGHBORS} (config.N_PCS_NEIGHBORS).")
    return suggestion


def batch_diagnostics(adata: ad.AnnData, n_pcs: int | None = None) -> pd.DataFrame:
    """Compare replicate-driven and condition-driven structure in PCA space.

    For each of the first `n_pcs` components, computes the fraction of variance
    explained by `condition` and by `replicate nested in condition` via a simple
    one-way decomposition. If replicate explains as much as condition, integration
    is worth considering; if condition dominates, correcting on sample would remove
    biology.
    """
    n_pcs = cfg.N_PCS_NEIGHBORS if n_pcs is None else n_pcs
    X = adata.obsm["X_pca"][:, :n_pcs]
    rows = []
    for factor in ("condition", "replicate", "sample"):
        labels = adata.obs[factor].astype(str).values
        for pc in range(n_pcs):
            values = X[:, pc]
            grand_mean = values.mean()
            ss_total = ((values - grand_mean) ** 2).sum()
            ss_between = sum(
                (labels == level).sum() * (values[labels == level].mean() - grand_mean) ** 2
                for level in np.unique(labels)
            )
            rows.append({
                "factor": factor,
                "PC": pc + 1,
                "variance_explained": ss_between / ss_total if ss_total > 0 else np.nan,
            })
    df = pd.DataFrame(rows)
    summary = df.groupby("factor")["variance_explained"].mean().round(4)
    print("Mean fraction of PC variance explained (first "
          f"{n_pcs} PCs):\n{summary.to_string()}")
    return df


def interpret_batch_diagnostics(diag: pd.DataFrame) -> list[str]:
    """Turn the diagnostic numbers into an explicit integrate / do-not-integrate call."""
    means = diag.groupby("factor")["variance_explained"].mean()
    cond, rep = means.get("condition", np.nan), means.get("replicate", np.nan)
    notes = [
        f"condition explains {cond:.1%} of PC variance on average; "
        f"replicate explains {rep:.1%}."
    ]
    if rep > cond:
        notes.append(
            "Replicate structure exceeds condition structure. Integration is "
            "justifiable -- if applied, set config.USE_INTEGRATION = True, record "
            "the method, and show the UMAP before and after."
        )
    else:
        notes.append(
            "Condition structure dominates, which is the expected design here. "
            "Integrating on 'sample' would remove the injury effect along with the "
            "batch effect, so the uncorrected embedding is kept as the primary "
            "result."
        )
    notes.append(
        "Caveat: condition and 10x run are fully confounded in this design "
        "(each timepoint has its own pair of runs), so no method can fully separate "
        "true injury response from run-to-run technical variation. State this as a "
        "limitation rather than claiming the batch effect was removed."
    )
    return notes


def run_integration(adata: ad.AnnData, method: str | None = None,
                    batch_key: str | None = None) -> str:
    """Apply batch correction and return the name of the representation to use.

    Only call this after the diagnostics justify it. Harmony is used because it
    corrects the PCA embedding rather than the expression matrix, leaving
    layers['lognorm'] untouched for downstream DE.
    """
    method = cfg.INTEGRATION_METHOD if method is None else method
    batch_key = cfg.INTEGRATION_BATCH_KEY if batch_key is None else batch_key

    if method == "harmony":
        try:
            sc.external.pp.harmony_integrate(adata, key=batch_key)
        except ImportError as exc:
            raise ImportError(
                "Harmony is not installed. Add it deliberately "
                "(`uv pip install harmonypy`) and record it in requirements.txt and "
                "the Methods section, or leave USE_INTEGRATION = False."
            ) from exc
        return "X_pca_harmony"

    raise ValueError(f"Unsupported integration method '{method}'.")


def run_neighbors_umap(
    adata: ad.AnnData,
    n_neighbors: int | None = None,
    n_pcs: int | None = None,
    use_rep: str | None = None,
    seed: int | None = None,
) -> ad.AnnData:
    """Build the kNN graph and the UMAP embedding."""
    n_neighbors = cfg.N_NEIGHBORS if n_neighbors is None else n_neighbors
    n_pcs = cfg.N_PCS_NEIGHBORS if n_pcs is None else n_pcs
    seed = cfg.RANDOM_SEED if seed is None else seed

    kwargs = {"n_neighbors": n_neighbors, "random_state": seed}
    if use_rep is not None:
        kwargs["use_rep"] = use_rep
    else:
        kwargs["n_pcs"] = n_pcs

    sc.pp.neighbors(adata, **kwargs)
    sc.tl.umap(adata, random_state=seed)
    print(f"Neighbour graph: k={n_neighbors}, "
          f"{'use_rep=' + use_rep if use_rep else f'n_pcs={n_pcs}'}; UMAP computed.")
    return adata


def run_leiden(
    adata: ad.AnnData,
    resolution: float | None = None,
    key_added: str | None = None,
    seed: int | None = None,
) -> ad.AnnData:
    """Leiden clustering at one resolution."""
    resolution = cfg.LEIDEN_RESOLUTION if resolution is None else resolution
    key_added = cfg.LEIDEN_KEY if key_added is None else key_added
    seed = cfg.RANDOM_SEED if seed is None else seed
    try:
        sc.tl.leiden(
            adata, resolution=resolution, key_added=key_added, random_state=seed,
            flavor="igraph", n_iterations=2, directed=False,
        )
    except TypeError:
        # Scanpy < 1.10 has no igraph flavour; fall back to the leidenalg backend.
        sc.tl.leiden(adata, resolution=resolution, key_added=key_added, random_state=seed)
    n_clusters = adata.obs[key_added].nunique()
    print(f"Leiden resolution {resolution}: {n_clusters} clusters -> obs['{key_added}']")
    return adata


def resolution_sweep(adata: ad.AnnData, resolutions: list[float] | None = None) -> pd.DataFrame:
    """Cluster at several resolutions to show the reported one is not a lucky pick.

    The reported clustering remains cfg.LEIDEN_RESOLUTION; the sweep is evidence
    for the Methods section, not a search for the prettiest answer.
    """
    resolutions = cfg.LEIDEN_RESOLUTION_SWEEP if resolutions is None else resolutions
    rows = []
    for res in resolutions:
        key = f"leiden_res_{res}"
        if key not in adata.obs:
            run_leiden(adata, resolution=res, key_added=key)
        counts = adata.obs[key].value_counts()
        rows.append({
            "resolution": res,
            "n_clusters": int(counts.size),
            "smallest_cluster": int(counts.min()),
            "largest_cluster": int(counts.max()),
            "n_clusters_under_50_cells": int((counts < 50).sum()),
        })
    return pd.DataFrame(rows)


def cluster_composition(adata: ad.AnnData, cluster_key: str | None = None) -> pd.DataFrame:
    """Cells per cluster broken down by sample and condition."""
    cluster_key = cfg.LEIDEN_KEY if cluster_key is None else cluster_key
    df = (
        adata.obs.groupby([cluster_key, "condition", "sample"], observed=True)
        .size().reset_index(name="n_cells")
    )
    totals = df.groupby(cluster_key, observed=True)["n_cells"].transform("sum")
    df["pct_of_cluster"] = (100 * df["n_cells"] / totals).round(2)
    return df


def flag_single_sample_clusters(
    adata: ad.AnnData, cluster_key: str | None = None, threshold: float = 90.0
) -> pd.DataFrame:
    """Find clusters dominated by one sample -- candidate technical artefacts."""
    cluster_key = cfg.LEIDEN_KEY if cluster_key is None else cluster_key
    comp = cluster_composition(adata, cluster_key)
    top = comp.loc[comp.groupby(cluster_key, observed=True)["pct_of_cluster"].idxmax()]
    flagged = top[top["pct_of_cluster"] >= threshold].copy()
    if len(flagged):
        print(
            f"{len(flagged)} cluster(s) are >={threshold}% one sample: "
            f"{flagged[cluster_key].tolist()}. Check whether these are a real "
            "condition-specific population or a single-run artefact before "
            "annotating them."
        )
    else:
        print(f"No cluster is >={threshold}% derived from a single sample.")
    return flagged
