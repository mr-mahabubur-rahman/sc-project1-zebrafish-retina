"""Marker-based cell-type annotation.

The functions here produce *evidence*. The annotation itself is a human decision
recorded in a dictionary in the notebook, because a cluster label is a scientific
claim that has to be defensible at an oral defence.

Two rules are enforced in code rather than left to discipline:
  1. A marker that is not in adata.var_names is reported as unavailable. Nothing
     is substituted for it.
  2. `build_annotation_table` requires a rationale string for every cluster, and
     any cluster left at "Unresolved" keeps that label in every figure.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from . import config as cfg
from .io_utils import check_genes_available


# --------------------------------------------------------------------------- #
# Marker discovery
# --------------------------------------------------------------------------- #

def rank_markers(
    adata: ad.AnnData,
    groupby: str | None = None,
    method: str | None = None,
    layer: str = "lognorm",
) -> ad.AnnData:
    """Rank genes per cluster with a Wilcoxon rank-sum test on log-normalised data.

    Explicitly uses the lognorm layer (via .raw when the object was subset to
    HVGs): running DE on the scaled matrix would test z-scores, not expression,
    and would silently restrict the result to the 2,000 HVGs.
    """
    groupby = cfg.LEIDEN_KEY if groupby is None else groupby
    method = cfg.DE_METHOD if method is None else method
    if groupby not in adata.obs:
        raise KeyError(f"'{groupby}' is not a column of adata.obs.")

    use_raw = adata.raw is not None and layer not in adata.layers
    sc.tl.rank_genes_groups(
        adata, groupby=groupby, method=method,
        layer=None if use_raw else layer, use_raw=use_raw,
        key_added="rank_genes_groups",
    )
    source = ".raw (lognorm, all genes)" if use_raw else f"layers['{layer}']"
    print(f"Ranked markers for '{groupby}' with {method} on {source}.")
    return adata


def markers_to_dataframe(adata: ad.AnnData, top_n: int | None = None) -> pd.DataFrame:
    """Flatten rank_genes_groups into a tidy table for tables/marker_genes.csv."""
    top_n = cfg.DE_TOP_N if top_n is None else top_n
    if "rank_genes_groups" not in adata.uns:
        raise KeyError("Run rank_markers() first.")
    df = sc.get.rank_genes_groups_df(adata, group=None)
    df = df.rename(columns={
        "group": "cluster", "names": "gene", "scores": "score",
        "logfoldchanges": "logfoldchange", "pvals": "pval", "pvals_adj": "pval_adj",
    })
    df = (
        df.sort_values(["cluster", "score"], ascending=[True, False])
        .groupby("cluster", observed=True).head(top_n).reset_index(drop=True)
    )
    return df[["cluster", "gene", "score", "logfoldchange", "pval", "pval_adj"]]


def top_markers_per_cluster(marker_df: pd.DataFrame, n: int = 10) -> pd.DataFrame:
    """Compact 'cluster -> top n genes' view for eyeballing during annotation."""
    return (
        marker_df.groupby("cluster", observed=True)
        .head(n).groupby("cluster", observed=True)["gene"]
        .apply(lambda s: ", ".join(s)).reset_index(name=f"top_{n}_genes")
    )


# --------------------------------------------------------------------------- #
# Panel scoring
# --------------------------------------------------------------------------- #

def available_panel(adata: ad.AnnData, panel: dict[str, list[str]]) -> dict[str, list[str]]:
    """Filter a marker panel down to genes actually present, reporting the losses."""
    out: dict[str, list[str]] = {}
    for cell_type, genes in panel.items():
        present, _ = check_genes_available(adata, genes, label=cell_type)
        if present:
            out[cell_type] = present
        else:
            print(f"[{cell_type}] no marker from this panel is present; the cell type "
                  "cannot be called from this panel.")
    return out


def score_panels(
    adata: ad.AnnData, panel: dict[str, list[str]], prefix: str = "score"
) -> ad.AnnData:
    """Add one sc.tl.score_genes column per cell type.

    Scores are a summary of several markers, which is exactly what we want for
    annotation -- but they are relative to a random background gene set, so a
    positive score is evidence, not proof.
    """
    usable = available_panel(adata, panel)
    for cell_type, genes in usable.items():
        key = f"{prefix}_{cell_type.replace(' ', '_').replace('/', '_')}"
        sc.tl.score_genes(adata, gene_list=genes, score_name=key, random_state=cfg.RANDOM_SEED)
    return adata


def panel_score_matrix(
    adata: ad.AnnData, cluster_key: str | None = None, prefix: str = "score"
) -> pd.DataFrame:
    """Mean panel score per cluster: the main quantitative input to annotation."""
    cluster_key = cfg.LEIDEN_KEY if cluster_key is None else cluster_key
    score_cols = [c for c in adata.obs.columns if c.startswith(f"{prefix}_")]
    if not score_cols:
        raise KeyError("No panel scores found. Run score_panels() first.")
    mat = adata.obs.groupby(cluster_key, observed=True)[score_cols].mean()
    mat.columns = [c[len(prefix) + 1:].replace("_", " ") for c in mat.columns]
    return mat.round(3)


def best_and_runner_up(score_matrix: pd.DataFrame) -> pd.DataFrame:
    """For each cluster: best-scoring panel, runner-up, and the margin between them.

    A small margin is the signal to look harder rather than to accept the top hit.
    """
    rows = []
    for cluster, row in score_matrix.iterrows():
        ordered = row.sort_values(ascending=False)
        best, second = ordered.index[0], ordered.index[1] if len(ordered) > 1 else None
        margin = float(ordered.iloc[0] - ordered.iloc[1]) if len(ordered) > 1 else np.nan
        rows.append({
            "cluster": cluster,
            "best_panel": best,
            "best_score": float(ordered.iloc[0]),
            "runner_up_panel": second,
            "runner_up_score": float(ordered.iloc[1]) if len(ordered) > 1 else np.nan,
            "margin": margin,
            "suggested_confidence": (
                "low" if (np.isnan(margin) or margin < 0.05 or ordered.iloc[0] <= 0)
                else "medium" if margin < 0.15 else "high"
            ),
        })
    return pd.DataFrame(rows)


def mean_marker_expression(
    adata: ad.AnnData, genes: list[str], cluster_key: str | None = None,
    layer: str = "lognorm",
) -> pd.DataFrame:
    """Mean log-normalised expression of specific genes per cluster."""
    cluster_key = cfg.LEIDEN_KEY if cluster_key is None else cluster_key
    present, _ = check_genes_available(adata, genes, label="mean_marker_expression")
    if not present:
        return pd.DataFrame()
    source = adata.raw.to_adata() if (adata.raw is not None and layer not in adata.layers) else adata
    X = source[:, present].layers[layer] if layer in source.layers else source[:, present].X
    dense = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
    df = pd.DataFrame(dense, columns=present, index=adata.obs_names)
    df[cluster_key] = adata.obs[cluster_key].values
    return df.groupby(cluster_key, observed=True).mean().round(3)


# --------------------------------------------------------------------------- #
# Annotation record
# --------------------------------------------------------------------------- #

def build_annotation_table(
    adata: ad.AnnData,
    annotations: dict[str, dict],
    cluster_key: str | None = None,
) -> pd.DataFrame:
    """Turn the analyst's annotation dictionary into tables/cluster_annotations.csv.

    Each entry must supply: cell_type, positive_markers, conflicting_markers,
    confidence ('high'/'medium'/'low'), rationale. Every cluster in the data must
    appear; a cluster with no defensible call is entered as 'Unresolved'.
    """
    cluster_key = cfg.LEIDEN_KEY if cluster_key is None else cluster_key
    clusters = list(adata.obs[cluster_key].cat.categories)
    missing = [c for c in clusters if c not in annotations]
    if missing:
        raise ValueError(
            f"No annotation entry for cluster(s) {missing}. Every cluster needs an "
            "entry; use cell_type='Unresolved' with a rationale where the evidence "
            "does not support a call."
        )

    required = {"cell_type", "positive_markers", "conflicting_markers", "confidence", "rationale"}
    valid_confidence = {"high", "medium", "low"}
    rows = []
    for cluster in clusters:
        entry = annotations[cluster]
        gaps = required - set(entry)
        if gaps:
            raise ValueError(f"Cluster {cluster} annotation is missing {sorted(gaps)}.")
        if entry["confidence"] not in valid_confidence:
            raise ValueError(
                f"Cluster {cluster}: confidence must be one of {sorted(valid_confidence)}."
            )
        if not str(entry["rationale"]).strip():
            raise ValueError(f"Cluster {cluster}: rationale must not be empty.")
        rows.append({
            "cluster": cluster,
            "n_cells": int((adata.obs[cluster_key] == cluster).sum()),
            "cell_type": entry["cell_type"],
            "positive_markers": ", ".join(entry["positive_markers"]),
            "conflicting_markers": ", ".join(entry["conflicting_markers"]),
            "confidence": entry["confidence"],
            "rationale": entry["rationale"],
        })
    return pd.DataFrame(rows)


def apply_annotations(
    adata: ad.AnnData, annotation_table: pd.DataFrame, cluster_key: str | None = None
) -> ad.AnnData:
    """Write the cell_type column and order its categories to match the colour map."""
    cluster_key = cfg.LEIDEN_KEY if cluster_key is None else cluster_key
    mapping = dict(zip(annotation_table["cluster"].astype(str), annotation_table["cell_type"]))
    adata.obs["cell_type"] = adata.obs[cluster_key].astype(str).map(mapping)

    if adata.obs["cell_type"].isna().any():
        raise ValueError("Some cells did not receive a cell_type; check the cluster keys.")

    observed = list(pd.unique(adata.obs["cell_type"]))
    ordered = [ct for ct in cfg.CELL_TYPE_COLORS if ct in observed]
    ordered += [ct for ct in observed if ct not in ordered]
    adata.obs["cell_type"] = pd.Categorical(adata.obs["cell_type"], categories=ordered, ordered=True)

    unknown = [ct for ct in observed if ct not in cfg.CELL_TYPE_COLORS]
    if unknown:
        print(f"NOTE: no colour defined for {unknown}. Add them to "
              "config.CELL_TYPE_COLORS so the colour stays constant across figures.")
    print(adata.obs["cell_type"].value_counts().to_string())
    return adata


# --------------------------------------------------------------------------- #
# Composition
# --------------------------------------------------------------------------- #

def cell_type_proportions(
    adata: ad.AnnData, group_key: str = "cell_type"
) -> pd.DataFrame:
    """Per-sample cell-type proportions, keeping replicate identity.

    Proportions are computed within each sample, so every replicate is one
    independent observation. Pooling cells across replicates first would weight the
    result by however many cells each run happened to capture.
    """
    counts = (
        adata.obs.groupby(["sample", "condition", "replicate", group_key], observed=True)
        .size().reset_index(name="cell_count")
    )
    totals = counts.groupby("sample", observed=True)["cell_count"].transform("sum")
    counts["proportion"] = counts["cell_count"] / totals
    counts["percent"] = (100 * counts["proportion"]).round(2)
    return counts.sort_values(["condition", "sample", group_key])


def proportions_by_condition(proportions: pd.DataFrame, group_key: str = "cell_type") -> pd.DataFrame:
    """Mean +/- range of the replicate-level proportions within each condition.

    With n = 2 replicates per condition a mean and a range are honest; a standard
    error or a p-value would not be.
    """
    grouped = proportions.groupby(["condition", group_key], observed=True)["percent"]
    df = grouped.agg(["mean", "min", "max", "count"]).reset_index()
    df = df.rename(columns={"mean": "mean_percent", "min": "min_percent",
                            "max": "max_percent", "count": "n_replicates"})
    df["mean_percent"] = df["mean_percent"].round(2)
    return df


PROPORTION_CAVEAT = (
    "Cell-type proportions from scRNA-seq measure what was captured, not what was "
    "present in the tissue. Differential survival, dissociation efficiency, cell "
    "size and adhesion, capture efficiency and sequencing depth all shift these "
    "numbers. Bise et al. make exactly this point about their own data: the higher "
    "rod fraction in MNU-treated retinas is attributed to easier release of rods "
    "from a damaged outer nuclear layer, not to more rods being present. Any "
    "proportion change reported here carries the same caveat."
)
