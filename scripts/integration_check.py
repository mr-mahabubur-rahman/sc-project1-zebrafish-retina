"""Integration sensitivity analysis.

The main pipeline does not apply batch correction, for the reason set out in
`docs/method_comparison.md`: each timepoint is its own pair of 10x runs, so
condition and sequencing batch are completely confounded, and correcting on sample
identity risks removing the injury response along with the batch effect. The
variance decomposition in step 04 supports that decision — condition explains
roughly nine-fold more principal component variance than replicate.

That is an argument, not a demonstration. This module runs the pipeline again WITH
Harmony correction, writes the result to a separate directory, and compares the
two on the quantities the report actually depends on:

  * how many clusters each produces, and how well they agree
  * the EGFP time course (expected to be unchanged, since positivity is a raw
    count call that does not depend on clustering)
  * the rod paralog table (the one place integration could plausibly change a
    conclusion)
  * cell-type proportions

Neither result replaces the other. The uncorrected analysis remains the primary
one; this is reported as a sensitivity analysis, and whichever way it comes out is
worth stating.

Requires harmonypy:  uv pip install harmonypy
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from . import config as cfg
from . import clustering, egfp_analysis
from .io_utils import check_genes_available

# Separate output directory so the primary results are never overwritten.
INT_RESULTS = cfg.REPO_ROOT / "results_integrated"
INT_TABLES = cfg.TABLES_DIR / "integration_check"


def _ensure_dirs() -> None:
    INT_RESULTS.mkdir(parents=True, exist_ok=True)
    INT_TABLES.mkdir(parents=True, exist_ok=True)


def save_int_table(df: pd.DataFrame, name: str, index: bool = False) -> Path:
    """Write a comparison table into tables/integration_check/."""
    _ensure_dirs()
    path = INT_TABLES / name
    df.to_csv(path, index=index)
    print(f"Wrote {path.relative_to(cfg.REPO_ROOT)}  ({len(df)} rows)")
    return path


# --------------------------------------------------------------------------- #
# 1. Run the corrected pipeline
# --------------------------------------------------------------------------- #

def run_integrated(adata: ad.AnnData, batch_key: str = "sample",
                   resolution: float | None = None) -> ad.AnnData:
    """Apply Harmony to the PCA embedding, then rebuild the graph and clustering.

    Harmony corrects the PCA representation rather than the expression matrix, so
    .raw and layers['lognorm'] are untouched and any differential expression run
    afterwards uses the same values as the uncorrected analysis. Only the
    neighbourhood structure changes.
    """
    resolution = cfg.LEIDEN_RESOLUTION if resolution is None else resolution

    if "X_pca" not in adata.obsm:
        adata = clustering.run_pca(adata)

    # harmonypy is called directly rather than through sc.external.pp.
    # harmony_integrate: that wrapper transposes the result, which was correct for
    # harmonypy < 2.0 (PCs x cells) but breaks on 2.0+ (cells x PCs) with an
    # obsm shape error. Orienting by shape here works on either version.
    try:
        import harmonypy
    except ImportError as exc:
        raise ImportError(
            "harmonypy is not installed. Install it deliberately with\n"
            "    uv pip install harmonypy\n"
            "and record it in requirements.txt and the Methods section."
        ) from exc

    if batch_key not in adata.obs:
        raise KeyError(f"Batch key '{batch_key}' is not a column of adata.obs.")

    ho = harmonypy.run_harmony(adata.obsm["X_pca"], adata.obs, [batch_key],
                               random_state=cfg.RANDOM_SEED)
    Z = np.asarray(ho.Z_corr)
    if Z.shape[0] == adata.n_obs:
        corrected = Z
    elif Z.shape[1] == adata.n_obs:
        corrected = Z.T
    else:
        raise ValueError(
            f"Harmony returned an array of shape {Z.shape}, which matches neither "
            f"{adata.n_obs} cells nor the PCA dimensions. Check the harmonypy "
            "version against the wrapper expectations."
        )
    adata.obsm["X_pca_harmony"] = np.ascontiguousarray(corrected)

    print(f"Harmony applied on '{batch_key}' (harmonypy returned "
          f"{Z.shape[0]}x{Z.shape[1]}); corrected embedding of shape "
          f"{adata.obsm['X_pca_harmony'].shape} in obsm['X_pca_harmony'].")

    sc.pp.neighbors(adata, n_neighbors=cfg.N_NEIGHBORS, use_rep="X_pca_harmony",
                    random_state=cfg.RANDOM_SEED)
    sc.tl.umap(adata, random_state=cfg.RANDOM_SEED)
    key = f"leiden_int_{resolution}"
    try:
        sc.tl.leiden(adata, resolution=resolution, key_added=key,
                     random_state=cfg.RANDOM_SEED, flavor="igraph",
                     n_iterations=2, directed=False)
    except TypeError:
        sc.tl.leiden(adata, resolution=resolution, key_added=key,
                     random_state=cfg.RANDOM_SEED)
    adata.uns["integration"] = {"method": "harmony", "batch_key": batch_key,
                                "resolution": resolution, "cluster_key": key}
    print(f"Leiden on the corrected embedding: {adata.obs[key].nunique()} clusters "
          f"-> obs['{key}']")
    return adata


# --------------------------------------------------------------------------- #
# 2. Compare the two clusterings
# --------------------------------------------------------------------------- #

def cluster_agreement(adata_int: ad.AnnData, uncorrected_key: str | None = None,
                      integrated_key: str | None = None) -> dict:
    """Adjusted Rand index and cluster counts for the two partitions.

    Both partitions must be present as columns of the same object, so the cells
    are matched. ARI near 1 means integration barely changed the partition; a low
    value means the structure was substantially rearranged.
    """
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

    uncorrected_key = uncorrected_key or cfg.LEIDEN_KEY
    integrated_key = integrated_key or adata_int.uns["integration"]["cluster_key"]
    for k in (uncorrected_key, integrated_key):
        if k not in adata_int.obs:
            raise KeyError(f"'{k}' is not a column of adata.obs.")

    a = adata_int.obs[uncorrected_key].astype(str)
    b = adata_int.obs[integrated_key].astype(str)
    out = {
        "n_clusters_uncorrected": int(a.nunique()),
        "n_clusters_integrated": int(b.nunique()),
        "adjusted_rand_index": round(float(adjusted_rand_score(a, b)), 4),
        "normalised_mutual_info": round(float(normalized_mutual_info_score(a, b)), 4),
        "n_cells": int(adata_int.n_obs),
    }
    print(f"Clusters: {out['n_clusters_uncorrected']} uncorrected vs "
          f"{out['n_clusters_integrated']} integrated")
    print(f"Adjusted Rand index: {out['adjusted_rand_index']}  "
          f"(1.0 = identical partition, 0 = no better than chance)")
    return out


def batch_variance_before_after(adata_int: ad.AnnData, n_pcs: int | None = None) -> pd.DataFrame:
    """Variance decomposition on the uncorrected and corrected embeddings.

    This is the quantity that motivated skipping integration in the first place;
    running it on both shows how much the correction actually removed.
    """
    n_pcs = cfg.N_PCS_NEIGHBORS if n_pcs is None else n_pcs
    rows = []
    for label, rep in [("uncorrected", "X_pca"), ("harmony", "X_pca_harmony")]:
        if rep not in adata_int.obsm:
            continue
        X = adata_int.obsm[rep][:, :n_pcs]
        for factor in ("condition", "replicate", "sample"):
            labels = adata_int.obs[factor].astype(str).values
            frac = []
            for pc in range(n_pcs):
                v = X[:, pc]
                gm = v.mean()
                ss_tot = ((v - gm) ** 2).sum()
                ss_bet = sum((labels == lev).sum() * (v[labels == lev].mean() - gm) ** 2
                             for lev in np.unique(labels))
                frac.append(ss_bet / ss_tot if ss_tot > 0 else np.nan)
            rows.append({"embedding": label, "factor": factor,
                         "mean_variance_explained": round(float(np.nanmean(frac)), 4)})
    df = pd.DataFrame(rows)
    print(df.pivot(index="factor", columns="embedding",
                   values="mean_variance_explained").to_string())
    return df


# --------------------------------------------------------------------------- #
# 3. Transfer cell type labels
# --------------------------------------------------------------------------- #

def transfer_labels(adata_int: ad.AnnData, annotated: ad.AnnData,
                    integrated_key: str | None = None,
                    label_key: str = "cell_type") -> ad.AnnData:
    """Assign each integrated cluster the majority cell type from the annotation.

    For a sensitivity analysis the question is whether integration changes the
    structure, not whether it changes the analyst's judgement, so labels are
    carried across by barcode rather than re-derived. The purity column records
    how cleanly each integrated cluster maps onto a single existing type; low
    purity means integration merged populations that were previously distinct.
    """
    integrated_key = integrated_key or adata_int.uns["integration"]["cluster_key"]
    shared = adata_int.obs_names.intersection(annotated.obs_names)
    if len(shared) == 0:
        raise ValueError("No shared barcodes; the two objects are not comparable.")
    if len(shared) < adata_int.n_obs:
        print(f"NOTE: {adata_int.n_obs - len(shared)} cells absent from the "
              "annotated object and left unlabelled.")

    src = annotated.obs.loc[shared, label_key].astype(str)
    adata_int.obs["_prev"] = pd.Series(src, index=shared).reindex(adata_int.obs_names)

    rows = []
    mapping = {}
    for cl, grp in adata_int.obs.groupby(integrated_key, observed=True):
        counts = grp["_prev"].value_counts(dropna=True)
        if counts.empty:
            mapping[cl] = "Unresolved"
            continue
        top = counts.index[0]
        mapping[cl] = top
        rows.append({"integrated_cluster": cl, "n_cells": int(len(grp)),
                     "majority_cell_type": top,
                     "purity": round(float(counts.iloc[0] / counts.sum()), 3),
                     "second_type": counts.index[1] if len(counts) > 1 else "",
                     "second_frac": round(float(counts.iloc[1] / counts.sum()), 3)
                                     if len(counts) > 1 else 0.0})
    adata_int.obs["cell_type_transferred"] = (
        adata_int.obs[integrated_key].astype(str).map(mapping))
    purity = pd.DataFrame(rows).sort_values("n_cells", ascending=False)
    low = purity[purity["purity"] < 0.8]
    if len(low):
        print(f"{len(low)} integrated cluster(s) below 0.8 purity — these merge "
              "populations that were separate before correction:")
        print(low[["integrated_cluster", "n_cells", "majority_cell_type",
                   "purity", "second_type"]].to_string(index=False))
    else:
        print("Every integrated cluster maps cleanly (purity >= 0.8) onto one "
              "previously annotated cell type.")
    return adata_int, purity


# --------------------------------------------------------------------------- #
# 4. The three comparisons the report depends on
# --------------------------------------------------------------------------- #

def compare_egfp(adata_int: ad.AnnData, annotated: ad.AnnData) -> pd.DataFrame:
    """EGFP-positive percentage per condition, both analyses.

    Expected to be identical: positivity is a raw transgene count call that does
    not depend on clustering or on the embedding. A difference here would indicate
    a bug, not a biological effect.
    """
    out = []
    for label, obj in [("uncorrected", annotated), ("integrated", adata_int)]:
        if "egfp_positive" not in obj.obs:
            obj, _ = egfp_analysis.add_egfp_columns(obj)
        per_sample = (obj.obs.groupby(["condition", "sample"], observed=True)
                      ["egfp_positive"].agg(["sum", "size"]).reset_index())
        per_sample["pct"] = 100 * per_sample["sum"] / per_sample["size"]
        for cond, grp in per_sample.groupby("condition", observed=True):
            out.append({"analysis": label, "condition": str(cond),
                        "mean_pct_positive": round(float(grp["pct"].mean()), 2),
                        "min_pct": round(float(grp["pct"].min()), 2),
                        "max_pct": round(float(grp["pct"].max()), 2)})
    df = pd.DataFrame(out)
    wide = df.pivot(index="condition", columns="analysis", values="mean_pct_positive")
    wide = wide.reindex(cfg.CONDITION_ORDER)
    wide["difference"] = (wide["integrated"] - wide["uncorrected"]).round(3)
    print("\nEGFP-positive % per condition:")
    print(wide.to_string())
    if wide["difference"].abs().max() < 1e-9:
        print("Identical, as expected — positivity does not depend on the embedding.")
    else:
        print("DIFFERENT. Positivity should not depend on clustering; investigate "
              "before interpreting.")
    return df


def compare_egfp_enrichment(adata_int: ad.AnnData, annotated: ad.AnnData) -> pd.DataFrame:
    """EGFP-positive percentage per cell type, both analyses.

    Unlike the time course, this CAN change: cells move between clusters, so the
    denominator of each cell type changes even though the positivity call does not.
    """
    rows = []
    for label, obj, key in [("uncorrected", annotated, "cell_type"),
                            ("integrated", adata_int, "cell_type_transferred")]:
        if key not in obj.obs:
            continue
        g = obj.obs.groupby(key, observed=True)["egfp_positive"]
        d = pd.DataFrame({"n_cells": g.size(), "n_pos": g.sum()}).reset_index()
        d = d.rename(columns={key: "cell_type"})
        d["pct_positive"] = (100 * d["n_pos"] / d["n_cells"]).round(2)
        d["analysis"] = label
        rows.append(d)
    df = pd.concat(rows, ignore_index=True)
    wide = df.pivot(index="cell_type", columns="analysis", values="pct_positive")
    if {"uncorrected", "integrated"}.issubset(wide.columns):
        wide["difference"] = (wide["integrated"] - wide["uncorrected"]).round(2)
        print("\nEGFP-positive % per cell type:")
        print(wide.sort_values("uncorrected", ascending=False).to_string())
    return df


PARALOGS = ["rho", "rhol", "pde6ga", "pde6gb", "guca1a", "guca1b"]


def compare_rod_paralogs(adata_int: ad.AnnData, annotated: ad.AnnData,
                         resolution: float = 0.3) -> tuple[pd.DataFrame, dict]:
    """Re-run the rod paralog test on the integrated data.

    This is the one comparison where integration could plausibly change a
    conclusion. The test is the same as in the main analysis: an inversion
    requires the pairs to split in OPPOSITE directions across subclusters.
    """
    rod_types = [c for c in annotated.obs["cell_type"].astype(str).unique()
                 if "rod" in c.lower()]
    if not rod_types:
        raise ValueError("No annotated rod population to subset.")

    mask = adata_int.obs["cell_type_transferred"].astype(str).isin(rod_types)
    rods = adata_int[mask].copy()
    print(f"\nRod subset from the integrated object: {rods.n_obs} cells "
          f"({', '.join(rod_types)})")

    sub = egfp_analysis.recluster_subset(rods, resolution=resolution)
    source = sub.raw.to_adata() if sub.raw is not None else sub
    present, missing = check_genes_available(source, PARALOGS, "rod paralogs (integrated)")
    if missing:
        print(f"  absent from this annotation: {missing}")

    X = source[:, present].X
    dense = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
    df = pd.DataFrame(dense, columns=present, index=sub.obs_names)
    df["subcluster"] = sub.obs["subcluster"].values
    means = df.groupby("subcluster", observed=True).mean().round(3)
    means.insert(0, "n_cells", sub.obs.groupby("subcluster", observed=True).size())

    print("\nParalog expression across integrated rod subclusters:")
    print(means.to_string())

    verdict = {}
    for a, b in [("rho", "rhol"), ("pde6ga", "pde6gb"), ("guca1a", "guca1b")]:
        if a in means and b in means:
            higher = means[a] > means[b]
            if higher.all():
                verdict[f"{a}/{b}"] = f"{a} > {b} in every subcluster — no inversion"
            elif (~higher).all():
                verdict[f"{a}/{b}"] = f"{b} > {a} in every subcluster — no inversion"
            else:
                verdict[f"{a}/{b}"] = "SPLIT across subclusters — inversion present"
    print("\nInversion test on the integrated data:")
    for k, v in verdict.items():
        print(f"  {k}: {v}")
    n_split = sum("inversion present" in v for v in verdict.values())
    print(f"\n{n_split} of {len(verdict)} pairs split. The published result requires "
          "all three to split in a coordinated, reciprocal way.")
    if n_split >= 2:
        print("This differs from the uncorrected analysis (1 of 3) and is worth "
              "reporting: it would locate the failure in the batch structure rather "
              "than in normalisation.")
    else:
        print("Consistent with the uncorrected analysis: integration does not "
              "recover the published rod split.")
    return means.reset_index(), verdict


def compare_proportions(adata_int: ad.AnnData, annotated: ad.AnnData) -> pd.DataFrame:
    """Cell-type proportions per condition, both analyses."""
    rows = []
    for label, obj, key in [("uncorrected", annotated, "cell_type"),
                            ("integrated", adata_int, "cell_type_transferred")]:
        pct = (pd.crosstab(obj.obs[key], obj.obs["condition"], normalize="columns") * 100)
        pct = pct.reindex(columns=cfg.CONDITION_ORDER).round(2)
        pct = pct.reset_index().rename(columns={key: "cell_type"})
        pct["analysis"] = label
        rows.append(pct)
    df = pd.concat(rows, ignore_index=True)
    print("\nCell-type proportions (%) per condition, both analyses — "
          "see tables/integration_check/proportions_comparison.csv")
    return df


# --------------------------------------------------------------------------- #
# 5. Summary
# --------------------------------------------------------------------------- #

def summarise(agreement: dict, egfp: pd.DataFrame, verdict: dict,
              purity: pd.DataFrame) -> pd.DataFrame:
    """One table stating what integration changed and what it did not."""
    wide = egfp.pivot(index="condition", columns="analysis", values="mean_pct_positive")
    egfp_max_diff = float((wide["integrated"] - wide["uncorrected"]).abs().max())
    n_split = sum("inversion present" in v for v in verdict.values())
    impure = int((purity["purity"] < 0.8).sum())

    rows = [
        {"quantity": "Number of clusters",
         "uncorrected": agreement["n_clusters_uncorrected"],
         "integrated": agreement["n_clusters_integrated"],
         "changed": agreement["n_clusters_uncorrected"] != agreement["n_clusters_integrated"]},
        {"quantity": "Adjusted Rand index between partitions",
         "uncorrected": "—", "integrated": agreement["adjusted_rand_index"],
         "changed": agreement["adjusted_rand_index"] < 0.9},
        {"quantity": "Integrated clusters mixing >1 cell type (purity <0.8)",
         "uncorrected": "—", "integrated": impure, "changed": impure > 0},
        {"quantity": "EGFP time course, max difference (pp)",
         "uncorrected": "—", "integrated": round(egfp_max_diff, 3),
         "changed": egfp_max_diff > 0.01},
        {"quantity": "Rod paralog pairs splitting (of 3)",
         "uncorrected": 1, "integrated": n_split, "changed": n_split != 1},
    ]
    df = pd.DataFrame(rows)
    print("\n" + "=" * 66)
    print("INTEGRATION SENSITIVITY — SUMMARY")
    print("=" * 66)
    print(df.to_string(index=False))
    print("\nRead this as: which conclusions in the report would change if the "
          "batch correction decision had gone the other way. Anything with "
          "changed=False is robust to that choice.")
    return df


# --------------------------------------------------------------------------- #
# 6. Figures
# --------------------------------------------------------------------------- #

import matplotlib.pyplot as plt  # noqa: E402

INT_FIGURES = cfg.FIGURES_DIR / "integration_check"


def _save_fig(fig, name: str):
    INT_FIGURES.mkdir(parents=True, exist_ok=True)
    for ext in cfg.FIGURE_FORMATS:
        fig.savefig(INT_FIGURES / f"{name}.{ext}", format=ext,
                    dpi=cfg.FIGURE_DPI, bbox_inches="tight")
    print(f"Saved {name} -> figures/integration_check/")
    return fig


def figure_embeddings(adata_int: ad.AnnData, annotated: ad.AnnData,
                      name: str = "int_fig1_umap_before_after") -> plt.Figure:
    """Both embeddings coloured by sample and by condition.

    The top row shows what the correction was asked to do — mix the samples. The
    bottom row shows what it also did: flatten the separation between conditions,
    which is the biological signal. Read the two rows together.
    """
    from .plotting import set_style
    set_style()

    unc = annotated.obsm["X_umap"]
    cor = adata_int.obsm["X_umap"]
    order = [c for c in cfg.EXPECTED_SAMPLES
             if c in set(annotated.obs["sample"].astype(str))]
    sample_colors = dict(zip(order, plt.get_cmap("tab10").colors))

    fig, axes = plt.subplots(2, 2, figsize=(12, 11))
    for col, (coords, obs, label) in enumerate([
            (unc, annotated.obs, "Uncorrected"),
            (cor, adata_int.obs, "Harmony-corrected")]):
        ax = axes[0, col]
        for s in order:
            m = (obs["sample"].astype(str) == s).values
            ax.scatter(coords[m, 0], coords[m, 1], s=1.2, alpha=0.55,
                       color=sample_colors[s], label=s, rasterized=True)
        ax.set_title(f"{label} — by sample", fontsize=11)
        if col == 1:
            ax.legend(markerscale=8, fontsize=7, frameon=False,
                      loc="center left", bbox_to_anchor=(1.01, 0.5))

        ax = axes[1, col]
        for c in cfg.CONDITION_ORDER:
            m = (obs["condition"].astype(str) == c).values
            ax.scatter(coords[m, 0], coords[m, 1], s=1.2, alpha=0.55,
                       color=cfg.CONDITION_COLORS[c],
                       label=cfg.CONDITION_LABELS[c], rasterized=True)
        ax.set_title(f"{label} — by condition", fontsize=11)
        if col == 1:
            ax.legend(markerscale=8, fontsize=8, frameon=False,
                      loc="center left", bbox_to_anchor=(1.01, 0.5))

    for ax in axes.ravel():
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlabel("UMAP1", fontsize=8); ax.set_ylabel("UMAP2", fontsize=8)
    fig.suptitle("Integration sensitivity | embeddings before and after correction",
                 fontsize=13)
    fig.tight_layout()
    return _save_fig(fig, name)


def figure_variance(var_df: pd.DataFrame,
                    name: str = "int_fig2_variance_before_after") -> plt.Figure:
    """Variance explained by each factor, before and after correction.

    The point of the panel: replicate is the purely technical term and is what a
    batch correction should remove. If condition falls further than replicate, the
    correction is removing biology.
    """
    from .plotting import set_style
    set_style()
    wide = var_df.pivot(index="factor", columns="embedding",
                        values="mean_variance_explained")
    wide = wide.reindex(["condition", "replicate", "sample"])

    x = np.arange(len(wide))
    w = 0.36
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))

    ax = axes[0]
    ax.bar(x - w/2, wide["uncorrected"], w, label="Uncorrected", color="#4D4D4D")
    ax.bar(x + w/2, wide["harmony"], w, label="Harmony-corrected", color="#B2182B")
    ax.set_xticks(x); ax.set_xticklabels(wide.index)
    ax.set_ylabel("mean variance explained\n(first 20 PCs)")
    ax.set_title("A  Variance by factor")
    ax.legend(frameon=False)
    for i, f in enumerate(wide.index):
        ax.text(i + w/2, wide.loc[f, "harmony"], f"{wide.loc[f,'harmony']:.4f}",
                ha="center", va="bottom", fontsize=7)
        ax.text(i - w/2, wide.loc[f, "uncorrected"], f"{wide.loc[f,'uncorrected']:.4f}",
                ha="center", va="bottom", fontsize=7)

    ax = axes[1]
    pct = 100 * (1 - wide["harmony"] / wide["uncorrected"])
    colors = ["#B2182B" if f == "condition" else "#4D4D4D" for f in wide.index]
    ax.bar(wide.index, pct, color=colors)
    for i, f in enumerate(wide.index):
        ax.text(i, pct[f], f"{pct[f]:.1f}%", ha="center", va="bottom", fontsize=9)
    ax.set_ylabel("% of variance removed")
    ax.set_title("B  Proportion removed by correction")
    ax.axhline(0, color="black", lw=0.8)
    ax.set_ylim(0, max(105, float(pct.max()) * 1.18))

    fig.suptitle("Integration sensitivity | what the correction removed", fontsize=13)
    fig.tight_layout()
    # Placed as a figure-level note rather than inside the axes, where it would
    # sit on top of the bars.
    fig.text(0.5, -0.03,
             "A batch correction should remove replicate variance (technical), "
             "not condition variance (biological). If condition falls further "
             "than replicate, the correction is removing the signal.",
             ha="center", va="top", fontsize=8.5, style="italic", color="#555555")
    return _save_fig(fig, name)


def figure_purity(adata_int: ad.AnnData, integrated_key: str | None = None,
                  label_key: str = "_prev",
                  name: str = "int_fig3_cluster_confusion") -> plt.Figure:
    """Where each annotated cell type ends up after correction.

    Rows are the original cell types, columns the integrated clusters, values the
    row-normalised percentage. A cell type confined to one column survived intact;
    one spread across several was dispersed by the correction.
    """
    from .plotting import set_style
    set_style()
    integrated_key = integrated_key or adata_int.uns["integration"]["cluster_key"]
    if label_key not in adata_int.obs:
        raise KeyError("Run transfer_labels() first (it creates obs['_prev']).")

    ct = pd.crosstab(adata_int.obs[label_key], adata_int.obs[integrated_key])
    ct = ct.loc[ct.sum(axis=1).sort_values(ascending=False).index]
    pct = 100 * ct.div(ct.sum(axis=1), axis=0)

    fig, ax = plt.subplots(figsize=(0.44 * pct.shape[1] + 5, 0.42 * pct.shape[0] + 2.6))
    im = ax.imshow(pct.values, cmap="Reds", aspect="auto", vmin=0, vmax=100)
    ax.set_xticks(range(pct.shape[1]))
    ax.set_xticklabels(pct.columns, fontsize=8)
    ax.set_yticks(range(pct.shape[0]))
    ax.set_yticklabels([f"{i}  (n={ct.loc[i].sum():,})" for i in pct.index], fontsize=8)
    ax.set_xlabel("integrated cluster")
    ax.set_title("Integration sensitivity | where each annotated cell type went\n"
                 "(row-normalised %; a row confined to one column survived intact)",
                 fontsize=11)
    for i in range(pct.shape[0]):
        for j in range(pct.shape[1]):
            v = pct.values[i, j]
            if v >= 10:
                ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=6.5,
                        color="white" if v > 55 else "black")
    fig.colorbar(im, ax=ax, label="% of the cell type", shrink=0.75)
    fig.tight_layout()
    return _save_fig(fig, name)


def figure_key_comparisons(egfp_df: pd.DataFrame, paralogs: pd.DataFrame,
                           uncorrected_paralogs: pd.DataFrame | None = None,
                           name: str = "int_fig4_key_results") -> plt.Figure:
    """The two conclusions that must not depend on the correction.

    Left: the reporter time course under both analyses — expected to be identical.
    Right: the rod paralog ratios after correction; an inversion would require the
    bars to cross the 1.0 line in opposite directions for different pairs.
    """
    from .plotting import set_style
    set_style()
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6))

    ax = axes[0]
    wide = egfp_df.pivot(index="condition", columns="analysis",
                         values="mean_pct_positive").reindex(cfg.CONDITION_ORDER)
    x = np.arange(len(wide))
    ax.bar(x - 0.18, wide["uncorrected"], 0.36, label="Uncorrected", color="#4D4D4D")
    ax.bar(x + 0.18, wide["integrated"], 0.36, label="Harmony-corrected", color="#B2182B")
    for i, c in enumerate(wide.index):
        ax.text(i, max(wide.loc[c]) + 0.15, f"{wide.loc[c,'uncorrected']:.2f}",
                ha="center", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels([cfg.CONDITION_LABELS[c] for c in wide.index])
    ax.set_ylabel("% EGFP-positive cells")
    ax.set_title("A  Reporter time course — identical under both")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1]
    pairs = [("rho", "rhol"), ("pde6ga", "pde6gb"), ("guca1a", "guca1b")]
    sub = paralogs.set_index("subcluster") if "subcluster" in paralogs else paralogs
    width = 0.8 / max(len(sub), 1)
    for i, (a, b) in enumerate(pairs):
        if a not in sub or b not in sub:
            continue
        ratios = (sub[a] / sub[b]).values
        for j, r in enumerate(ratios):
            ax.bar(i + (j - len(ratios)/2 + 0.5) * width, r, width * 0.92,
                   color="#B2182B" if r < 1 else "#4D4D4D")
    ax.axhline(1.0, color="black", ls="--", lw=1)
    ax.set_xticks(range(len(pairs)))
    ax.set_xticklabels([f"{a}/{b}" for a, b in pairs], style="italic")
    ax.set_ylabel("expression ratio")
    ax.set_title("B  Rod paralog ratios after correction\n"
                 "(bars per subcluster; an inversion needs pairs to cross 1.0 "
                 "in opposite directions)", fontsize=10)
    ax.text(0.02, 0.95, "red = ratio below 1", transform=ax.transAxes,
            fontsize=8, va="top", color="#B2182B")

    fig.suptitle("Integration sensitivity | conclusions that do not depend on the "
                 "correction", fontsize=13)
    fig.tight_layout()
    return _save_fig(fig, name)
