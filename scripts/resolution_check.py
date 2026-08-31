"""Resolution sensitivity: does the rod split appear at the published granularity?

Bise et al. clustered at resolution 0.2 with Seurat's SNN algorithm and obtained
17 clusters, among them two rod populations distinguished by reciprocal use of
three paralog pairs. The primary analysis here clustered at resolution 0.6 and
tested those pairs by subclustering *within* an already-annotated rod population —
a different procedure, and one of the four candidate explanations in the report for
why the split did not reproduce.

This module tests that explanation directly. It re-clusters the whole dataset at
resolution 0.2, identifies the rod clusters from marker expression rather than by
carrying labels across, and applies the same paralog test. If two rod clusters
emerge with reciprocal paralog usage, clustering granularity explains the failure.
If they do not, granularity is eliminated and normalisation remains the leading
candidate.

Written as a sensitivity analysis: the resolution-0.6 result stays primary.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from . import config as cfg
from . import clustering
from .io_utils import check_genes_available

RES_TABLES = cfg.TABLES_DIR / "resolution_check"
RES_FIGURES = cfg.FIGURES_DIR / "resolution_check"

PARALOGS = [("rho", "rhol"), ("pde6ga", "pde6gb"), ("guca1a", "guca1b")]
ROD_MARKERS = ["rho", "rhol", "nr2e3", "gnat1", "sagb", "pde6b"]


def _ensure_dirs() -> None:
    RES_TABLES.mkdir(parents=True, exist_ok=True)
    RES_FIGURES.mkdir(parents=True, exist_ok=True)


def save_res_table(df: pd.DataFrame, name: str, index: bool = False) -> Path:
    _ensure_dirs()
    path = RES_TABLES / name
    df.to_csv(path, index=index)
    print(f"Wrote {path.relative_to(cfg.REPO_ROOT)}  ({len(df)} rows)")
    return path


def _dense(adata: ad.AnnData, genes: list[str]) -> pd.DataFrame:
    """Dense expression frame from .raw, which holds the full gene set."""
    source = adata.raw.to_adata() if adata.raw is not None else adata
    present = [g for g in genes if g in source.var_names]
    if not present:
        return pd.DataFrame(index=adata.obs_names)
    X = source[:, present].X
    dense = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
    return pd.DataFrame(dense, columns=present, index=adata.obs_names)


# --------------------------------------------------------------------------- #
# 1. Cluster at the published resolution
# --------------------------------------------------------------------------- #

def cluster_at_resolution(adata: ad.AnnData, resolution: float = 0.2) -> ad.AnnData:
    """Re-cluster the whole dataset at a given Leiden resolution.

    The neighbour graph is reused if present, so this changes only the community
    detection step — the embedding and feature selection are untouched.
    """
    key = f"leiden_res_{resolution}"
    if key in adata.obs:
        print(f"'{key}' already present; reusing it.")
    else:
        if "neighbors" not in adata.uns:
            adata = clustering.run_neighbors_umap(adata)
        adata = clustering.run_leiden(adata, resolution=resolution, key_added=key)
    counts = adata.obs[key].value_counts().sort_index()
    print(f"Resolution {resolution}: {counts.size} clusters "
          f"(published analysis reported 17 at this resolution)")
    return adata


# nr2e3 is a rod-specific transcription factor. It is expressed at far lower
# absolute levels than rho, which makes it much less prone to the ambient
# contamination that inflates rho in every droplet (see report §2.6), and
# therefore a better discriminator of genuine rod identity.
ROD_SPECIFIC = "nr2e3"


def identify_rod_clusters(adata: ad.AnnData, cluster_key: str,
                          min_rho: float | None = None,
                          rho_quantile: float = 0.60,
                          require_specific: bool = True
                          ) -> tuple[list[str], pd.DataFrame]:
    """Find rod clusters from marker expression rather than by carrying labels.

    Identifying rods independently at this resolution matters: carrying the
    resolution-0.6 annotation across would impose the earlier partition on the new
    one, and the question is precisely whether the new partition differs.

    An absolute rho threshold does NOT work in this dataset. Rods dominate the
    tissue and lyse readily during dissociation, so ambient rho transcripts
    contaminate every droplet and every cluster carries a substantial mean —
    the same effect that made gene-set panel scores unusable in the main analysis.
    Selection is therefore relative (rho above a quantile of the cluster means)
    and additionally requires the rod-specific transcription factor nr2e3 to be
    above its own cluster-mean median.

    Pass `min_rho` to override with an absolute threshold, having inspected the
    printed table; record any such choice in the Methods.
    """
    df = _dense(adata, ROD_MARKERS)
    if df.empty:
        raise ValueError("No rod markers available in .raw.")
    df[cluster_key] = adata.obs[cluster_key].values
    means = df.groupby(cluster_key, observed=True).mean().round(3)
    means.insert(0, "n_cells", adata.obs.groupby(cluster_key, observed=True).size())

    if "rho" not in means:
        raise ValueError("rho is not available; rod clusters cannot be identified.")

    print("\nRod marker expression per cluster at this resolution:")
    print(means.to_string())

    if min_rho is not None:
        keep = means["rho"] >= min_rho
        criterion = f"mean rho >= {min_rho} (absolute threshold, supplied)"
    else:
        cut = float(means["rho"].quantile(rho_quantile))
        keep = means["rho"] >= cut
        criterion = (f"mean rho >= {cut:.3f} "
                     f"({rho_quantile:.0%} quantile of the cluster means)")
        spread = means["rho"].max() - means["rho"].min()
        print(f"\nrho ranges {means['rho'].min():.2f}-{means['rho'].max():.2f} "
              f"across clusters (spread {spread:.2f}). Every cluster carries "
              "ambient rho, so selection is relative rather than absolute.")

    if require_specific and ROD_SPECIFIC in means:
        med = float(means[ROD_SPECIFIC].median())
        keep = keep & (means[ROD_SPECIFIC] >= med)
        criterion += f" AND {ROD_SPECIFIC} >= {med:.3f} (median of cluster means)"
    elif require_specific:
        print(f"NOTE: {ROD_SPECIFIC} unavailable; selection rests on rho alone and "
              "is therefore more exposed to ambient contamination.")

    rod_clusters = means.index[keep].astype(str).tolist()
    print(f"\nCriterion: {criterion}")
    print(f"Rod clusters: {rod_clusters}  ({len(rod_clusters)} of {len(means)})")

    if len(rod_clusters) == len(means):
        print("WARNING: every cluster was selected, so the criterion has not "
              "discriminated. Inspect the table above and set min_rho explicitly.")
    elif len(rod_clusters) < 2:
        print("Fewer than two rod clusters at this resolution, so the published "
              "two-population structure does not arise here.")
    return rod_clusters, means.reset_index()


# --------------------------------------------------------------------------- #
# 2. The paralog test
# --------------------------------------------------------------------------- #

def paralog_test(adata: ad.AnnData, cluster_key: str,
                 rod_clusters: list[str]) -> tuple[pd.DataFrame, dict]:
    """Apply the published paralog test across the rod clusters at this resolution.

    The published result requires the three pairs to split in a coordinated,
    reciprocal way: one cluster high in rho, pde6ga and guca1a, another high in
    rhol, pde6gb and guca1b.
    """
    genes = [g for pair in PARALOGS for g in pair]
    source = adata.raw.to_adata() if adata.raw is not None else adata
    present, missing = check_genes_available(source, genes, "paralogs at this resolution")
    if missing:
        print(f"  absent from this annotation: {missing}")

    sub = adata[adata.obs[cluster_key].astype(str).isin(rod_clusters)]
    df = _dense(sub, present)
    df[cluster_key] = sub.obs[cluster_key].astype(str).values
    means = df.groupby(cluster_key, observed=True).mean().round(3)
    means.insert(0, "n_cells", sub.obs.groupby(cluster_key, observed=True).size())

    print(f"\nParalog expression across the {len(rod_clusters)} rod cluster(s):")
    print(means.to_string())

    verdict = {}
    for a, b in PARALOGS:
        if a in means and b in means and len(means) > 1:
            higher = means[a] > means[b]
            if higher.all():
                verdict[f"{a}/{b}"] = f"{a} > {b} in every cluster — no inversion"
            elif (~higher).all():
                verdict[f"{a}/{b}"] = f"{b} > {a} in every cluster — no inversion"
            else:
                verdict[f"{a}/{b}"] = "SPLIT across clusters — inversion present"
        elif len(means) <= 1:
            verdict[f"{a}/{b}"] = "only one rod cluster — the test cannot be applied"

    print("\nInversion test:")
    for k, v in verdict.items():
        print(f"  {k}: {v}")

    n_split = sum("inversion present" in v for v in verdict.values())
    print(f"\n{n_split} of {len(verdict)} pairs split.")
    if n_split == 3:
        print("All three pairs split. If they split in a COORDINATED way — the same "
              "cluster high in rho, pde6ga and guca1a throughout — the published "
              "result is reproduced at this resolution, and clustering granularity "
              "explains the earlier failure. Check the direction column below.")
    elif n_split == 0:
        print("No pair splits. Clustering granularity does not explain the failure; "
              "normalisation remains the leading candidate (report §4.2).")
    else:
        print("A partial, discordant split — the same outcome as the primary "
              "analysis. Clustering granularity does not explain the failure.")

    # Coordination check: does one cluster carry rho, pde6ga AND guca1a together?
    if len(means) > 1:
        coord = []
        for cl in means.index:
            row = means.loc[cl]
            pattern = []
            for a, b in PARALOGS:
                if a in row and b in row:
                    pattern.append("A" if row[a] > row[b] else "B")
            coord.append({"cluster": cl, "n_cells": int(row["n_cells"]),
                          "pattern": "".join(pattern)})
        cdf = pd.DataFrame(coord)
        print("\nPer-cluster pattern (A = rho/pde6ga/guca1a side, B = the paralog side):")
        print(cdf.to_string(index=False))
        print("The published structure would show one cluster reading AAA and "
              "another BBB.")
        if {"AAA", "BBB"}.issubset(set(cdf["pattern"])):
            print(">>> Both AAA and BBB present: the published pattern IS reproduced "
                  "at this resolution.")
        else:
            print(">>> Neither a clean AAA/BBB pair is present: the published "
                  "pattern is not reproduced at this resolution.")
        verdict["_patterns"] = dict(zip(cdf["cluster"], cdf["pattern"]))

    return means.reset_index(), verdict


# --------------------------------------------------------------------------- #
# 3. Figure
# --------------------------------------------------------------------------- #

def figure_resolution_comparison(means_02: pd.DataFrame,
                                 means_06: pd.DataFrame | None = None,
                                 name: str = "res_fig1_paralogs_by_resolution"):
    """Paralog ratios at resolution 0.2, optionally beside the primary result."""
    import matplotlib.pyplot as plt
    from .plotting import set_style
    set_style()
    _ensure_dirs()

    panels = [("Resolution 0.2 (published granularity)", means_02)]
    if means_06 is not None:
        panels.append(("Resolution 0.6 (primary analysis)", means_06))

    fig, axes = plt.subplots(1, len(panels), figsize=(6.5 * len(panels), 4.6),
                             squeeze=False)
    for ax, (title, means) in zip(axes[0], panels):
        idx = means.columns[0]
        m = means.set_index(idx)
        width = 0.8 / max(len(m), 1)
        for i, (a, b) in enumerate(PARALOGS):
            if a not in m or b not in m:
                continue
            ratios = (m[a] / m[b]).values
            for j, r in enumerate(ratios):
                ax.bar(i + (j - len(ratios) / 2 + 0.5) * width, r, width * 0.9,
                       color="#B2182B" if r < 1 else "#4D4D4D")
        ax.axhline(1.0, color="black", ls="--", lw=1)
        ax.set_xticks(range(len(PARALOGS)))
        ax.set_xticklabels([f"{a}/{b}" for a, b in PARALOGS], style="italic")
        ax.set_ylabel("expression ratio")
        ax.set_title(title, fontsize=10)

    fig.suptitle("Rod paralog ratios by clustering resolution\n"
                 "(one bar per rod cluster; the published result needs one cluster "
                 "above 1.0 for all three pairs and another below)", fontsize=11)
    fig.tight_layout()
    for ext in cfg.FIGURE_FORMATS:
        fig.savefig(RES_FIGURES / f"{name}.{ext}", format=ext,
                    dpi=cfg.FIGURE_DPI, bbox_inches="tight")
    print(f"Saved {name} -> figures/resolution_check/")
    return fig


def summarise(n_clusters: int, rod_clusters: list[str], verdict: dict) -> pd.DataFrame:
    """One row stating whether clustering granularity explains the failure."""
    n_split = sum(1 for k, v in verdict.items()
                  if k != "_patterns" and "inversion present" in v)
    patterns = verdict.get("_patterns", {})
    reproduced = {"AAA", "BBB"}.issubset(set(patterns.values()))
    df = pd.DataFrame([{
        "resolution": 0.2,
        "n_clusters": n_clusters,
        "n_clusters_published": 17,
        "n_rod_clusters": len(rod_clusters),
        "paralog_pairs_splitting": n_split,
        "coordinated_AAA_BBB_pattern": reproduced,
        "granularity_explains_failure": reproduced,
    }])
    print("\n" + "=" * 64)
    print("RESOLUTION SENSITIVITY — SUMMARY")
    print("=" * 64)
    print(df.T.to_string(header=False))
    print()
    if reproduced:
        print("The published rod split IS recovered at resolution 0.2. Clustering "
              "granularity explains the failure reported in §3.10, and the report "
              "should be revised to say so.")
    else:
        print("The published rod split is NOT recovered at resolution 0.2. "
              "Clustering granularity is eliminated as an explanation, alongside "
              "batch structure (§4.7). Normalisation remains the leading candidate.")
    return df
