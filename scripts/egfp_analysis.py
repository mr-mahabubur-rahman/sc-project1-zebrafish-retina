"""careg:EGFP reporter dynamics, Muller glia sub-states, rod and cone analysis.

Positivity rule, stated once so it can be quoted in the Methods section: a cell is
EGFP-positive if its RAW COUNT for the transgene feature is at least
config.EGFP_POSITIVE_MIN_COUNTS (default 1). Rationale: 3' UMI data are sparse and
a single-copy transgene yields low counts, so any higher cut would be arbitrary
without evidence from the observed distribution. The continuous log-normalised
value is reported alongside the binary call everywhere, so no conclusion rests on
the threshold alone.

Interpretive guard: EGFP marks *reporter activity*, not cell identity. A cell is
Muller glia because it expresses MG identity genes, not because it is EGFP+.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from . import config as cfg
from .io_utils import check_genes_available, find_egfp_feature
from .preprocessing import get_expression


# --------------------------------------------------------------------------- #
# EGFP positivity
# --------------------------------------------------------------------------- #

def add_egfp_columns(
    adata: ad.AnnData, egfp_name: str | None = None, min_counts: int | None = None
) -> tuple[ad.AnnData, str | None]:
    """Add egfp_counts, egfp_lognorm and egfp_positive to .obs.

    Returns the object and the resolved feature name (None if the transgene is
    absent, in which case no columns are added and downstream EGFP steps must be
    reported as not performed).
    """
    egfp_name = egfp_name or find_egfp_feature(adata)
    if egfp_name is None:
        return adata, None
    min_counts = cfg.EGFP_POSITIVE_MIN_COUNTS if min_counts is None else min_counts

    # The gene may have been dropped from .X/.layers by HVG subsetting, so fall
    # back in order: preserved obs column -> counts layer -> .raw.
    if "egfp_counts" in adata.obs:
        counts = adata.obs["egfp_counts"].to_numpy()
    elif egfp_name in adata.var_names and "counts" in adata.layers:
        counts = get_expression(adata, egfp_name, layer="counts")
    elif adata.raw is not None and egfp_name in adata.raw.var_names:
        raw_col = adata.raw[:, egfp_name].X
        counts = (np.asarray(raw_col.todense()).ravel()
                  if hasattr(raw_col, "todense") else np.asarray(raw_col).ravel())
        print("NOTE: EGFP counts recovered from .raw (log-normalised). The positivity "
              "call is still >0 vs 0, which is unaffected by the transform, but the "
              "count distribution plot will show normalised values.")
    else:
        raise KeyError(
            f"EGFP feature '{egfp_name}' is not reachable in .obs, .layers['counts'] "
            "or .raw. Run preprocessing.preserve_gene_counts_in_obs() before "
            "subsetting to HVGs."
        )
    adata.obs["egfp_counts"] = counts
    adata.obs["egfp_positive"] = counts >= min_counts
    try:
        adata.obs["egfp_lognorm"] = get_expression(adata, egfp_name, layer="lognorm")
    except KeyError:
        adata.obs["egfp_lognorm"] = np.log1p(counts)
        print("NOTE: lognorm layer unavailable for EGFP; used log1p(counts) for display only.")

    n_pos = int(adata.obs["egfp_positive"].sum())
    print(f"EGFP feature '{egfp_name}': {n_pos}/{adata.n_obs} cells positive "
          f"({100 * n_pos / adata.n_obs:.2f}%) at >= {min_counts} raw count(s). "
          f"Max counts in one cell: {counts.max():.0f}.")
    return adata, egfp_name


def egfp_count_distribution(adata: ad.AnnData) -> pd.DataFrame:
    """Histogram of raw EGFP counts, the evidence behind the threshold choice.

    If most positive cells carry a single count, a >0 threshold is the only
    defensible one; if there is a clear bimodal gap, a higher cut can be justified
    and must then be recorded as a deviation.
    """
    if "egfp_counts" not in adata.obs:
        raise KeyError("Run add_egfp_columns() first.")
    counts = adata.obs["egfp_counts"].astype(int)
    dist = counts.value_counts().sort_index().reset_index()
    dist.columns = ["egfp_counts", "n_cells"]
    dist["pct_of_all_cells"] = (100 * dist["n_cells"] / adata.n_obs).round(3)
    return dist


def egfp_summary(
    adata: ad.AnnData, group_keys: list[str] | None = None
) -> pd.DataFrame:
    """Per sample x condition x cell type EGFP summary for tables/egfp_summary.csv."""
    if "egfp_positive" not in adata.obs:
        raise KeyError("Run add_egfp_columns() first.")
    group_keys = group_keys or ["sample", "condition", "cell_type"]
    available = [k for k in group_keys if k in adata.obs]
    if len(available) < len(group_keys):
        print(f"Grouping only by {available}; missing {set(group_keys) - set(available)}.")

    grouped = adata.obs.groupby(available, observed=True)
    df = pd.DataFrame({
        "n_cells": grouped.size(),
        "egfp_positive": grouped["egfp_positive"].sum(),
        "mean_egfp": grouped["egfp_lognorm"].mean(),
        "median_egfp": grouped["egfp_lognorm"].median(),
    }).reset_index()
    df["egfp_positive_fraction"] = (df["egfp_positive"] / df["n_cells"]).round(4)
    df["egfp_positive_percent"] = (100 * df["egfp_positive_fraction"]).round(2)
    df[["mean_egfp", "median_egfp"]] = df[["mean_egfp", "median_egfp"]].round(4)
    return df


def egfp_by_condition(adata: ad.AnnData) -> pd.DataFrame:
    """EGFP+ fraction per condition, with the replicate-level values kept visible."""
    per_sample = (
        adata.obs.groupby(["condition", "sample"], observed=True)["egfp_positive"]
        .agg(["sum", "size"]).reset_index()
    )
    per_sample["percent_positive"] = (100 * per_sample["sum"] / per_sample["size"]).round(2)
    summary = (
        per_sample.groupby("condition", observed=True)["percent_positive"]
        .agg(["mean", "min", "max"]).round(2).reset_index()
        .rename(columns={"mean": "mean_percent_positive",
                         "min": "min_percent_positive", "max": "max_percent_positive"})
    )
    return per_sample.merge(summary, on="condition")


def egfp_enrichment_by_cluster(
    adata: ad.AnnData, cluster_key: str | None = None
) -> pd.DataFrame:
    """Which clusters concentrate EGFP+ cells, as an enrichment over the dataset rate.

    Enrichment > 1 means the cluster holds more EGFP+ cells than its size would give
    by chance. Reported alongside raw counts, because a high ratio on ten cells is
    not evidence.
    """
    cluster_key = cluster_key or ("cell_type" if "cell_type" in adata.obs else cfg.LEIDEN_KEY)
    overall = adata.obs["egfp_positive"].mean()
    grouped = adata.obs.groupby(cluster_key, observed=True)["egfp_positive"]
    df = pd.DataFrame({"n_cells": grouped.size(), "n_egfp_positive": grouped.sum()}).reset_index()
    df["percent_positive"] = (100 * df["n_egfp_positive"] / df["n_cells"]).round(2)
    df["enrichment_vs_dataset"] = (df["n_egfp_positive"] / df["n_cells"] / overall).round(2)
    df["share_of_all_egfp_cells"] = (
        100 * df["n_egfp_positive"] / df["n_egfp_positive"].sum()
    ).round(2)
    return df.sort_values("percent_positive", ascending=False)


def egfp_de_within_group(
    adata: ad.AnnData, group_value: str, group_key: str = "cell_type",
    min_cells: int = 30,
) -> pd.DataFrame:
    """EGFP+ vs EGFP- differential expression inside one cell type.

    This is the analogue of the paper's comparison of careg:EGFP-positive against
    negative Muller glia. Returns an empty frame with an explanation if either side
    has too few cells for the test to mean anything.
    """
    subset = adata[adata.obs[group_key] == group_value].copy()
    n_pos = int(subset.obs["egfp_positive"].sum())
    n_neg = int(subset.n_obs - n_pos)
    if min(n_pos, n_neg) < min_cells:
        print(
            f"Skipping DE for '{group_value}': {n_pos} EGFP+ and {n_neg} EGFP- cells; "
            f"fewer than {min_cells} on one side. Report this as underpowered rather "
            "than reporting the gene list."
        )
        return pd.DataFrame()

    subset.obs["egfp_group"] = np.where(subset.obs["egfp_positive"], "EGFP_pos", "EGFP_neg")
    use_raw = subset.raw is not None and "lognorm" not in subset.layers
    sc.tl.rank_genes_groups(
        subset, groupby="egfp_group", groups=["EGFP_pos"], reference="EGFP_neg",
        method=cfg.DE_METHOD, layer=None if use_raw else "lognorm", use_raw=use_raw,
    )
    df = sc.get.rank_genes_groups_df(subset, group="EGFP_pos")
    df = df.rename(columns={"names": "gene", "scores": "score",
                            "logfoldchanges": "logfoldchange",
                            "pvals": "pval", "pvals_adj": "pval_adj"})
    df.insert(0, "group", group_value)
    df.insert(1, "n_egfp_pos", n_pos)
    df.insert(2, "n_egfp_neg", n_neg)
    return df


# --------------------------------------------------------------------------- #
# Muller glia deep dive
# --------------------------------------------------------------------------- #

def subset_muller_glia(
    adata: ad.AnnData, cell_types: list[str] | None = None, group_key: str = "cell_type"
) -> ad.AnnData:
    """Subset to annotated Muller glia. Identity comes from the annotation, not EGFP."""
    cell_types = cell_types or [
        ct for ct in adata.obs[group_key].cat.categories if "muller" in str(ct).lower()
    ]
    if not cell_types:
        raise ValueError(
            "No annotated Muller glia cluster to subset. If MG were not identified, "
            "say so in the report; do not substitute EGFP+ cells for them."
        )
    mg = adata[adata.obs[group_key].isin(cell_types)].copy()
    print(f"Muller glia subset: {mg.n_obs} cells from {cell_types}")
    if mg.n_obs < 200:
        print("NOTE: fewer than 200 MG cells. Sub-state analysis on this many cells is "
              "exploratory; report it as such.")
    return mg


def recluster_subset(
    adata: ad.AnnData, resolution: float = 0.3, n_pcs: int = 15,
    n_neighbors: int = 15, n_top_genes: int = 2000,
) -> ad.AnnData:
    """Re-run HVG/PCA/UMAP/Leiden within a subset, starting from stored counts.

    Feature selection is redone inside the subset: HVGs chosen across the whole
    retina are dominated by between-cell-type variation and are the wrong features
    for resolving states within one cell type.
    """
    # Rebuild from .raw, which holds log-normalised values for EVERY gene. Working
    # from the HVG-subset object instead would restrict the sub-analysis to the
    # 2,000 genes chosen for global structure -- which is exactly the wrong feature
    # set for resolving states inside one cell type, and would silently drop
    # markers such as meig1, ppdpfa and pcna.
    if adata.raw is not None:
        source = adata.raw.to_adata()
        source.obs = adata.obs.copy()
        print(f"Sub-analysis rebuilt from .raw: {source.n_vars} genes available "
              f"(vs {adata.n_vars} in the HVG-subset object).")
    elif "lognorm" in adata.layers:
        source = ad.AnnData(X=adata.layers["lognorm"].copy(),
                            obs=adata.obs.copy(), var=adata.var.copy())
    else:
        raise KeyError(
            "Subset has neither .raw nor a 'lognorm' layer, so the full gene set is "
            "unreachable. Re-subset from the annotated checkpoint."
        )
    sub = source
    sub.layers["lognorm"] = sub.X.copy()
    if "counts" in adata.layers and adata.n_vars == sub.n_vars:
        sub.layers["counts"] = adata.layers["counts"].copy()
    sub.raw = sub
    sc.pp.highly_variable_genes(
        sub, n_top_genes=min(n_top_genes, sub.n_vars - 1),
        batch_key="sample" if sub.obs["sample"].nunique() > 1 else None,
    )
    sub = sub[:, sub.var["highly_variable"]].copy()
    sc.pp.scale(sub, max_value=cfg.SCALE_MAX_VALUE)
    sc.tl.pca(sub, n_comps=min(50, sub.n_vars - 1, sub.n_obs - 1),
              svd_solver="arpack", random_state=cfg.RANDOM_SEED)
    sc.pp.neighbors(sub, n_neighbors=min(n_neighbors, sub.n_obs - 1),
                    n_pcs=min(n_pcs, sub.obsm["X_pca"].shape[1]),
                    random_state=cfg.RANDOM_SEED)
    sc.tl.umap(sub, random_state=cfg.RANDOM_SEED)
    try:
        sc.tl.leiden(sub, resolution=resolution, key_added="subcluster",
                     random_state=cfg.RANDOM_SEED, flavor="igraph",
                     n_iterations=2, directed=False)
    except TypeError:
        sc.tl.leiden(sub, resolution=resolution, key_added="subcluster",
                     random_state=cfg.RANDOM_SEED)
    print(f"Subset re-clustered at resolution {resolution}: "
          f"{sub.obs['subcluster'].nunique()} subclusters, {sub.n_obs} cells.")
    return sub


def mg_state_evidence(
    mg: ad.AnnData, subcluster_key: str = "subcluster"
) -> pd.DataFrame:
    """Evidence table for calling resting / activated / proliferative MG sub-states.

    Reports mean expression of identity, activation and proliferation markers plus
    EGFP and condition composition per subcluster. It does NOT assign the three
    states: if the markers do not separate, the honest answer is that the data do
    not resolve three populations.
    """
    genes = cfg.MG_IDENTITY_MARKERS + cfg.ACTIVATION_MARKERS
    source = mg.raw.to_adata() if mg.raw is not None else mg
    present, _ = check_genes_available(source, genes, label="MG state panel")

    X = source[:, present].X
    dense = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
    df = pd.DataFrame(dense, columns=present, index=mg.obs_names)
    df[subcluster_key] = mg.obs[subcluster_key].values
    expr = df.groupby(subcluster_key, observed=True).mean().round(3)

    extra = pd.DataFrame({
        "n_cells": mg.obs.groupby(subcluster_key, observed=True).size(),
    })
    if "egfp_positive" in mg.obs:
        extra["pct_egfp_positive"] = (
            100 * mg.obs.groupby(subcluster_key, observed=True)["egfp_positive"].mean()
        ).round(2)
    comp = (
        mg.obs.groupby([subcluster_key, "condition"], observed=True).size()
        .unstack(fill_value=0)
    )
    comp = (100 * comp.div(comp.sum(axis=1), axis=0)).round(1)
    comp.columns = [f"pct_{c}" for c in comp.columns]
    return extra.join(comp).join(expr)


MG_STATE_GUIDANCE = (
    "Reading the MG evidence table: a resting state should be high on identity "
    "markers (rlbp1a, glula/glulb, slc1a3b) and low on pcna/mki67/ascl1a; an "
    "activated state keeps identity markers while gaining ascl1a/her4.1/EGFP; a "
    "progenitor-like state loses identity markers and gains pcna/mki67. If no "
    "subcluster shows the third pattern, report two states, not three. Note also "
    "that the paper found EGFP mostly in activated MG and largely ABSENT from the "
    "proliferating progenitor clusters, so EGFP+ and pcna+ need not coincide."
)


# --------------------------------------------------------------------------- #
# Photoreceptors
# --------------------------------------------------------------------------- #

def rod_heterogeneity_evidence(
    adata: ad.AnnData, rod_types: list[str] | None = None, group_key: str = "cell_type",
    resolution: float = 0.3,
) -> tuple[ad.AnnData, pd.DataFrame]:
    """Re-cluster rods and score the paper's immature/mature axis without assuming it.

    Returns the rod subset and a table of the paralog pairs that distinguished the
    two rod clusters in Bise et al. (rho/rhol, pde6ga/pde6gb, guca1a/guca1b) plus
    the immature candidates meig1, ppdpfa, rom1b. Whether the split exists is a
    result to report either way.
    """
    rod_types = rod_types or [
        ct for ct in adata.obs[group_key].cat.categories if "rod" in str(ct).lower()
    ]
    if not rod_types:
        raise ValueError("No annotated rod cluster to analyse.")
    rods = adata[adata.obs[group_key].isin(rod_types)].copy()
    rods = recluster_subset(rods, resolution=resolution)

    panel = ["rho", "rhol", "pde6ga", "pde6gb", "guca1a", "guca1b",
             "meig1", "ppdpfa", "rom1b", "gnat1", "gngt1", "nr2e3",
             "sagb", "rcvrna", "gnb1a", "crx", "prom1"]
    # Check availability against the FULL gene set (.raw), not the HVG subset:
    # meig1, rom1b and guca1b are rarely selected as highly variable, and checking
    # the subset first silently dropped exactly the genes the paper's rod
    # comparison depends on.
    source = rods.raw.to_adata() if rods.raw is not None else rods
    present, _ = check_genes_available(source, panel, label="rod maturity panel")
    X = source[:, present].X
    dense = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
    df = pd.DataFrame(dense, columns=present, index=rods.obs_names)
    df["subcluster"] = rods.obs["subcluster"].values
    expr = df.groupby("subcluster", observed=True).mean().round(3)

    comp = rods.obs.groupby(["subcluster", "condition"], observed=True).size().unstack(fill_value=0)
    comp_pct = (100 * comp.div(comp.sum(axis=1), axis=0)).round(1)
    comp_pct.columns = [f"pct_{c}" for c in comp_pct.columns]
    n_cells = pd.DataFrame({"n_cells": rods.obs.groupby("subcluster", observed=True).size()})
    return rods, n_cells.join(comp_pct).join(expr)


def cone_subtype_evidence(
    adata: ad.AnnData, cone_types: list[str] | None = None, group_key: str = "cell_type",
    resolution: float = 0.3,
) -> tuple[ad.AnnData, pd.DataFrame]:
    """Re-cluster cones and report opsin and arrestin usage per subcluster."""
    cone_types = cone_types or [
        ct for ct in adata.obs[group_key].cat.categories if "cone" in str(ct).lower()
    ]
    if not cone_types:
        raise ValueError("No annotated cone cluster to analyse.")
    cones = adata[adata.obs[group_key].isin(cone_types)].copy()
    cones = recluster_subset(cones, resolution=resolution)

    panel = ["opn1sw1", "opn1sw2", "opn1mw1", "opn1mw2", "opn1lw1", "opn1lw2",
             "arr3a", "arr3b", "tbx2a", "cngb3.2", "guca1e", "gnat2", "pde6c"]
    source = cones.raw.to_adata() if cones.raw is not None else cones
    present, _ = check_genes_available(source, panel, label="cone subtype panel")
    X = source[:, present].X
    dense = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
    df = pd.DataFrame(dense, columns=present, index=cones.obs_names)
    df["subcluster"] = cones.obs["subcluster"].values
    expr = df.groupby("subcluster", observed=True).mean().round(3)

    comp = cones.obs.groupby(["subcluster", "condition"], observed=True).size().unstack(fill_value=0)
    comp_pct = (100 * comp.div(comp.sum(axis=1), axis=0)).round(1)
    comp_pct.columns = [f"pct_{c}" for c in comp_pct.columns]
    n_cells = pd.DataFrame({"n_cells": cones.obs.groupby("subcluster", observed=True).size()})
    return cones, n_cells.join(comp_pct).join(expr)


def injury_response_de(
    adata: ad.AnnData, cell_type: str, condition: str,
    reference: str = "ctrl", group_key: str = "cell_type", min_cells: int = 30,
) -> pd.DataFrame:
    """DE within one cell type, one post-injury timepoint vs control.

    Mirrors the paper's per-cell-type comparisons. Note the limitation to state in
    the report: with two replicates per condition, a cell-level Wilcoxon test
    treats cells as independent samples and so overstates significance. Treat the
    ranking as descriptive, and say so.
    """
    subset = adata[
        (adata.obs[group_key] == cell_type)
        & (adata.obs["condition"].isin([condition, reference]))
    ].copy()
    counts = subset.obs["condition"].value_counts()
    if counts.get(condition, 0) < min_cells or counts.get(reference, 0) < min_cells:
        print(f"Skipping {cell_type} {condition} vs {reference}: "
              f"{counts.get(condition, 0)} vs {counts.get(reference, 0)} cells.")
        return pd.DataFrame()

    subset.obs["_grp"] = subset.obs["condition"].astype(str)
    use_raw = subset.raw is not None and "lognorm" not in subset.layers
    sc.tl.rank_genes_groups(
        subset, groupby="_grp", groups=[condition], reference=reference,
        method=cfg.DE_METHOD, layer=None if use_raw else "lognorm", use_raw=use_raw,
    )
    df = sc.get.rank_genes_groups_df(subset, group=condition)
    df = df.rename(columns={"names": "gene", "scores": "score",
                            "logfoldchanges": "logfoldchange",
                            "pvals": "pval", "pvals_adj": "pval_adj"})
    df.insert(0, "cell_type", cell_type)
    df.insert(1, "comparison", f"{condition}_vs_{reference}")
    df.insert(2, "n_cells_test", int(counts.get(condition, 0)))
    df.insert(3, "n_cells_reference", int(counts.get(reference, 0)))
    return df
