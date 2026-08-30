"""Publication-quality plotting with one shared style and deterministic filenames.

Every figure function returns a matplotlib Figure and writes both a PNG (300 dpi)
and a vector PDF into the correct figures/ subdirectory. Colours come from
config.CELL_TYPE_COLORS and config.CONDITION_COLORS so a cell type keeps the same
colour in Figures 2, 3, 4 and 5.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from . import config as cfg


# --------------------------------------------------------------------------- #
# Style and saving
# --------------------------------------------------------------------------- #

def set_style() -> None:
    """Apply the project-wide matplotlib and scanpy defaults."""
    mpl.rcParams.update({
        "figure.dpi": 110,
        "savefig.dpi": cfg.FIGURE_DPI,
        "savefig.bbox": "tight",
        "font.size": cfg.BASE_FONT_SIZE,
        "axes.titlesize": cfg.BASE_FONT_SIZE + 1,
        "axes.labelsize": cfg.BASE_FONT_SIZE,
        "xtick.labelsize": cfg.BASE_FONT_SIZE - 1,
        "ytick.labelsize": cfg.BASE_FONT_SIZE - 1,
        "legend.fontsize": cfg.BASE_FONT_SIZE - 1,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,   # editable text in the vector output
        "ps.fonttype": 42,
        "figure.facecolor": "white",
    })
    sc.settings.verbosity = 1
    try:
        sc.set_figure_params(dpi=110, facecolor="white", frameon=False)
    except AttributeError:  # scanpy < 1.11
        sc.settings.set_figure_params(dpi=110, facecolor="white", frameon=False)
    sc.settings.autoshow = False



def title_clear_of_content(fig: plt.Figure, text: str, position: str = "top",
                           pad: float = 0.03, fontsize: int = 12) -> None:
    """Place a title outside everything already drawn, measured rather than guessed.

    Scanpy's dotplot and matrixplot carry rotated group-bracket labels whose height
    depends on the label text -- "Retinal ganglion cells" reaches much higher than
    "Bipolar cells". A fixed y offset therefore either overlaps the labels or
    floats far above them. This renders once, measures the tight bounding box of
    all artists, and places the title just beyond it.

    Pass position="bottom" to put the title under the x tick labels instead.
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    bbox = fig.get_tightbbox(renderer)          # inches, from the figure origin
    height = fig.get_figheight()
    if position == "top":
        fig.text(0.5, bbox.y1 / height + pad, text,
                 ha="center", va="bottom", fontsize=fontsize)
    else:
        fig.text(0.5, bbox.y0 / height - pad, text,
                 ha="center", va="top", fontsize=fontsize)


def save_figure(fig: plt.Figure, name: str, subdir: str) -> list[Path]:
    """Write a figure as PNG and PDF into figures/<subdir>/ and return the paths."""
    if subdir not in cfg.FIGURE_DIRS:
        raise KeyError(f"Unknown figure subdir '{subdir}'. Known: {sorted(cfg.FIGURE_DIRS)}")
    out_dir = cfg.FIGURE_DIRS[subdir]
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for ext in cfg.FIGURE_FORMATS:
        path = out_dir / f"{name}.{ext}"
        fig.savefig(path, format=ext, dpi=cfg.FIGURE_DPI, bbox_inches="tight")
        paths.append(path)
    print(f"Saved {name} -> {out_dir.relative_to(cfg.REPO_ROOT)}/ ({', '.join(cfg.FIGURE_FORMATS)})")
    return paths


def cell_type_palette(categories) -> list[str]:
    """Colour list aligned to the given cell-type categories, grey for unknowns."""
    return [cfg.CELL_TYPE_COLORS.get(str(c), "#BBBBBB") for c in categories]


def condition_palette(categories) -> list[str]:
    return [cfg.CONDITION_COLORS.get(str(c), "#BBBBBB") for c in categories]


def apply_palettes(adata: ad.AnnData) -> ad.AnnData:
    """Store the project palettes in .uns so scanpy plots pick them up automatically."""
    if "cell_type" in adata.obs and hasattr(adata.obs["cell_type"], "cat"):
        adata.uns["cell_type_colors"] = cell_type_palette(adata.obs["cell_type"].cat.categories)
    if "condition" in adata.obs and hasattr(adata.obs["condition"], "cat"):
        adata.uns["condition_colors"] = condition_palette(adata.obs["condition"].cat.categories)
    return adata


# --------------------------------------------------------------------------- #
# Figure 1 -- QC
# --------------------------------------------------------------------------- #

def figure_01_qc(adata: ad.AnnData, filter_log: pd.DataFrame | None = None,
                 name: str = "figure_01_qc_overview") -> plt.Figure:
    """QC panel: four metric violins by sample, a counts-vs-genes scatter, filter summary."""
    metrics = ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo"]
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    axes = axes.ravel()

    samples = list(adata.obs["sample"].cat.categories)
    colors = [cfg.CONDITION_COLORS.get(
        str(adata.obs.loc[adata.obs["sample"] == s, "condition"].iloc[0]), "#999999")
        for s in samples]

    for ax, metric in zip(axes[:4], metrics):
        data = [adata.obs.loc[adata.obs["sample"] == s, metric].values for s in samples]
        parts = ax.violinplot(data, showextrema=False, widths=0.85)
        for body, color in zip(parts["bodies"], colors):
            body.set_facecolor(color)
            body.set_alpha(0.75)
        ax.boxplot(data, widths=0.12, showfliers=False,
                   medianprops={"color": "black", "linewidth": 1})
        ax.set_xticks(range(1, len(samples) + 1))
        ax.set_xticklabels(samples, rotation=45, ha="right")
        ax.set_title(metric.replace("_", " "))
        ax.set_ylabel(metric)
        if metric == "total_counts":
            ax.set_yscale("log")
    axes[2].axhline(cfg.MAX_PCT_MT, color="red", ls="--", lw=1,
                    label=f"threshold {cfg.MAX_PCT_MT}%")
    axes[2].legend(frameon=False)
    axes[0].axhline(cfg.MIN_GENES_PER_CELL, color="red", ls="--", lw=1,
                    label=f"threshold {cfg.MIN_GENES_PER_CELL}")
    axes[0].legend(frameon=False)

    ax = axes[4]
    sctr = ax.scatter(adata.obs["total_counts"], adata.obs["n_genes_by_counts"],
                      c=adata.obs["pct_counts_mt"], s=1.5, alpha=0.35,
                      cmap="viridis", rasterized=True)
    ax.set_xscale("log")
    ax.set_xlabel("total counts (log scale)")
    ax.set_ylabel("genes detected")
    ax.set_title("Counts vs genes, coloured by % mitochondrial")
    fig.colorbar(sctr, ax=ax, label="% mito")

    ax = axes[5]
    if filter_log is not None and len(filter_log):
        ax.barh(filter_log["step"], filter_log["cells"], color="#4D4D4D")
        ax.set_xlabel("cells remaining")
        ax.set_title("Cells retained at each filtering step")
        ax.invert_yaxis()
    else:
        ax.axis("off")

    fig.suptitle("Figure 1 | Quality control across the eight samples", fontsize=12)
    fig.tight_layout()
    save_figure(fig, name, "qc")
    return fig


# --------------------------------------------------------------------------- #
# Figure 2 -- global atlas
# --------------------------------------------------------------------------- #

def figure_02_umap(adata: ad.AnnData, cluster_key: str | None = None,
                   name: str = "figure_02_umap_atlas") -> plt.Figure:
    """UMAP coloured by cluster (A), condition (B) and sample (C)."""
    cluster_key = cluster_key or cfg.LEIDEN_KEY
    apply_palettes(adata)
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
    sc.pl.umap(adata, color=cluster_key, ax=axes[0], show=False,
               legend_loc="on data", legend_fontsize=7, frameon=False,
               title="A  Leiden clusters")
    sc.pl.umap(adata, color="condition", ax=axes[1], show=False, frameon=False,
               title="B  Condition", palette=condition_palette(adata.obs["condition"].cat.categories))
    sc.pl.umap(adata, color="sample", ax=axes[2], show=False, frameon=False,
               title="C  Sample (replicate identity retained)")
    fig.suptitle("Figure 2 | Global cell atlas", fontsize=12)
    fig.tight_layout()
    save_figure(fig, name, "umap")
    return fig


def figure_02_split_by_condition(adata: ad.AnnData, cluster_key: str | None = None,
                                 name: str = "figure_02_umap_split_by_condition") -> plt.Figure:
    """One UMAP panel per condition -- the clearest view of batch vs biology."""
    cluster_key = cluster_key or cfg.LEIDEN_KEY
    conditions = list(adata.obs["condition"].cat.categories)
    fig, axes = plt.subplots(1, len(conditions), figsize=(4.5 * len(conditions), 4.5))
    axes = np.atleast_1d(axes)
    coords = adata.obsm["X_umap"]
    for ax, cond in zip(axes, conditions):
        mask = (adata.obs["condition"] == cond).values
        ax.scatter(coords[:, 0], coords[:, 1], s=1, c="#E8E8E8", rasterized=True)
        ax.scatter(coords[mask, 0], coords[mask, 1], s=1.5,
                   c=cfg.CONDITION_COLORS.get(str(cond), "#B2182B"), rasterized=True)
        ax.set_title(f"{cfg.CONDITION_LABELS.get(str(cond), cond)}  (n={int(mask.sum())})")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
    fig.suptitle("Figure 2D | Cells by timepoint on the shared embedding", fontsize=12)
    fig.tight_layout()
    save_figure(fig, name, "umap")
    return fig


# --------------------------------------------------------------------------- #
# Figure 3 -- annotation
# --------------------------------------------------------------------------- #

def figure_03_markers(adata: ad.AnnData, marker_panel: dict[str, list[str]],
                      cluster_key: str | None = None,
                      name: str = "figure_03_marker_dotplot") -> plt.Figure:
    """Dotplot of the canonical panel across clusters."""
    cluster_key = cluster_key or cfg.LEIDEN_KEY
    usable = {ct: [g for g in genes if g in adata.var_names]
              for ct, genes in marker_panel.items()}
    usable = {ct: genes for ct, genes in usable.items() if genes}
    dp = sc.pl.dotplot(adata, var_names=usable, groupby=cluster_key,
                       standard_scale="var", show=False, return_fig=True)
    fig = dp.get_axes()["mainplot_ax"].get_figure()
    title_clear_of_content(fig, "Figure 3A | Canonical marker expression")
    save_figure(fig, name, "markers")
    return fig


def figure_03_heatmap(adata: ad.AnnData, marker_panel: dict[str, list[str]],
                      cluster_key: str | None = None,
                      name: str = "figure_03_marker_heatmap") -> plt.Figure:
    """Heatmap of the same panel, mean-scaled per gene."""
    cluster_key = cluster_key or cfg.LEIDEN_KEY
    usable = {ct: [g for g in genes if g in adata.var_names]
              for ct, genes in marker_panel.items()}
    usable = {ct: genes for ct, genes in usable.items() if genes}
    hm = sc.pl.matrixplot(adata, var_names=usable, groupby=cluster_key,
                          standard_scale="var", cmap="viridis",
                          show=False, return_fig=True)
    fig = hm.get_axes()["mainplot_ax"].get_figure()
    title_clear_of_content(fig, "Figure 3B | Mean scaled marker expression")
    save_figure(fig, name, "markers")
    return fig


def figure_03_annotated_umap(adata: ad.AnnData,
                             name: str = "figure_03_annotated_umap") -> plt.Figure:
    """UMAP coloured by the manually assigned cell type."""
    apply_palettes(adata)
    fig, ax = plt.subplots(figsize=(8, 6.5))
    sc.pl.umap(adata, color="cell_type", ax=ax, show=False, frameon=False,
               palette=cell_type_palette(adata.obs["cell_type"].cat.categories),
               title="Figure 3C | Annotated cell types")
    fig.tight_layout()
    save_figure(fig, name, "markers")
    return fig


# --------------------------------------------------------------------------- #
# Figure 4 -- EGFP and Muller glia
# --------------------------------------------------------------------------- #

def figure_04_egfp(adata: ad.AnnData, egfp_by_cond: pd.DataFrame,
                   name: str = "figure_04_egfp_dynamics") -> plt.Figure:
    """Panels A-D: EGFP on UMAP, distribution by timepoint, positive fraction, by cell type."""
    apply_palettes(adata)
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.25)

    ax_a = fig.add_subplot(gs[0, 0])
    sc.pl.umap(adata, color="egfp_lognorm", ax=ax_a, show=False, frameon=False,
               cmap="Reds", title="A  EGFP expression (log-normalised)")

    ax_b = fig.add_subplot(gs[0, 1])
    conditions = list(adata.obs["condition"].cat.categories)
    data = [adata.obs.loc[adata.obs["condition"] == c, "egfp_lognorm"].values
            for c in conditions]
    parts = ax_b.violinplot(data, showextrema=False, widths=0.85)
    for body, cond in zip(parts["bodies"], conditions):
        body.set_facecolor(cfg.CONDITION_COLORS.get(str(cond), "#999999"))
        body.set_alpha(0.8)
    ax_b.set_xticks(range(1, len(conditions) + 1))
    ax_b.set_xticklabels([cfg.CONDITION_LABELS.get(str(c), c) for c in conditions])
    ax_b.set_ylabel("EGFP (log-normalised)")
    ax_b.set_title("B  EGFP expression across the injury time course")

    ax_c = fig.add_subplot(gs[1, 0])
    per_sample = egfp_by_cond.drop_duplicates(subset=["sample"])
    means = (egfp_by_cond.drop_duplicates(subset=["condition"])
             .set_index("condition")["mean_percent_positive"])
    x = np.arange(len(conditions))
    ax_c.bar(x, [means.get(c, np.nan) for c in conditions],
             color=[cfg.CONDITION_COLORS.get(str(c), "#999999") for c in conditions],
             alpha=0.85, label="condition mean")
    for i, cond in enumerate(conditions):
        pts = per_sample.loc[per_sample["condition"] == cond, "percent_positive"]
        ax_c.scatter(np.full(len(pts), i), pts, color="black", zorder=3, s=25,
                     label="replicate" if i == 0 else None)
    ax_c.set_xticks(x)
    ax_c.set_xticklabels([cfg.CONDITION_LABELS.get(str(c), c) for c in conditions])
    ax_c.set_ylabel("% EGFP-positive cells")
    ax_c.set_title("C  EGFP-positive fraction (replicates shown individually)")
    ax_c.legend(frameon=False)

    ax_d = fig.add_subplot(gs[1, 1])
    group_key = "cell_type" if "cell_type" in adata.obs else cfg.LEIDEN_KEY
    by_type = (adata.obs.groupby(group_key, observed=True)["egfp_positive"]
               .mean().mul(100).sort_values())
    ax_d.barh(by_type.index.astype(str), by_type.values,
              color=[cfg.CELL_TYPE_COLORS.get(str(c), "#999999") for c in by_type.index])
    ax_d.set_xlabel("% EGFP-positive cells")
    ax_d.set_title("D  EGFP-positive fraction by cell type")

    fig.suptitle("Figure 4 | careg:EGFP reporter dynamics", fontsize=12)
    save_figure(fig, name, "egfp_mg")
    return fig


def figure_04_activation_markers(adata: ad.AnnData, genes: list[str] | None = None,
                                 name: str = "figure_04_egfp_vs_proliferation") -> plt.Figure:
    """Panel E: EGFP alongside proliferation and activation markers, by cell type."""
    genes = genes or (["egfp_lognorm"] + cfg.ACTIVATION_MARKERS)
    var_genes = [g for g in genes if g in adata.var_names]
    obs_genes = [g for g in genes if g in adata.obs.columns]
    plot_vars = obs_genes + var_genes
    if not plot_vars:
        raise ValueError("None of the requested activation markers are available.")
    group_key = "cell_type" if "cell_type" in adata.obs else cfg.LEIDEN_KEY
    dp = sc.pl.dotplot(adata, var_names=plot_vars, groupby=group_key,
                       standard_scale="var", show=False, return_fig=True)
    fig = dp.get_axes()["mainplot_ax"].get_figure()
    title_clear_of_content(fig, "Figure 4E | EGFP and activation/proliferation markers")
    save_figure(fig, name, "egfp_mg")
    return fig


def figure_04_mg_umap(mg: ad.AnnData, name: str = "figure_04_mg_substates") -> plt.Figure:
    """Muller glia sub-analysis: subclusters, condition, EGFP, proliferation."""
    panels = ["subcluster", "condition"]
    if "egfp_lognorm" in mg.obs:
        panels.append("egfp_lognorm")
    for gene in ("pcna", "mki67", "ascl1a"):
        if gene in mg.var_names or (mg.raw is not None and gene in mg.raw.var_names):
            panels.append(gene)
    n = len(panels)
    ncols = 3
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.5 * ncols, 4.8 * nrows))
    axes = np.atleast_1d(axes).ravel()
    for ax, panel in zip(axes, panels):
        sc.pl.umap(mg, color=panel, ax=ax, show=False, frameon=False, title=panel)
    for ax in axes[n:]:
        ax.axis("off")
    fig.suptitle("Figure 4F | Muller glia sub-states", fontsize=12)
    fig.tight_layout()
    save_figure(fig, name, "egfp_mg")
    return fig


# --------------------------------------------------------------------------- #
# Figure 5 -- composition
# --------------------------------------------------------------------------- #

def figure_05_proportions(proportions: pd.DataFrame,
                          by_condition: pd.DataFrame,
                          name: str = "figure_05_cell_type_proportions") -> plt.Figure:
    """Stacked composition per sample (A) and per-cell-type trajectories with replicates (B)."""
    fig, axes = plt.subplots(1, 2, figsize=(17, 6))

    wide = proportions.pivot_table(index="sample", columns="cell_type",
                                   values="percent", observed=True, fill_value=0)
    sample_order = [s for s in cfg.EXPECTED_SAMPLES if s in wide.index]
    wide = wide.loc[sample_order] if sample_order else wide
    bottom = np.zeros(len(wide))
    for cell_type in wide.columns:
        axes[0].bar(wide.index.astype(str), wide[cell_type], bottom=bottom,
                    color=cfg.CELL_TYPE_COLORS.get(str(cell_type), "#BBBBBB"),
                    label=str(cell_type), width=0.8)
        bottom += wide[cell_type].values
    axes[0].set_ylabel("% of cells in sample")
    axes[0].set_title("A  Composition per sample (replicates side by side)")
    axes[0].tick_params(axis="x", rotation=45)
    axes[0].legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False, fontsize=7)

    ax = axes[1]
    x_pos = {c: i for i, c in enumerate(cfg.CONDITION_ORDER)}
    for cell_type, grp in by_condition.groupby("cell_type", observed=True):
        grp = grp.set_index("condition").reindex(cfg.CONDITION_ORDER).reset_index()
        color = cfg.CELL_TYPE_COLORS.get(str(cell_type), "#BBBBBB")
        ax.plot([x_pos[c] for c in grp["condition"]], grp["mean_percent"],
                marker="o", color=color, label=str(cell_type), lw=1.8)
        ax.fill_between([x_pos[c] for c in grp["condition"]],
                        grp["min_percent"], grp["max_percent"], color=color, alpha=0.18)
    reps = proportions.copy()
    ax.scatter([x_pos[str(c)] for c in reps["condition"]], reps["percent"],
               s=10, color="black", alpha=0.5, zorder=3, label="replicate values")
    ax.set_xticks(range(len(cfg.CONDITION_ORDER)))
    ax.set_xticklabels([cfg.CONDITION_LABELS[c] for c in cfg.CONDITION_ORDER])
    ax.set_ylabel("% of cells")
    ax.set_title("B  Cell-type proportion across the time course\n"
                 "(line = condition mean, band = replicate range)")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False, fontsize=7)

    fig.suptitle("Figure 5 | Cell-type dynamics -- capture-biased, see caption", fontsize=12)
    fig.tight_layout()
    save_figure(fig, name, "proportions")
    return fig


def figure_hvg(adata: ad.AnnData, name: str = "figure_s1_hvg") -> plt.Figure:
    """Supplementary HVG diagnostic."""
    sc.pl.highly_variable_genes(adata, show=False)
    fig = plt.gcf()
    fig.suptitle("Supplementary | Highly variable gene selection", fontsize=11)
    save_figure(fig, name, "qc")
    return fig


def figure_pca_variance(variance_table: pd.DataFrame, n_pcs_used: int,
                        name: str = "figure_s2_pca_variance") -> plt.Figure:
    """Supplementary elbow plot with the chosen n_pcs marked."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    axes[0].plot(variance_table["PC"], variance_table["variance_ratio"], marker="o", ms=3)
    axes[0].axvline(n_pcs_used, color="red", ls="--", lw=1, label=f"n_pcs = {n_pcs_used}")
    axes[0].set_yscale("log")
    axes[0].set_xlabel("principal component")
    axes[0].set_ylabel("variance ratio (log)")
    axes[0].set_title("Variance explained per PC")
    axes[0].legend(frameon=False)

    axes[1].plot(variance_table["PC"], 100 * variance_table["cumulative_variance_ratio"],
                 marker="o", ms=3, color="#4D4D4D")
    axes[1].axvline(n_pcs_used, color="red", ls="--", lw=1)
    axes[1].set_xlabel("principal component")
    axes[1].set_ylabel("cumulative variance explained (%)")
    axes[1].set_title("Cumulative variance")
    fig.tight_layout()
    save_figure(fig, name, "umap")
    return fig
