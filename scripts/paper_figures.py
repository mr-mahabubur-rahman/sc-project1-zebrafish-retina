"""Render the paper's scRNA-seq figure panels from this reanalysis.

Covers the four figures in Bise et al. 2023 that have a single-cell counterpart:

    Figure 5  C  UMAP atlas            D  canonical marker dot plot
              E  cluster percentage per timepoint
    Figure 6  B  rod paralog dot plot  C/D GO terms   E  gene heatmap
              F  volcano, injured vs control
    Figure 7  B  opsin heatmap         D  cone gene heatmap
              E  volcano               F  GO terms
    Figure 8  A  EGFP+ cells per condition on UMAP    B  EGFP+ counts
              C  EGFP+ percentage per cluster         E  EGFP+/- MG heatmap

Figures 1-4 and 9-11 of the paper are immunofluorescence, TOR signalling and
rapamycin experiments. No transcriptomic analysis can produce them, and nothing
here attempts to.

Honesty constraint carried through every function: a UMAP is not reproducible as
an image. Coordinates depend on normalisation, feature selection and random
initialisation, so panel C here will not superimpose on the paper's panel C even
though it shows the same populations. What is comparable is the *content* -- which
cell types exist, their relative abundance, which genes mark them, and the
percentages. Compare numbers, not pixels.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

from . import config as cfg
from .io_utils import check_genes_available
from .plotting import save_figure

# The paper's figures live in their own directory so they are never confused with
# the project's own Figures 1-5.
PAPER_FIG_DIR = cfg.FIGURES_DIR / "paper_figure_reproduction"


def _save(fig: plt.Figure, name: str) -> plt.Figure:
    """Write a panel as PNG + PDF into figures/paper_figure_reproduction/."""
    PAPER_FIG_DIR.mkdir(parents=True, exist_ok=True)
    for ext in cfg.FIGURE_FORMATS:
        fig.savefig(PAPER_FIG_DIR / f"{name}.{ext}", format=ext,
                    dpi=cfg.FIGURE_DPI, bbox_inches="tight")
    print(f"Saved {name} -> figures/paper_figure_reproduction/")
    return fig


def _dense(adata: AnnData, genes: list[str], use_raw: bool = True) -> pd.DataFrame:
    """Dense expression frame for the given genes, from .raw when available."""
    source = adata.raw.to_adata() if (use_raw and adata.raw is not None) else adata
    present = [g for g in genes if g in source.var_names]
    if not present:
        return pd.DataFrame(index=adata.obs_names)
    X = source[:, present].X
    dense = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
    return pd.DataFrame(dense, columns=present, index=adata.obs_names)


# --------------------------------------------------------------------------- #
# Figure 5 -- atlas
# --------------------------------------------------------------------------- #

def figure_5c_umap(adata: AnnData, group_key: str = "cell_type") -> plt.Figure:
    """Paper Fig 5C: UMAP of the integrated data, coloured by cell population."""
    fig, ax = plt.subplots(figsize=(7.5, 6))
    sc.pl.umap(adata, color=group_key, ax=ax, show=False, frameon=False,
               legend_fontsize=7, title="Fig 5C equivalent | cell populations")
    ax.text(0.01, -0.06,
            "Embedding coordinates are method-dependent and are not expected to "
            "match the published UMAP.",
            transform=ax.transAxes, fontsize=6.5, color="#555555")
    fig.tight_layout()
    return _save(fig, "fig5C_umap_atlas")


def figure_5d_marker_dotplot(adata: AnnData, group_key: str = "cell_type",
                             markers: dict | None = None) -> plt.Figure:
    """Paper Fig 5D: canonical marker dot plot, dot size = fraction expressing."""
    markers = markers or cfg.MARKERS_PAPER
    source = adata.raw.to_adata() if adata.raw is not None else adata
    source.obs = adata.obs.copy()
    usable = {}
    for cell_type, genes in markers.items():
        present = [g for g in genes if g in source.var_names]
        if present:
            usable[cell_type] = present
        else:
            print(f"[Fig 5D] no marker available for '{cell_type}'; omitted.")
    dp = sc.pl.dotplot(source, var_names=usable, groupby=group_key,
                       standard_scale="var", show=False, return_fig=True,
                       title="Fig 5D equivalent | canonical markers")
    fig = dp.get_axes()["mainplot_ax"].get_figure()
    return _save(fig, "fig5D_marker_dotplot")


def figure_5e_cluster_percentages(adata: AnnData, group_key: str = "cell_type") -> plt.Figure:
    """Paper Fig 5E: percentage of cells in each population, per timepoint.

    The paper pools cells within a condition. We additionally overlay the two
    replicate values, because pooling hides the fact that one replicate can supply
    twice the cells of the other.
    """
    pooled = (pd.crosstab(adata.obs[group_key], adata.obs["condition"],
                          normalize="columns") * 100)
    pooled = pooled[cfg.CONDITION_ORDER]

    per_sample = (
        adata.obs.groupby(["sample", "condition", group_key], observed=True)
        .size().reset_index(name="n")
    )
    per_sample["pct"] = 100 * per_sample["n"] / per_sample.groupby(
        "sample", observed=True)["n"].transform("sum")

    populations = list(pooled.index)
    x = np.arange(len(populations))
    width = 0.2
    fig, ax = plt.subplots(figsize=(max(9, 0.75 * len(populations)), 5))
    for i, cond in enumerate(cfg.CONDITION_ORDER):
        ax.bar(x + (i - 1.5) * width, pooled[cond].values, width,
               label=cfg.CONDITION_LABELS[cond],
               color=cfg.CONDITION_COLORS[cond], alpha=0.9)
        pts = per_sample[per_sample["condition"] == cond]
        for j, pop in enumerate(populations):
            vals = pts.loc[pts[group_key] == pop, "pct"]
            ax.scatter(np.full(len(vals), x[j] + (i - 1.5) * width), vals,
                       s=8, color="black", zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels(populations, rotation=45, ha="right")
    ax.set_ylabel("% of cells in condition")
    ax.set_title("Fig 5E equivalent | population percentage per timepoint\n"
                 "(bars = pooled; dots = individual replicates)")
    ax.legend(frameon=False, ncol=4)
    fig.tight_layout()
    return _save(fig, "fig5E_cluster_percentages")


# --------------------------------------------------------------------------- #
# Figure 6 -- rods
# --------------------------------------------------------------------------- #

PARALOG_PAIRS = [("rho", "rhol"), ("pde6ga", "pde6gb"), ("guca1a", "guca1b")]


def figure_6b_paralog_dotplot(rods: AnnData, group_key: str = "subcluster") -> plt.Figure:
    """Paper Fig 6B: the three paralog pairs that separate immature from mature rods.

    In the paper the two rod clusters show an INVERSE pattern: one is
    rho / pde6ga / guca1a high, the other rhol / pde6gb / guca1b high. Whether that
    inversion appears here is the decisive test of the immature/mature split, so
    this panel is printed with an explicit verdict rather than left to the eye.
    """
    genes = [g for pair in PARALOG_PAIRS for g in pair]
    present, missing = check_genes_available(
        rods.raw.to_adata() if rods.raw is not None else rods, genes, "Fig 6B paralogs")
    if missing:
        print(f"[Fig 6B] absent from this annotation: {missing}. "
              "The corresponding pair cannot be tested.")
    if not present:
        raise ValueError("None of the paralog genes are present; Fig 6B cannot be drawn.")

    df = _dense(rods, present)
    df[group_key] = rods.obs[group_key].values
    means = df.groupby(group_key, observed=True).mean()

    source = rods.raw.to_adata() if rods.raw is not None else rods
    source.obs = rods.obs.copy()
    dp = sc.pl.dotplot(source, var_names=present, groupby=group_key,
                       standard_scale="var", show=False, return_fig=True,
                       title="Fig 6B equivalent | rod paralog pairs")
    fig = dp.get_axes()["mainplot_ax"].get_figure()
    _save(fig, "fig6B_rod_paralog_dotplot")

    print("\nParalog inversion test (mean log-normalised expression):")
    print(means.round(3).to_string())
    verdict = []
    for a, b in PARALOG_PAIRS:
        if a in means and b in means:
            higher = (means[a] > means[b])
            if higher.all():
                verdict.append(f"  {a} > {b} in EVERY subcluster - no inversion")
            elif (~higher).all():
                verdict.append(f"  {b} > {a} in EVERY subcluster - no inversion")
            else:
                verdict.append(f"  {a}/{b} SPLIT across subclusters - inversion present")
    print("\n".join(verdict) if verdict else "  no testable pair available")
    print(
        "\nAn inversion requires the pairs to split in OPPOSITE directions across "
        "subclusters. If every pair points the same way in every subcluster, the "
        "paper's immature/mature rod distinction is not reproduced here, and that "
        "is the result to report."
    )
    return fig


def figure_6e_gene_heatmap(rods: AnnData, group_key: str = "subcluster",
                           genes: list[str] | None = None) -> plt.Figure:
    """Paper Fig 6E: heatmap of genes differing between rod populations."""
    genes = genes or ["rho", "rhol", "pde6ga", "pde6gb", "guca1a", "guca1b",
                      "meig1", "ppdpfa", "rom1b", "gnat1", "gngt1", "nr2e3",
                      "sagb", "rcvrna", "gnb1a", "crx", "prom1", "rbp4l"]
    df = _dense(rods, genes)
    if df.empty:
        raise ValueError("No genes from the rod panel are present.")
    df[group_key] = rods.obs[group_key].values
    means = df.groupby(group_key, observed=True).mean()
    z = (means - means.mean()) / means.std().replace(0, np.nan)

    fig, ax = plt.subplots(figsize=(0.45 * z.shape[1] + 3, 0.5 * z.shape[0] + 2.5))
    im = ax.imshow(z.values, cmap="viridis", aspect="auto")
    ax.set_xticks(range(z.shape[1]))
    ax.set_xticklabels(z.columns, rotation=90, style="italic")
    ax.set_yticks(range(z.shape[0]))
    ax.set_yticklabels([f"subcluster {i}" for i in z.index])
    ax.set_title("Fig 6E equivalent | rod subcluster expression (z-scored)")
    fig.colorbar(im, ax=ax, label="z-score of mean expression", shrink=0.7)
    fig.tight_layout()
    return _save(fig, "fig6E_rod_gene_heatmap")


# --------------------------------------------------------------------------- #
# Volcano (Fig 6F and 7E share this)
# --------------------------------------------------------------------------- #

def volcano(de_df: pd.DataFrame, title: str, filename: str,
            n_label: int = 12, lfc_cut: float = 0.25,
            padj_cut: float = 0.05) -> plt.Figure:
    """Volcano plot of a rank_genes_groups result, in the style of Fig 6F / 7E."""
    if de_df is None or len(de_df) == 0:
        raise ValueError(f"No differential expression rows supplied for '{title}'.")
    df = de_df.copy()
    df["padj"] = df["pval_adj"].replace(0, np.nextafter(0, 1))
    df["neglog10"] = -np.log10(df["padj"])
    df["direction"] = np.where(
        (df["pval_adj"] < padj_cut) & (df["logfoldchange"] >= lfc_cut), "up",
        np.where((df["pval_adj"] < padj_cut) & (df["logfoldchange"] <= -lfc_cut),
                 "down", "ns"))

    colors = {"up": "#B2182B", "down": "#2166AC", "ns": "#CCCCCC"}
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    for group, sub in df.groupby("direction"):
        ax.scatter(sub["logfoldchange"], sub["neglog10"], s=9,
                   c=colors[group], label=group, alpha=0.75, edgecolors="none")
    ax.axhline(-np.log10(padj_cut), ls="--", lw=0.8, color="#666666")
    ax.axvline(lfc_cut, ls="--", lw=0.8, color="#666666")
    ax.axvline(-lfc_cut, ls="--", lw=0.8, color="#666666")

    signif = df[df["direction"] != "ns"].nlargest(n_label, "neglog10")
    for _, row in signif.iterrows():
        ax.annotate(row["gene"], (row["logfoldchange"], row["neglog10"]),
                    fontsize=6.5, style="italic",
                    xytext=(3, 3), textcoords="offset points")

    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10 adjusted p-value")
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()

    n_up = int((df["direction"] == "up").sum())
    n_down = int((df["direction"] == "down").sum())
    print(f"{title}: {n_up} up, {n_down} down "
          f"(|log2FC| >= {lfc_cut}, padj < {padj_cut})")
    print("  Reminder: these p-values come from a cell-level test with n = 2 "
          "replicates per condition, so they overstate significance. Treat the "
          "ranking as descriptive.")
    return _save(fig, filename)


# --------------------------------------------------------------------------- #
# Figure 7 -- cones
# --------------------------------------------------------------------------- #

OPSIN_GENES = ["opn1sw1", "opn1sw2", "opn1mw1", "opn1mw2", "opn1mw3", "opn1mw4",
               "opn1lw1", "opn1lw2"]
CONE_GENES = ["arr3a", "arr3b", "spock3", "efna1b", "tgfa", "tbx2a",
              "cngb3.2", "guca1e", "kcnv2b", "gnat2", "pde6c"]


def figure_7b_opsin_heatmap(cones: AnnData, group_key: str = "subcluster") -> plt.Figure:
    """Paper Fig 7B: opsin-1 expression across cone subclusters (UV vs non-UV)."""
    df = _dense(cones, OPSIN_GENES)
    if df.empty:
        raise ValueError("No opsin genes present; Fig 7B cannot be drawn.")
    df[group_key] = cones.obs[group_key].values
    means = df.groupby(group_key, observed=True).mean()

    fig, ax = plt.subplots(figsize=(0.6 * means.shape[1] + 3, 0.5 * means.shape[0] + 2.5))
    im = ax.imshow(means.values, cmap="RdBu_r", aspect="auto")
    ax.set_xticks(range(means.shape[1]))
    ax.set_xticklabels(means.columns, rotation=90, style="italic")
    ax.set_yticks(range(means.shape[0]))
    ax.set_yticklabels([f"subcluster {i}" for i in means.index])
    ax.set_title("Fig 7B equivalent | opsin-1 genes across cone subclusters")
    fig.colorbar(im, ax=ax, label="mean log-normalised expression", shrink=0.7)
    fig.tight_layout()
    _save(fig, "fig7B_opsin_heatmap")

    if "opn1sw1" in means:
        top = means["opn1sw1"].idxmax()
        print(f"\nHighest opn1sw1 (UV opsin): subcluster {top}")
        print(means[["opn1sw1"] + [c for c in means.columns if c != "opn1sw1"]]
              .round(3).to_string())
        print("\nThe paper separates a UV cone cluster (opn1sw1-high, arr3b+, arr3a-) "
              "from non-UV cones (arr3a+). Check the arr3 paralogs in Fig 7D before "
              "concluding the subtypes are resolved.")
    return fig


def figure_7d_cone_heatmap(cones: AnnData, group_key: str = "subcluster") -> plt.Figure:
    """Paper Fig 7D: selected cone-subtype genes, including the arr3 paralogs."""
    df = _dense(cones, CONE_GENES)
    if df.empty:
        raise ValueError("No cone subtype genes present.")
    df[group_key] = cones.obs[group_key].values
    means = df.groupby(group_key, observed=True).mean()
    z = (means - means.mean()) / means.std().replace(0, np.nan)

    fig, ax = plt.subplots(figsize=(0.6 * z.shape[1] + 3, 0.5 * z.shape[0] + 2.5))
    im = ax.imshow(z.values, cmap="RdBu_r", aspect="auto")
    ax.set_xticks(range(z.shape[1]))
    ax.set_xticklabels(z.columns, rotation=90, style="italic")
    ax.set_yticks(range(z.shape[0]))
    ax.set_yticklabels([f"subcluster {i}" for i in z.index])
    ax.set_title("Fig 7D equivalent | cone subtype genes (z-scored)")
    fig.colorbar(im, ax=ax, label="z-score", shrink=0.7)
    fig.tight_layout()
    _save(fig, "fig7D_cone_gene_heatmap")

    if "arr3a" in means and "arr3b" in means:
        print("\narr3 paralog usage (the paper reports a near-inverse pattern):")
        print(means[["arr3a", "arr3b"]].round(3).to_string())
    return fig


# --------------------------------------------------------------------------- #
# Figure 8 -- EGFP
# --------------------------------------------------------------------------- #

def figure_8a_egfp_umap(adata: AnnData) -> plt.Figure:
    """Paper Fig 8A: EGFP-positive cells highlighted on the UMAP, per condition."""
    if "egfp_positive" not in adata.obs:
        raise KeyError("Run egfp_analysis.add_egfp_columns() first.")
    conditions = list(adata.obs["condition"].cat.categories)
    coords = adata.obsm["X_umap"]
    fig, axes = plt.subplots(1, len(conditions), figsize=(4.2 * len(conditions), 4.2))
    axes = np.atleast_1d(axes)
    for ax, cond in zip(axes, conditions):
        in_cond = (adata.obs["condition"] == cond).values
        pos = in_cond & adata.obs["egfp_positive"].values
        ax.scatter(coords[:, 0], coords[:, 1], s=1, c="#E5E5E5", rasterized=True)
        ax.scatter(coords[pos, 0], coords[pos, 1], s=6, c="#1B7837", rasterized=True)
        ax.set_title(f"{cfg.CONDITION_LABELS.get(str(cond), cond)}\n"
                     f"{int(pos.sum())} EGFP+ / {int(in_cond.sum())} cells "
                     f"({100 * pos.sum() / in_cond.sum():.2f}%)", fontsize=9)
        ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle("Fig 8A equivalent | EGFP-positive cells per condition", fontsize=11)
    fig.tight_layout()
    return _save(fig, "fig8A_egfp_umap_by_condition")


def figure_8b_egfp_counts(adata: AnnData) -> plt.Figure:
    """Paper Fig 8B: number of EGFP-positive cells per timepoint."""
    counts = (adata.obs.groupby("condition", observed=True)["egfp_positive"]
              .agg(["sum", "size"]).reindex(cfg.CONDITION_ORDER))
    per_sample = (adata.obs.groupby(["condition", "sample"], observed=True)
                  ["egfp_positive"].sum().reset_index())

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    x = np.arange(len(cfg.CONDITION_ORDER))
    axes[0].bar(x, counts["sum"].values,
                color=[cfg.CONDITION_COLORS[c] for c in cfg.CONDITION_ORDER])
    for i, cond in enumerate(cfg.CONDITION_ORDER):
        vals = per_sample.loc[per_sample["condition"] == cond, "egfp_positive"]
        axes[0].scatter(np.full(len(vals), i), vals, color="black", s=20, zorder=3)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels([cfg.CONDITION_LABELS[c] for c in cfg.CONDITION_ORDER])
    axes[0].set_ylabel("number of EGFP-positive cells")
    axes[0].set_title("A  EGFP+ cell count")

    pct = 100 * counts["sum"] / counts["size"]
    axes[1].bar(x, pct.values,
                color=[cfg.CONDITION_COLORS[c] for c in cfg.CONDITION_ORDER])
    axes[1].set_xticks(x)
    axes[1].set_xticklabels([cfg.CONDITION_LABELS[c] for c in cfg.CONDITION_ORDER])
    axes[1].set_ylabel("% EGFP-positive")
    axes[1].set_title("B  EGFP+ percentage")
    for i, v in enumerate(pct.values):
        axes[1].text(i, v, f"{v:.2f}%", ha="center", va="bottom", fontsize=8)

    fig.suptitle("Fig 8B equivalent | EGFP-positive cells across the time course",
                 fontsize=11)
    fig.tight_layout()
    _save(fig, "fig8B_egfp_counts")
    print("\nPaper values for comparison: 0.54% (ctrl), 5.40% (3dp), "
          "2.64% (7dp), 2.66% (10dp).")
    return fig


def figure_8c_egfp_per_cluster(adata: AnnData, group_key: str = "cell_type") -> plt.Figure:
    """Paper Fig 8C: proportion of EGFP-positive cells in each population."""
    grouped = adata.obs.groupby(group_key, observed=True)["egfp_positive"]
    df = pd.DataFrame({"n": grouped.size(), "pos": grouped.sum()})
    df["pct"] = 100 * df["pos"] / df["n"]
    df = df.sort_values("pct", ascending=False)

    fig, ax = plt.subplots(figsize=(max(8, 0.6 * len(df)), 4.5))
    ax.bar(df.index.astype(str), df["pct"],
           color=[cfg.CELL_TYPE_COLORS.get(str(c), "#999999") for c in df.index])
    for i, (n, pos) in enumerate(zip(df["n"], df["pos"])):
        ax.text(i, df["pct"].iloc[i], f"{int(pos)}/{int(n)}",
                ha="center", va="bottom", fontsize=6)
    ax.set_ylabel("% EGFP-positive cells")
    ax.set_xticklabels(df.index.astype(str), rotation=45, ha="right")
    ax.set_title("Fig 8C equivalent | EGFP-positive proportion per population\n"
                 "(labels show positive/total; small denominators are unreliable)")
    fig.tight_layout()
    _save(fig, "fig8C_egfp_per_cluster")
    print("\nPaper: EGFP transcripts in nearly 10% of Muller glia, an outstanding "
          "proportion among other cell types.")
    return fig


def figure_8e_egfp_mg_heatmap(adata: AnnData, mg_de: pd.DataFrame,
                              group_key: str = "cell_type",
                              n_genes: int = 25) -> plt.Figure:
    """Paper Fig 8E: genes differing between EGFP+ and EGFP- Muller glia.

    The paper reports 42 upregulated genes. If the comparison here was underpowered
    the honest output is no panel at all, so this raises rather than drawing an
    empty heatmap.
    """
    if mg_de is None or len(mg_de) == 0:
        raise ValueError(
            "No EGFP+/- Muller glia DE result. If egfp_de_within_group() reported "
            "too few cells on one side, report the comparison as underpowered "
            "rather than producing this panel."
        )
    top = mg_de.nlargest(n_genes, "score")["gene"].tolist()
    mg_types = [c for c in adata.obs[group_key].cat.categories
                if "muller" in str(c).lower()]
    mg = adata[adata.obs[group_key].isin(mg_types)]

    df = _dense(mg, top)
    if df.empty:
        raise ValueError("None of the top DE genes are retrievable from .raw.")
    df["egfp"] = np.where(mg.obs["egfp_positive"].values, "EGFP+", "EGFP-")
    means = df.groupby("egfp", observed=True).mean().reindex(["EGFP-", "EGFP+"])
    z = (means - means.mean()) / means.std().replace(0, np.nan)

    fig, ax = plt.subplots(figsize=(0.42 * z.shape[1] + 3, 3))
    im = ax.imshow(z.values, cmap="RdBu_r", aspect="auto")
    ax.set_xticks(range(z.shape[1]))
    ax.set_xticklabels(z.columns, rotation=90, style="italic", fontsize=7)
    ax.set_yticks(range(z.shape[0]))
    ax.set_yticklabels(z.index)
    n_pos = int(mg.obs["egfp_positive"].sum())
    ax.set_title(f"Fig 8E equivalent | EGFP+ vs EGFP- Muller glia\n"
                 f"({n_pos} EGFP+ of {mg.n_obs} MG cells)")
    fig.colorbar(im, ax=ax, label="z-score", shrink=0.8)
    fig.tight_layout()
    _save(fig, "fig8E_egfp_mg_heatmap")
    print(f"\nPaper: 42 genes upregulated in careg:EGFP-positive Muller glia, "
          f"including six3b, glulb, gfap, apoeb, mdka, crabp1a, id1, mmp9.")
    return fig


# --------------------------------------------------------------------------- #
# GO enrichment (Fig 6C/D and 7F) -- optional, needs network access
# --------------------------------------------------------------------------- #

def go_enrichment_bar(genes: list[str], title: str, filename: str,
                      gene_sets: str = "GO_Biological_Process_2018",
                      top_n: int = 8) -> plt.Figure | None:
    """Paper Fig 6C/D and 7F: GO term bar chart, -log10 Fisher p-value.

    Uses gseapy's Enrichr client with organism='Fish'. Requires network access; if
    unavailable this returns None and says so rather than failing the notebook.

    Not equivalent to the paper's method: they used topGO with a Fisher exact test
    against a custom background. Enrichr uses its own background and gene-set
    versions, so term lists will differ even on identical input. State this in
    Methods.
    """
    if not genes:
        print(f"[{title}] no input genes; skipped.")
        return None
    try:
        import gseapy as gp
    except ImportError:
        print("[GO] gseapy is not installed; skipping enrichment.")
        return None

    try:
        res = gp.enrichr(gene_list=list(genes), gene_sets=gene_sets,
                         organism="Fish", outdir=None)
        table = res.results
    except Exception as exc:  # network, API change, empty result
        print(f"[GO] Enrichment unavailable ({type(exc).__name__}: {exc}). "
              "Report GO analysis as not performed rather than omitting silently.")
        return None

    if table is None or len(table) == 0:
        print(f"[{title}] Enrichr returned no terms.")
        return None

    table = table.nsmallest(top_n, "P-value").iloc[::-1]
    fig, ax = plt.subplots(figsize=(7.5, 0.45 * len(table) + 2))
    ax.barh(table["Term"].str.slice(0, 55), -np.log10(table["P-value"]),
            color="#F4C430", edgecolor="#8B6914")
    ax.set_xlabel("-log10 p-value")
    ax.set_title(title)
    fig.tight_layout()
    return _save(fig, filename)


# --------------------------------------------------------------------------- #
# Highlight panels: Fig 6A (rods), 7A (cones), 8D (Muller glia)
# --------------------------------------------------------------------------- #

def highlight_umap(adata: AnnData, cell_types: list[str], title: str,
                   filename: str, group_key: str = "cell_type",
                   color: str = "#B2182B") -> plt.Figure:
    """One population in colour, the rest in grey — the paper's 6A / 7A / 8D style.

    Bise et al. use this layout to show where a subset sits in the whole atlas
    before zooming into it. Purely a locator panel: it carries no quantity that can
    be compared numerically, so it supports the figure rather than the argument.
    """
    present = [c for c in cell_types if c in set(adata.obs[group_key].astype(str))]
    missing = [c for c in cell_types if c not in present]
    if missing:
        print(f"[{title}] not among the annotated populations: {missing}")
    if not present:
        raise ValueError(f"None of {cell_types} are present in obs['{group_key}'].")

    mask = adata.obs[group_key].astype(str).isin(present).values
    coords = adata.obsm["X_umap"]

    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    ax.scatter(coords[~mask, 0], coords[~mask, 1], s=2, c="#DDDDDD",
               rasterized=True, label="other cells")
    ax.scatter(coords[mask, 0], coords[mask, 1], s=4, c=color,
               rasterized=True, label=", ".join(present))
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
    ax.set_title(f"{title}\n{int(mask.sum())} of {adata.n_obs} cells "
                 f"({100 * mask.sum() / adata.n_obs:.1f}%)")
    ax.legend(frameon=False, fontsize=7, markerscale=2.5, loc="best")
    fig.tight_layout()
    return _save(fig, filename)


def figure_6a_rod_umap(adata: AnnData, group_key: str = "cell_type") -> plt.Figure:
    """Paper Fig 6A: rod clusters in colour, all other cells grey."""
    rod_types = [str(c) for c in adata.obs[group_key].cat.categories
                 if "rod" in str(c).lower()]
    return highlight_umap(adata, rod_types, "Fig 6A equivalent | rod populations",
                          "fig6A_rod_umap", group_key, color=cfg.CELL_TYPE_COLORS.get(
                              "Rods", "#0072B2"))


def figure_7a_cone_umap(adata: AnnData, group_key: str = "cell_type") -> plt.Figure:
    """Paper Fig 7A: cone clusters in colour, all other cells grey."""
    cone_types = [str(c) for c in adata.obs[group_key].cat.categories
                  if "cone" in str(c).lower()]
    return highlight_umap(adata, cone_types, "Fig 7A equivalent | cone populations",
                          "fig7A_cone_umap", group_key, color=cfg.CELL_TYPE_COLORS.get(
                              "Cones", "#009E73"))


def figure_8d_mg_umap(adata: AnnData, group_key: str = "cell_type") -> plt.Figure:
    """Paper Fig 8D: the Muller glia cluster highlighted in the integrated data."""
    mg_types = [str(c) for c in adata.obs[group_key].cat.categories
                if "muller" in str(c).lower()]
    return highlight_umap(adata, mg_types, "Fig 8D equivalent | Muller glia cluster",
                          "fig8D_muller_glia_umap", group_key,
                          color=cfg.CELL_TYPE_COLORS.get("Muller glia", "#E69F00"))


NOT_REPRODUCIBLE = {
    "Fig 5A": "experimental design schematic (MNU timeline, dark adaptation) - a drawing, not data",
    "Fig 5B": "immunofluorescence: careg:EGFP absent after inactivated MNU - microscopy",
    "Fig 7C": "colour scale bar for panels 7B and 7D - drawn as a colorbar on each heatmap here",
    "Fig 1-4": "immunofluorescence of retinal sections",
    "Fig 9-11": "TOR signalling (p-rpS6) and rapamycin experiments - protein and drug work",
}
