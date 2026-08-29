"""Loading, validation, metadata and checkpoint helpers.

Design rules enforced here:
  * sample directories are discovered, never hard-coded to an absolute path;
  * every expected 10x file is checked before reading, with a readable error;
  * the EGFP feature name is detected from the data, never assumed;
  * `condition` is a pandas Categorical in biological time order from the moment
    it is created, so no downstream plot can silently reorder the time course.
"""

from __future__ import annotations

import json
import platform
from pathlib import Path
from typing import Iterable

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from . import config as cfg


# --------------------------------------------------------------------------- #
# Discovery and validation
# --------------------------------------------------------------------------- #

def discover_samples(data_dir: Path | None = None) -> list[Path]:
    """Return sample directories that contain a 10x matrix folder.

    Raises
    ------
    FileNotFoundError
        If `data_dir` does not exist or contains no usable sample directory.
    """
    data_dir = Path(data_dir) if data_dir is not None else cfg.DATA_DIR
    if not data_dir.exists():
        raise FileNotFoundError(
            f"Data directory not found: {data_dir}\n"
            "Run the notebooks from the repository root, and place the eight 10x "
            "sample folders under data/ as described in README.md."
        )

    sample_dirs = sorted(
        p for p in data_dir.iterdir()
        if p.is_dir() and (p / cfg.MATRIX_SUBDIR).is_dir()
    )
    if not sample_dirs:
        raise FileNotFoundError(
            f"No sample directories with a '{cfg.MATRIX_SUBDIR}/' subfolder were found "
            f"in {data_dir}. Found instead: {[p.name for p in data_dir.iterdir()]}"
        )
    return sample_dirs


def validate_sample_dir(sample_dir: Path) -> None:
    """Check that a sample directory holds all three 10x files."""
    matrix_dir = sample_dir / cfg.MATRIX_SUBDIR
    missing = [f for f in cfg.REQUIRED_10X_FILES if not (matrix_dir / f).exists()]
    if missing:
        raise FileNotFoundError(
            f"Sample '{sample_dir.name}' is missing {missing} in {matrix_dir}. "
            "Cell Ranger writes barcodes.tsv.gz, features.tsv.gz and matrix.mtx.gz; "
            "check whether the folder was unpacked or renamed."
        )


def report_expected_samples(sample_dirs: Iterable[Path]) -> pd.DataFrame:
    """Compare discovered sample names against the eight expected ones."""
    found = {p.name for p in sample_dirs}
    rows = [
        {"sample": s, "expected": True, "found": s in found}
        for s in cfg.EXPECTED_SAMPLES
    ]
    rows += [
        {"sample": s, "expected": False, "found": True}
        for s in sorted(found - set(cfg.EXPECTED_SAMPLES))
    ]
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Metadata
# --------------------------------------------------------------------------- #

def parse_sample_name(sample_name: str) -> dict[str, str]:
    """Split a sample folder name into condition and replicate.

    'ctrl1' -> condition 'ctrl', replicate '1'
    '10dp2' -> condition '10dp', replicate '2'
    """
    if sample_name.startswith("ctrl"):
        condition, replicate = "ctrl", sample_name[len("ctrl"):] or "1"
    elif "dp" in sample_name:
        days, replicate = sample_name.split("dp", 1)
        condition, replicate = f"{days}dp", replicate or "1"
    else:
        raise ValueError(
            f"Cannot parse sample name '{sample_name}'. Expected 'ctrl<rep>' or "
            "'<days>dp<rep>', e.g. 'ctrl1' or '10dp2'."
        )

    if condition not in cfg.CONDITION_ORDER:
        raise ValueError(
            f"Sample '{sample_name}' gives condition '{condition}', which is not in "
            f"CONDITION_ORDER {cfg.CONDITION_ORDER}."
        )
    return {"condition": condition, "replicate": replicate}


def attach_metadata(adata: ad.AnnData, sample_name: str) -> ad.AnnData:
    """Add sample / condition / replicate / timepoint / batch columns in place."""
    parsed = parse_sample_name(sample_name)
    adata.obs["sample"] = sample_name
    adata.obs["condition"] = parsed["condition"]
    adata.obs["replicate"] = parsed["replicate"]
    adata.obs["timepoint"] = cfg.CONDITION_LABELS[parsed["condition"]]
    adata.obs["days_post_injury"] = cfg.CONDITION_DAYS[parsed["condition"]]
    # One 10x run per sample, so sample and batch coincide here. Kept as a separate
    # column because "batch" is what integration tools consume, and the two would
    # differ if samples were ever multiplexed.
    adata.obs["batch"] = sample_name
    return adata


def set_categorical_order(adata: ad.AnnData) -> ad.AnnData:
    """Force biological ordering on the categorical obs columns."""
    adata.obs["condition"] = pd.Categorical(
        adata.obs["condition"], categories=cfg.CONDITION_ORDER, ordered=True
    )
    adata.obs["timepoint"] = pd.Categorical(
        adata.obs["timepoint"],
        categories=[cfg.CONDITION_LABELS[c] for c in cfg.CONDITION_ORDER],
        ordered=True,
    )
    sample_order = [s for s in cfg.EXPECTED_SAMPLES if s in set(adata.obs["sample"])]
    sample_order += sorted(set(adata.obs["sample"]) - set(sample_order))
    adata.obs["sample"] = pd.Categorical(
        adata.obs["sample"], categories=sample_order, ordered=True
    )
    return adata


# --------------------------------------------------------------------------- #
# Loading
# --------------------------------------------------------------------------- #

def load_sample(sample_dir: Path) -> ad.AnnData:
    """Read one 10x filtered matrix and attach its metadata."""
    validate_sample_dir(sample_dir)
    adata = sc.read_10x_mtx(
        sample_dir / cfg.MATRIX_SUBDIR, var_names="gene_symbols", cache=False
    )
    adata.var_names_make_unique()
    attach_metadata(adata, sample_dir.name)
    adata.obs_names = [f"{sample_dir.name}_{bc}" for bc in adata.obs_names]
    return adata


def load_all_samples(data_dir: Path | None = None) -> tuple[ad.AnnData, pd.DataFrame]:
    """Load every discovered sample and concatenate into one AnnData object.

    Returns
    -------
    adata : AnnData
        Concatenated object with unique obs_names and ordered categoricals.
    per_sample : DataFrame
        One row per sample: cells, genes, total counts, median counts/genes and
        whether an EGFP-like feature is present in that sample's feature table.
    """
    sample_dirs = discover_samples(data_dir)
    print(f"Discovered {len(sample_dirs)} samples: {[p.name for p in sample_dirs]}")

    adatas: dict[str, ad.AnnData] = {}
    rows: list[dict] = []

    for sample_dir in sample_dirs:
        a = load_sample(sample_dir)
        counts_per_cell = np.asarray(a.X.sum(axis=1)).ravel()
        genes_per_cell = np.asarray((a.X > 0).sum(axis=1)).ravel()
        egfp = find_egfp_feature(a, verbose=False)
        rows.append({
            "sample": sample_dir.name,
            "condition": a.obs["condition"].iloc[0],
            "replicate": a.obs["replicate"].iloc[0],
            "n_cells": a.n_obs,
            "n_genes": a.n_vars,
            "total_counts": float(counts_per_cell.sum()),
            "median_counts_per_cell": float(np.median(counts_per_cell)),
            "median_genes_per_cell": float(np.median(genes_per_cell)),
            "egfp_feature": egfp if egfp is not None else "ABSENT",
        })
        adatas[sample_dir.name] = a
        print(f"  {sample_dir.name:<7} {a.n_obs:>6} cells x {a.n_vars:>6} genes"
              f"   EGFP feature: {egfp if egfp else 'ABSENT'}")

    adata = ad.concat(adatas, label="sample_batch", index_unique=None, join="outer")
    adata.var_names_make_unique()
    set_categorical_order(adata)

    assert adata.n_obs > 0, "Concatenated object has no cells."
    assert adata.obs_names.is_unique, "Cell barcodes are not unique after concatenation."
    for col in ("sample", "condition", "replicate", "timepoint", "batch"):
        assert col in adata.obs, f"Metadata column '{col}' missing after concatenation."

    print(f"\nTotal dataset: {adata.n_obs} cells x {adata.n_vars} genes")
    return adata, pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# EGFP feature detection
# --------------------------------------------------------------------------- #

def find_egfp_feature(adata: ad.AnnData, verbose: bool = True) -> str | None:
    """Locate the transgene feature in var_names without assuming its spelling.

    Tries the configured candidates case-insensitively, then falls back to any
    feature whose name contains 'gfp' or 'careg'. Returns the exact var_name as it
    appears in the object, or None.
    """
    lookup = {name.lower(): name for name in adata.var_names}

    for candidate in cfg.EGFP_NAME_CANDIDATES:
        hit = lookup.get(candidate.lower())
        if hit is not None:
            if verbose:
                print(f"EGFP feature found in var_names as: '{hit}'")
            return hit

    fuzzy = [n for n in adata.var_names if "gfp" in n.lower() or "careg" in n.lower()]
    if len(fuzzy) == 1:
        if verbose:
            print(f"EGFP feature found by substring match: '{fuzzy[0]}'")
        return fuzzy[0]
    if len(fuzzy) > 1:
        raise ValueError(
            f"Several features look like the transgene: {fuzzy}. Inspect "
            "features.tsv.gz and set the correct name explicitly rather than guessing."
        )

    if verbose:
        print(
            "No EGFP-like feature in var_names. The transgene was probably not added "
            "to the Cell Ranger reference for this matrix. Every EGFP-dependent "
            "result (Figure 4, the careg analysis) is then unavailable and must be "
            "reported as such, not approximated with another gene."
        )
    return None


def require_egfp(adata: ad.AnnData) -> str:
    """Return the EGFP feature name or raise with an actionable message."""
    name = find_egfp_feature(adata, verbose=False)
    if name is None:
        raise KeyError(
            "EGFP feature not present in adata.var_names. Check features.tsv.gz for "
            "the transgene entry; if it is genuinely absent, the reporter analysis "
            "cannot be run on this matrix."
        )
    return name


def check_genes_available(
    adata: ad.AnnData, genes: Iterable[str], label: str = "panel"
) -> tuple[list[str], list[str]]:
    """Split a gene list into present / absent. Absent genes are reported, never replaced."""
    present = [g for g in genes if g in adata.var_names]
    missing = [g for g in genes if g not in adata.var_names]
    if missing:
        print(f"[{label}] not in this dataset ({len(missing)}): {missing}")
    return present, missing


# --------------------------------------------------------------------------- #
# Checkpoints, seeds, provenance
# --------------------------------------------------------------------------- #

def set_seeds(seed: int | None = None) -> int:
    """Seed every RNG the workflow touches. Returns the seed used."""
    import random
    seed = cfg.RANDOM_SEED if seed is None else seed
    random.seed(seed)
    np.random.seed(seed)
    return seed


def save_checkpoint(adata: ad.AnnData, key: str) -> Path:
    """Write an intermediate AnnData object to results/ and return its path."""
    if key not in cfg.H5AD:
        raise KeyError(f"Unknown checkpoint '{key}'. Known: {sorted(cfg.H5AD)}")
    path = cfg.H5AD[key]
    path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(path, compression="gzip")
    size_mb = path.stat().st_size / 1e6
    print(f"Saved checkpoint '{key}' -> {path.relative_to(cfg.REPO_ROOT)} ({size_mb:.1f} MB)")
    return path


def load_checkpoint(key: str) -> ad.AnnData:
    """Read an intermediate AnnData object, restoring categorical ordering."""
    if key not in cfg.H5AD:
        raise KeyError(f"Unknown checkpoint '{key}'. Known: {sorted(cfg.H5AD)}")
    path = cfg.H5AD[key]
    if not path.exists():
        raise FileNotFoundError(
            f"Checkpoint '{key}' not found at {path}. Run the earlier notebook first; "
            "results/ is git-ignored, so a fresh clone starts with none of these files."
        )
    adata = sc.read_h5ad(path)
    set_categorical_order(adata)
    print(f"Loaded '{key}': {adata.n_obs} cells x {adata.n_vars} genes")
    return adata


def capture_versions(extra: dict | None = None) -> dict:
    """Collect package versions for the reproducibility record."""
    import scipy
    import matplotlib

    versions = {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "scanpy": sc.__version__,
        "anndata": ad.__version__,
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "matplotlib": matplotlib.__version__,
        "random_seed": cfg.RANDOM_SEED,
    }
    try:
        import leidenalg
        versions["leidenalg"] = leidenalg.version
    except Exception:
        versions["leidenalg"] = "not installed"
    if extra:
        versions.update(extra)
    return versions


def write_versions(path: Path | None = None, extra: dict | None = None) -> Path:
    """Write the version record to docs/session_info.json."""
    path = path or (cfg.DOCS_DIR / "session_info.json")
    path.parent.mkdir(parents=True, exist_ok=True)
    versions = capture_versions(extra)
    path.write_text(json.dumps(versions, indent=2))
    for k, v in versions.items():
        print(f"  {k:<12} {v}")
    return path


def save_table(df: pd.DataFrame, name: str, index: bool = False) -> Path:
    """Write a DataFrame to tables/<name> and return the path."""
    cfg.TABLES_DIR.mkdir(parents=True, exist_ok=True)
    path = cfg.TABLES_DIR / name
    df.to_csv(path, index=index)
    print(f"Wrote {path.relative_to(cfg.REPO_ROOT)}  ({len(df)} rows)")
    return path
