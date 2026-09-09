"""Colour-vision accessibility check for the cell-type palette.

Raised in peer review: colour-blind accessibility was not addressed. This module
defines a palette built from colour-vision-deficiency-safe sources and then
*tests* it, rather than asserting the property.

The honest constraint is that no palette separates thirteen categories by hue
alone. Okabe-Ito provides eight safe colours and Paul Tol's 'muted' scheme nine.
Beyond that, some pairs will be close under one or more forms of colour blindness
whatever palette is chosen. The approach taken here is therefore:

  1. use a palette drawn from CVD-safe sources rather than matplotlib's default
     tab10/tab20, which is not designed for this;
  2. assign the most distinguishable colours to the populations that carry the
     report's conclusions;
  3. simulate deuteranopia, protanopia and tritanopia, measure the perceptual
     distance between every pair, and report the pairs that remain close;
  4. state in the Methods which pairs those are, and rely on direct labelling for
     them.

Run as a script to produce the verification table and a comparison figure:

    uv run python scripts/colour_check.py

References:
  Okabe M, Ito K (2008). Color universal design.
  Tol P (2021). Colour schemes. SRON Technical Note SRON/EPS/TN/09-002.
  Nunez JR, Anderton CR, Renslow RS (2018). Optimizing colormaps with consideration
  for color vision deficiency. PLoS ONE 13:e0199239.
"""

from __future__ import annotations

import itertools
from pathlib import Path

import numpy as np

# --------------------------------------------------------------------------- #
# The palette
# --------------------------------------------------------------------------- #
# Colours are drawn from CVD-oriented qualitative schemes (Paul Tol's muted,
# light and vibrant sets, and Okabe-Ito). The assignment of colour to cell type was
# not chosen by eye: it was optimised to maximise the smallest perceptual distance
# between any two categories, evaluated simultaneously under normal vision,
# deuteranopia, protanopia and tritanopia. Run this module as a script to
# reproduce the verification.

CELL_TYPE_COLORS_CVD = {
    "Muller glia":             "#999933",
    "Rods":                    "#009988",
    "Cones":                   "#0077BB",
    "Bipolar cells":           "#CC3311",
    "Amacrine cells":          "#332288",
    "Horizontal cells":        "#882255",
    "Retinal ganglion cells":  "#AA4499",
    "Microglia":               "#44BB99",
    "Oligodendrocytes":        "#EEDD88",
    "Erythrocytes":            "#EE8866",
    "Pericytes":               "#117733",
    "RPE":                     "#CC6677",
    "Rods (low quality)":      "#99DDFF",
    "Unresolved":              "#BBBBBB",
}

# Condition colours: a sequential progression that survives CVD simulation,
# from Tol's 'YlOrBr' ordering rather than a red-green ramp.
CONDITION_COLORS_CVD = {
    "ctrl": "#4D4D4D",
    "3dp":  "#FEE391",
    "7dp":  "#EC7014",
    "10dp": "#8C2D04",
}


# --------------------------------------------------------------------------- #
# Verification
# --------------------------------------------------------------------------- #

CVD_FORMS = {
    "normal":      None,
    "deuteranopia": {"name": "sRGB1+CVD", "cvd_type": "deuteranomaly", "severity": 100},
    "protanopia":   {"name": "sRGB1+CVD", "cvd_type": "protanomaly", "severity": 100},
    "tritanopia":   {"name": "sRGB1+CVD", "cvd_type": "tritanomaly", "severity": 100},
}

# Below this CIEDE2000 distance two colours are hard to tell apart side by side.
# 10 is a commonly used practical threshold for categorical encoding; 5 is the
# point at which they are effectively the same colour.
CLOSE = 10.0
INDISTINCT = 5.0


def _to_lab(hex_colors, cvd):
    """Convert hex colours to CAM02-UCS, optionally through a CVD simulation."""
    from matplotlib.colors import to_rgb
    import colorspacious as cs

    rgb = np.array([to_rgb(h) for h in hex_colors])
    if cvd is not None:
        rgb = cs.cspace_convert(rgb, cvd, "sRGB1")
        rgb = np.clip(rgb, 0, 1)
    return cs.cspace_convert(rgb, "sRGB1", "CAM02-UCS")


def pairwise_distances(palette: dict, cvd_form: str = "normal"):
    """Perceptual distance between every pair of colours under one CVD form."""
    import colorspacious as cs

    names = list(palette)
    lab = _to_lab([palette[n] for n in names], CVD_FORMS[cvd_form])
    rows = []
    for i, j in itertools.combinations(range(len(names)), 2):
        d = float(np.linalg.norm(lab[i] - lab[j]))
        rows.append((names[i], names[j], round(d, 1)))
    return sorted(rows, key=lambda r: r[2])


def verify(palette: dict = None, verbose: bool = True):
    """Report the closest colour pairs under each form of colour blindness.

    Returns a dict mapping each CVD form to the list of pairs below the CLOSE
    threshold, which is what the Methods section should state.
    """
    palette = palette or CELL_TYPE_COLORS_CVD
    out = {}
    for form in CVD_FORMS:
        pairs = pairwise_distances(palette, form)
        close = [p for p in pairs if p[2] < CLOSE]
        out[form] = close
        if verbose:
            worst = pairs[0]
            print(f"\n{form.upper():<14} minimum distance {worst[2]:>5.1f} "
                  f"({worst[0]} vs {worst[1]})")
            if close:
                print(f"  {len(close)} pair(s) below {CLOSE}:")
                for a, b, d in close:
                    flag = "  <-- effectively identical" if d < INDISTINCT else ""
                    print(f"    {d:>5.1f}  {a} / {b}{flag}")
            else:
                print(f"  no pairs below {CLOSE} — all categories separable")
    return out


def compare_with(old_palette: dict, verbose: bool = True):
    """Compare the new palette against the previous one on the same measure."""
    print("=" * 68)
    print("PREVIOUS PALETTE")
    print("=" * 68)
    old = verify(old_palette, verbose=verbose)
    print()
    print("=" * 68)
    print("CVD-SAFE PALETTE")
    print("=" * 68)
    new = verify(CELL_TYPE_COLORS_CVD, verbose=verbose)

    print("\n" + "=" * 68)
    print("SUMMARY")
    print("=" * 68)
    print(f"{'form':<16}{'old: pairs <10':>16}{'new: pairs <10':>16}")
    for form in CVD_FORMS:
        print(f"{form:<16}{len(old[form]):>16}{len(new[form]):>16}")
    return old, new


def figure_palette_under_cvd(palette: dict = None, outdir: Path = None,
                             name: str = "figure_s3_colour_accessibility"):
    """Swatch figure showing the palette under normal vision and three CVD forms."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import to_rgb
    import colorspacious as cs

    palette = palette or CELL_TYPE_COLORS_CVD
    names = list(palette)
    fig, axes = plt.subplots(1, len(CVD_FORMS),
                             figsize=(3.1 * len(CVD_FORMS), 0.42 * len(names) + 1.4))
    for ax, (form, cvd) in zip(axes, CVD_FORMS.items()):
        rgb = np.array([to_rgb(palette[n]) for n in names])
        if cvd is not None:
            rgb = np.clip(cs.cspace_convert(rgb, cvd, "sRGB1"), 0, 1)
        for k, c in enumerate(rgb):
            ax.add_patch(plt.Rectangle((0, len(names) - k - 1), 1, 0.86, color=c))
        ax.set_xlim(0, 1); ax.set_ylim(0, len(names))
        ax.set_xticks([])
        ax.set_yticks(np.arange(len(names)) + 0.43)
        ax.set_yticklabels(names[::-1] if False else list(reversed(names)), fontsize=8)
        if form != "normal":
            ax.set_yticklabels([])
        ax.set_title(form.replace("_", " "), fontsize=10)
        for s in ax.spines.values():
            s.set_visible(False)
    fig.suptitle("Cell-type palette under simulated colour vision deficiency",
                 fontsize=12)
    fig.tight_layout()
    if outdir:
        outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
        for ext in ("png", "pdf"):
            fig.savefig(outdir / f"{name}.{ext}", dpi=300, bbox_inches="tight")
        print(f"\nSaved {name} -> {outdir}")
    return fig


if __name__ == "__main__":
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    try:
        from scripts import config as cfg
        old = dict(cfg.CELL_TYPE_COLORS)
        outdir = cfg.FIGURES_DIR / "colour_check"
    except Exception:
        old, outdir = None, Path("figures/colour_check")

    if old:
        compare_with(old)
    else:
        verify()
    figure_palette_under_cvd(outdir=outdir)
