# Validation checklist

Work through this before tagging `v1.0-peerreview`. Anything unticked is either
fixed or stated as a limitation in the report — an acknowledged gap is defensible,
a silent one is not.

## Environment
- [ ] Kernel is "Python 3 (retina)", not system Python
- [ ] `docs/session_info.json` written and committed
- [ ] `requirements.txt` re-pinned with `uv pip freeze` after the final install
- [ ] Fresh clone + `uv pip install -r requirements.txt` reproduces the environment

## Data loading
- [ ] All 8 samples discovered; any missing one recorded
- [ ] Barcodes unique; no cell lost to a name collision
- [ ] `condition` prints as `['ctrl', '3dp', '7dp', '10dp']` — not alphabetical
- [ ] EGFP feature presence and exact spelling confirmed against `features.tsv.gz`
- [ ] Control EGFP rate is low (a high control rate means investigate, not proceed)

## QC
- [ ] `pct_counts_mt` is not identically zero
- [ ] Threshold diagnostics inspected **before** filtering
- [ ] Per-sample cell loss recorded; condition-dependent loss noted for the Discussion
- [ ] Any deviation from baseline thresholds is in `cfg.PARAM_OVERRIDES` with evidence

## Preprocessing
- [ ] `layers['counts']` still integer
- [ ] `.raw` covers all genes, not only HVGs
- [ ] Scaled `X` has mean ≈ 0, sd ≈ 1
- [ ] You can say in one sentence which matrix DE uses, and why

## Clustering
- [ ] Elbow plot inspected; `n_pcs` justified
- [ ] Batch diagnostics run; integrate / do-not-integrate decision written down
- [ ] No cluster >90% one sample, or those that are have been explained
- [ ] Resolution sweep reported

## Annotation
- [ ] Every cluster has an entry with a non-empty rationale
- [ ] No cell type called from a single marker
- [ ] Müller glia identified by identity markers, independently of EGFP
- [ ] Ambiguous clusters labelled `Unresolved`, not guessed
- [ ] RGC outcome recorded (found, or absent as in the paper)
- [ ] Every marker gene verified to exist in the zebrafish annotation

## EGFP
- [ ] Positivity threshold justified from the observed count distribution
- [ ] EGFP⁺ fraction reported per replicate, not only pooled
- [ ] EGFP enrichment by cell type reported with raw counts alongside ratios
- [ ] Any EGFP⁺ cells outside Müller glia examined rather than ignored

## Photoreceptors
- [ ] Rod paralog inversion (`rho`/`rhol`, `pde6ga`/`pde6gb`, `guca1a`/`guca1b`) either reproduced or reported as not reproduced
- [ ] Rod subclusters checked against counts/genes/mito before being called biological
- [ ] No recovery claim anywhere — the series stops at day 10

## Figures and tables
- [ ] Figures 1–5 exist as both PNG and PDF
- [ ] Cell-type colours identical across Figures 2, 3, 4 and 5
- [ ] Timepoints ordered ctrl → 3dp → 7dp → 10dp in every plot
- [ ] Every figure has a title, axis labels and a legend
- [ ] All tables in `tables/` written and committed

## Report
- [ ] Observed result / interpretation / paper comparison / limitation kept visibly separate
- [ ] Capture-bias caveat attached to Figure 5
- [ ] Pseudobulk limitation stated for the DE results
- [ ] Paper figures cited by their real numbers (atlas = Fig 5, EGFP = Fig 8)
- [ ] Methods state normalisation, HVG count, PCs, resolution, DE test and the integration decision
- [ ] AI Usage Disclosure complete with verbatim prompts

## Git
- [ ] No `.h5ad` and no `data/` in the commit (`git status` before `git add`)
- [ ] `git count-objects -vH` shows a sane repository size
- [ ] Tagged `v1.0-peerreview` and pushed, tag included
