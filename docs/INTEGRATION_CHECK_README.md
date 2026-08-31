# Step 08 — Integration sensitivity analysis

## What this is

The main pipeline does not apply batch correction, for the reason in §2.5 of the
report: condition and sequencing batch are completely confounded by design, so
correcting on `sample` risks removing the injury response along with the batch
effect.

That was an argument. This runs the analysis **with** Harmony and asks which
conclusions would change if the decision had gone the other way. It is a
sensitivity analysis, not a replacement — the uncorrected result stays primary.

## Install

Three files, then one package:

| File | Destination |
|---|---|
| `integration_check.py` | `scripts\` |
| `step_08_Integration_Sensitivity.ipynb` | `notebooks\` |
| `build_notebooks.py` | `tools\` (overwrites the outdated one) |

```powershell
cd D:\retina\project1_zebrafish_retina
uv pip install harmonypy
uv pip freeze | Out-File -FilePath requirements.txt -Encoding ascii
```

Note the `Out-File ... -Encoding ascii` — PowerShell's `>` writes UTF-16, which
git treats as binary.

## Run

Open `step_08_Integration_Sensitivity.ipynb`, select the **Python 3 (retina)**
kernel, Run All. It loads `results/03_preprocessed.h5ad` and
`results/05_annotated.h5ad`, so steps 01–05 must have been run.

Expect a few minutes; Harmony iterates to convergence.

## Outputs

Written to `tables/integration_check/` (small, committed) and
`results_integrated/08_integrated.h5ad` (large, git-ignored):

| File | Contents |
|---|---|
| `batch_variance_before_after.csv` | PC variance by condition/replicate/sample, both embeddings |
| `cluster_agreement.csv` | Cluster counts and adjusted Rand index |
| `integrated_cluster_purity.csv` | How cleanly each integrated cluster maps to one cell type |
| `egfp_timecourse_comparison.csv` | The reporter series, both analyses |
| `egfp_enrichment_comparison.csv` | EGFP% per cell type, both analyses |
| `rod_paralogs_integrated.csv` | The paralog test on corrected data |
| `proportions_comparison.csv` | Cell-type proportions, both analyses |
| `integration_sensitivity_summary.csv` | One table: what changed, what did not |

## How to read the result

**The EGFP time course must be identical.** Positivity is a raw transgene count
call that does not depend on clustering. If it moves, something is wrong — treat
that as a bug, not a finding.

**Watch the variance table carefully.** Harmony will reduce `replicate` variance,
which is the point. Watch what it does to `condition`: if condition falls by a
similar proportion, the correction is removing biology along with batch, and that
is the concrete form of the risk §2.5 describes.

**Cluster purity below 0.8** means integration merged populations that were
previously distinct. Which ones matters — merging two cone clusters is
unremarkable; absorbing the 97 retinal ganglion cells would be worth reporting.

**The rod paralog verdict is the substantive test.** The uncorrected analysis had
1 of 3 pairs splitting, discordantly. If integration gives a coordinated split
across all three, the failure was in batch structure rather than in normalisation
— a genuine finding. If it stays at 1, normalisation and the stress axis remain
the leading explanations, and §4.2 of the report is strengthened.

## What to write up

Add a short subsection to the Discussion recording:

1. That integration was **tested** as a sensitivity analysis, not adopted.
2. How much batch variance Harmony removed, and how much condition variance with it.
3. Whether the EGFP time course moved (it should not).
4. Whether Müller glia remained the reporter-enriched population.
5. Whether the rod paralog verdict changed.

Then update `docs/method_comparison.md`, add harmonypy to `requirements.txt`, and
record the harmonypy version in the Methods.

## A note on the implementation

`scanpy.external.pp.harmony_integrate` is **not** used. That wrapper transposes
Harmony's output, which was correct for harmonypy < 2.0 (PCs × cells) but fails on
2.0+ (cells × PCs) with an obsm shape error. This module calls harmonypy directly
and orients the result by shape, so it works on either version. Worth noting in
Appendix A.4 as a caught incompatibility.
