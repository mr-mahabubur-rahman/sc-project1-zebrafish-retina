# §4.7 — Integration tested as a sensitivity analysis

*Insert after §4.6 (Limitations) in the report, renumbering the Conclusions section
if necessary. Also add the summary table to Appendix B, and harmonypy to the
Methods software list.*

---

## 4.7  Batch correction tested, and not adopted

The decision not to integrate (§2.5) rested on an argument: condition and
sequencing batch are completely confounded, so correcting on `sample` risks
removing the injury response along with the batch effect. To test that argument
rather than rely on it, the pipeline was re-run with Harmony correction applied to
the PCA embedding, and the two analyses compared on the quantities this report
depends on. The uncorrected analysis remains the primary result; the corrected
object is retained separately in `results_integrated/`, and the comparison tables
in `tables/integration_check/`.

### The correction removed biological signal preferentially

**Table 17.** Mean variance explained across the first 20 principal components,
before and after Harmony correction on `sample`.

| Factor | Uncorrected | After Harmony | Reduction |
|---|---|---|---|
| condition | 0.0201 | 0.0041 | **79.6%** |
| replicate | 0.0023 | 0.0016 | **30.4%** |
| sample | 0.0244 | 0.0060 | 75.4% |

Replicate is the purely technical term, since replicates within a condition differ
only in handling; condition carries the injury response. A well-behaved batch
correction should therefore remove replicate variance preferentially.

**The opposite occurred.** Harmony removed 79.6% of the condition variance while
removing only 30.4% of the replicate variance, and the ratio between the two fell
from 8.7-fold to 2.6-fold. This is the confound made concrete: because each sample
belongs to exactly one condition, correcting on sample identity necessarily
removes condition differences, and here it did so far more effectively than it
removed the technical variation it was applied to eliminate.

Harmony also stopped before convergence at the default limit of ten iterations,
which is itself consistent with a batch structure that cannot be cleanly separated
from the biological one.

### What the correction damaged

The number of clusters was unchanged at 20, but the partitions differ
substantially in composition (adjusted Rand index 0.62), so correction rearranged
the structure rather than simplifying it. Three integrated clusters fell below 0.8
purity against the original annotation, and the affected populations are precisely
those on which this report's two novel observations rest:

| Integrated cluster | Cells | Majority type | Purity | Mixed with |
|---|---|---|---|---|
| 0 | 3,797 | Rods | 0.735 | Rods (low quality), 24.5% |
| 1 | 928 | Microglia | 0.700 | Unresolved, 27.7% |
| 19 | 39 | Retinal ganglion cells | 0.615 | Müller glia, 15.4% |

**Retinal ganglion cells were largely dispersed.** Of the 97 cells annotated as
RGCs in the primary analysis, only 39 remain in the majority integrated cluster —
40% retained, at 0.615 purity. The divergence from the published analysis reported
in §4.4 would substantially disappear under correction.

**The stressed rod population was dissolved.** Rods and low-quality rods merged
into a single cluster, eliminating the distinction on which §3.4 and §4.5 depend.
The consequence is visible in the composition estimates: the uncorrected rod
series peaks at 3 dpMNU and declines (34.10 → 26.70% at 7 dpMNU), whereas the
integrated series stays flat (38.62 → 38.37%), because the stressed population's
own dynamics have been absorbed into the aggregate.

Two populations survived correction cleanly: horizontal cells (0.926 purity, 651
of 712 cells) and cones (0.970). The pattern is consistent — abundant,
transcriptionally distinct populations were preserved, while small or subtle ones
were degraded.

### What the correction did not change

**Table 18.** Integration sensitivity summary.

| Quantity | Uncorrected | Integrated | Changed? |
|---|---|---|---|
| Number of clusters | 20 | 20 | No |
| Adjusted Rand index between partitions | — | 0.621 | Yes |
| Integrated clusters below 0.8 purity | — | 3 | Yes |
| **EGFP time course, maximum difference** | — | **0.00 pp** | **No** |
| **Rod paralog pairs splitting (of 3)** | **1** | **1** | **No** |

**The reporter time course is exactly invariant.** EGFP-positive percentages are
identical to two decimal places at every timepoint (0.56 / 5.39 / 2.45 / 2.80),
because positivity is a raw transgene count call that does not depend on the
embedding or on clustering. The central result of this report is therefore
unaffected by the batch-correction decision — not approximately, but exactly.

**Müller glia remain the reporter-enriched population**, at 12.12% positive after
correction against 11.34% before, still separated from every other population by a
clear margin. Cell-type enrichment values move by around one percentage point in
either direction, as expected when cells shift between clusters and denominators
change; only the small-population estimates move appreciably, and those were
already flagged as unreliable.

**The rod paralog verdict is unchanged.** After correction, one of three pairs
splits — the same pair, *guca1a*/*guca1b*, and still discordantly with the other
two. *rhol* remains between 1.3 and 1.6 while *rho* ranges 5.8–7.6, so no
*rhol*-high mature population exists under either treatment.

### Consequence for the interpretation in §4.2

Four explanations were offered for the failure to reproduce the immature/mature
rod distinction. **Batch structure can now be eliminated**: the split does not
appear when the batch effect is corrected, and it does not appear when it is not.

The remaining candidates are unchanged in rank:

1. **Normalisation** — now the leading explanation. SCTransform's regularised
   negative binomial residuals treat highly expressed genes differently from
   log-normalisation, and rod phototransduction transcripts are the most abundant
   in the tissue. Testing this would require
   `sc.experimental.pp.normalize_pearson_residuals`, the closest available
   analogue, and is proposed as future work.
2. **Clustering granularity** — untested; would require repeating the rod
   subclustering at the published resolution of 0.2.
3. **Cell number** — untested.
4. **Stress dominating rod variance** — supported by the mitochondrial content and
   reduced transcriptional complexity at matched sequencing depth (§3.10).

### Summary

Integration was tested and not adopted. The test justified the original decision
empirically: correction removed four-fifths of the biological signal while
removing less than a third of the technical variation, and dispersed the two small
populations on which this report's novel observations depend. It also demonstrated
that the two central conclusions — the reporter time course and the rod paralog
verdict — are wholly insensitive to that methodological choice.

---

## Additional edits required elsewhere in the report

**§2.2, software environment** — add: *harmonypy 0.2.0 was installed for the
sensitivity analysis in §4.7 and is not used in the primary pipeline.*

**§2.5, batch structure** — add a closing sentence: *This decision is tested
empirically in §4.7.*

**Appendix A.4, errors identified** — add:

> `scanpy.external.pp.harmony_integrate` transposes Harmony's output, which was
> correct for harmonypy versions below 2.0 (PCs × cells) but fails on 2.0 and
> above (cells × PCs) with an obsm shape error. The integration module calls
> harmonypy directly and orients the result by array shape, so it runs on either
> version. Detected when the wrapper raised a dimension mismatch on a 4,000-cell
> test object.

Also worth adding, if you have not already:

> Two GO terms returned for rods — "embryo development ending in birth or egg
> hatching" and "chordate embryonic development" — are driven by the same
> ribosomal protein genes annotated to developmental phenotypes. This is a known
> artefact of ribosomal gene annotation in zebrafish and is not evidence of a
> developmental programme.

**Appendix B, reproducibility** — add: *The integration sensitivity analysis
(step 08) writes to `results_integrated/` and `tables/integration_check/`. The
former is git-ignored; the comparison tables are committed.*

**Figure captions to update after re-running step 07:**

- **Figure 18 (rod volcano):** *Axis limits are set from the significant genes;
  6,133 non-significant genes with extreme fold changes at near-zero expression
  fall outside the plotted range. Every labelled upregulated gene except* rbp4l
  *and* guca1a *is mitochondrial, while the downregulated genes are
  phototransduction components.*
- **Figure 19 (rod GO):** *Hatched bars mark six terms whose p-values fell below
  floating-point resolution and are plotted at the smallest resolvable value
  (1.5×10⁻¹⁵); their true significance is higher than shown. Eight terms plotted
  from 179 input genes.*
- **Figure 23 (cone volcano):** *Axis limits set from the significant genes. The
  upregulated population is dense at modest fold changes while the downregulated
  side is sparser with larger effects.*
