# §4.2 revised — with both sensitivity analyses

*Replaces the existing §4.2 in the report. Two of the four candidate explanations
are now eliminated by direct experiment.*

---

## 4.2  Why the reporter reproduced and the rods did not

The two pipelines differ at every major processing step, and a substantial
difference in analysed cell number separates them. That the reporter kinetics
nonetheless agree to within 0.2 percentage points at every timepoint indicates the
measurement is robust to analytical choice. This is expected on mechanistic
grounds: EGFP positivity is a near-binary determination on a transgene with no
endogenous paralog, no ambiguity of gene identity, and a strong biological signal.
Such a measurement is insensitive to normalisation and clustering decisions, since
it depends only on whether a transcript was captured. The sensitivity analysis in
§4.7 confirms this directly — the time course is identical to two decimal places
under batch correction.

The rod distinction is the opposite case. It rests on graded differences in the
relative expression of paralogous genes, quantities directly affected by the
normalisation method.

### The cone result excludes a pipeline limitation

The same code, the same subclustering resolution and the same paralog-inversion
logic reproduced the *arr3a*/*arr3b* split exactly (§3.11). The pipeline is
therefore demonstrably capable of resolving photoreceptor subtype structure, and
the rod failure is specific to that comparison rather than a general limitation.

The distinction between the two cases is instructive. UV cone identity is
**categorical** — a cone expresses *arr3b* or *arr3a*, and the ratio differs by an
order of magnitude between populations. Rod maturation is **graded** — both
paralogs are expressed in every rod, and the classification depends on relative
levels. Categorical identity survives a change of normalisation; graded ratios in
highly expressed genes need not.

### Four explanations, two now eliminated

**1. Normalisation — untested, and now the only remaining candidate.**
SCTransform models expression as regularised negative binomial residuals
(Hafemeister and Satija, 2019), a different estimator from log-normalisation rather
than a variant of it. It weights highly expressed genes differently, and rod
phototransduction transcripts are the most abundant in retinal tissue — precisely
the regime where the two methods diverge most. Testing this would require
`sc.experimental.pp.normalize_pearson_residuals`, the closest available analogue to
SCTransform in Scanpy, and is proposed as future work.

**2. Clustering granularity — eliminated (§4.8).**
The published analysis clustered at resolution 0.2; the primary analysis here used
0.6 and tested the paralogs by subclustering within an annotated rod population.
Re-clustering the whole dataset at 0.2 and repeating the test on the rod clusters
that emerged gives the same verdict: one of three pairs splits, discordantly, and
no cluster carries the coordinated pattern the published result requires.

**3. Cell number — untested.**
A population defined by graded expression differences may not separate at the cell
numbers available here. This cannot be tested without additional data.

**4. Stress dominating rod variance — supported by the evidence.**
The rod subclusters differ in mitochondrial content and transcriptional complexity
*at matched sequencing depth* (§3.10), and the stressed population is 98.9%
injury-derived. If stress is the largest axis of variation among rods, it may
partition the cells before any maturation axis can.

**5. Batch structure — eliminated (§4.7).**
Applying Harmony correction does not recover the split; the verdict is unchanged at
one of three pairs, the same pair, and still discordant.

### Three independent tests, one answer

The paralog test was applied under three different analytical treatments:

| Treatment | Rod groups | Pairs splitting | Coordinated pattern |
|---|---|---|---|
| Primary (resolution 0.6, subclustered) | 4 subclusters | 1 of 3 | No |
| Harmony-corrected (§4.7) | 3–4 subclusters | 1 of 3 | No |
| Resolution 0.2, global clustering (§4.8) | 4 clusters | 1 of 3 | No |

In every case the pair that splits is *guca1a*/*guca1b*, and in every case it splits
discordantly with the other two. **In no treatment does any group carry the
*rhol*-high mature signature**: *rhol* never exceeds 1.7 while *rho* ranges 4.3–7.6,
under every condition tested.

The failure to reproduce is therefore not attributable to the two analytical
choices that most plausibly differed from the published work. What remains is
normalisation, and the possibility that stress structure in these data dominates
whatever maturation signal is present.

---

# §4.8 — new subsection

*Insert after §4.7 (the integration sensitivity analysis).*

## 4.8  Clustering resolution tested as a sensitivity analysis

Bise *et al.* clustered at resolution 0.2 and reported 17 clusters, among them two
rod populations. The primary analysis here clustered at 0.6 and tested the paralogs
by subclustering *within* an annotated rod population — a different procedure at a
different granularity, and one of the candidate explanations in §4.2. It was
therefore tested directly.

### Rod clusters must be identified relative to background

An absolute expression threshold does not identify rod clusters in this dataset. At
resolution 0.2, mean *rho* ranges from 3.09 to 7.49 across the twelve clusters with
no gap: every cluster carries substantial ambient *rho*, because rods dominate the
tissue and lyse readily during dissociation. This is the same effect that made
gene-set panel scores unusable in §2.6, and an initial attempt using an absolute
cut of *rho* ≥ 3.0 selected all twelve clusters.

Selection was therefore made relative — *rho* above the 60th percentile of the
cluster means — with an additional requirement on *nr2e3*, a rod-specific
transcription factor expressed at much lower absolute levels and correspondingly
less exposed to ambient contamination. *nr2e3* discriminates cleanly where *rho*
does not: 0.148–0.322 in the four selected clusters against 0.043–0.096 in the
eight rejected ones, a three- to seven-fold separation.

**Table 19.** Rod clusters identified at resolution 0.2.

| Cluster | Cells | *rho* | *nr2e3* | *gnat1* | Selected |
|---|---|---|---|---|---|
| 3 | 3,470 | 7.487 | 0.148 | 5.075 | **Yes** |
| 0 | 1,894 | 6.275 | 0.322 | 4.538 | **Yes** |
| 1 | 2,922 | 5.098 | 0.291 | 3.402 | **Yes** |
| 5 | 1,547 | 4.293 | 0.320 | 2.544 | **Yes** |
| 8 | 199 | 4.500 | 0.077 | 2.222 | No |
| 9 | 851 | 3.866 | 0.053 | 2.039 | No |
| 11 | 51 | 3.612 | 0.046 | 1.469 | No |
| 4 | 1,709 | 3.357 | 0.081 | 1.791 | No |
| 7 | 580 | 3.210 | 0.043 | 1.514 | No |
| 6 | 129 | 3.180 | 0.057 | 1.510 | No |
| 2 | 1,540 | 3.154 | 0.096 | 1.538 | No |
| 10 | 1,001 | 3.087 | 0.043 | 1.338 | No |

### The paralog test at the published granularity

**Table 20.** Paralog expression across the four rod clusters at resolution 0.2.

| Cluster | Cells | *rho* | *rhol* | *pde6ga* | *pde6gb* | *guca1a* | *guca1b* | Pattern |
|---|---|---|---|---|---|---|---|---|
| 0 | 1,894 | 6.275 | 1.689 | 4.359 | 4.792 | 1.671 | 1.899 | ABB |
| 1 | 2,922 | 5.098 | 1.006 | 3.000 | 3.475 | 0.919 | 1.143 | ABB |
| 3 | 3,470 | 7.487 | 1.256 | 3.894 | 5.594 | 2.777 | 1.029 | ABA |
| 5 | 1,547 | 4.293 | 0.507 | 2.228 | 2.891 | 0.755 | 0.634 | ABA |

*Pattern reads A where the first paralog exceeds the second and B where the second
exceeds the first, for the three pairs in column order. The published structure
requires one cluster reading AAA and another BBB.*

*rho* exceeds *rhol* in all four clusters. *pde6gb* exceeds *pde6ga* in all four.
Only *guca1a*/*guca1b* splits, two clusters each way, and discordantly with the
other pairs. **No cluster reads AAA and none reads BBB**, so the coordinated
reciprocal pattern the published result requires does not appear at the published
resolution.

### Resolution is not portable between algorithms

Clustering at resolution 0.2 produced **12 clusters against the 17 reported**.
Leiden and Seurat's SNN modularity optimisation solve related but distinct
problems, so the same nominal resolution parameter does not yield the same
partition. This is worth stating in its own right: matching a published resolution
value does not match a published granularity, and any reanalysis that assumes
otherwise is comparing different things.

### Conclusion

Clustering granularity does not explain the failure to reproduce the rod split. Of
the four explanations offered in §4.2, this and batch structure (§4.7) are now
eliminated by direct experiment; normalisation remains the leading candidate.

---

## Edits required elsewhere

**§2.2, software environment** — no change (step 09 adds no dependency).

**§3.10** — add after the paralog verdict:

> This verdict is unchanged when the analysis is repeated at the published
> clustering resolution of 0.2 (§4.8) and when batch correction is applied (§4.7).
> In no treatment does any rod group carry the *rhol*-high signature.

**Table 15 (panel status)** — the Figure 6B row can now read:

> Not reproduced — 2 of 3 pairs fail at resolution 0.6, 0.2 and under batch
> correction (§3.10, §4.7, §4.8)

**Appendix A.4** — add:

> An initial attempt to identify rod clusters at resolution 0.2 used an absolute
> threshold of mean *rho* ≥ 3.0 and selected all twelve clusters, because ambient
> *rho* contaminates every droplet — the same effect documented in §2.6 as making
> gene-set panel scores unusable. Detected because the summary reported twelve rod
> clusters out of twelve. Corrected by selecting relative to the distribution of
> cluster means and additionally requiring *nr2e3*, a rod-specific transcription
> factor much less exposed to ambient contamination. The function now warns
> explicitly when a criterion selects every cluster.

**Appendix B** — add:

> The resolution sensitivity analysis (step 09) writes to
> `tables/resolution_check/` and `figures/resolution_check/`.
