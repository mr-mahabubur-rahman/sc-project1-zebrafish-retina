# Independent Scanpy reanalysis of *careg:EGFP* reporter dynamics and Müller glia activation in the regenerating zebrafish retina after MNU-induced photoreceptor ablation

**Mahabubur Rahman**

Repository: https://github.com/mr-mahabubur-rahman/sc-project1-zebrafish-retina · Tag: `v1.0-peerreview`
Primary reference: Bise *et al.* (2023), *Frontiers in Molecular Neuroscience* 16:1160707
Data: NCBI GEO **GSE202212**

> **[VERIFY]** marks citations to confirm against the original source.
> Appendix A must be completed by the author.

---

## Abstract

Zebrafish regenerate lost photoreceptors by reprogramming Müller glia into
proliferative progenitors, a capacity absent in mammals. Bise *et al.* (2023)
identified the regeneration-responsive element *careg* as a reporter of activated
Müller glia after methylnitrosourea (MNU) photoreceptor ablation. Their
single-cell RNA-sequencing data (control and 3, 7 and 10 days post-MNU, two
replicates each) were reanalysed here using an independent Python/Scanpy pipeline
rather than the published Seurat workflow, to establish which conclusions are
robust to analytical choice. Of 20,097 cells, 15,893 (79.1%) passed quality
control and resolved into thirteen retinal populations. The reporter result
reproduced closely: EGFP-positive cells rose from 0.56% of control cells to 5.40%
at 3 dpMNU, declining to 2.45% and 2.80% at 7 and 10 dpMNU — within 0.2 percentage
points of the published values throughout. Müller glia were enriched 3.75-fold and
held 42.6% of all reporter-positive cells. Because EGFP-positive cells proved to
be sequenced eleven-fold more deeply than negative cells, positivity was
re-examined within depth quartiles; at matched depth the induction was stronger,
not weaker, rising from 0% of control cells to 54.3% at 3 dpMNU. All 23 genes
tested from the published EGFP-positive signature were recovered as significantly
upregulated. Müller glia subclustering resolved a reporter-active, *ascl1a*- and
*mki67*-positive population with reduced glial identity markers. The published
immature/mature rod distinction did not reproduce: the diagnostic paralog
inversion failed for two of three gene pairs, while the equivalent test in cones
reproduced exactly, identifying the rod result as specific rather than a pipeline
limitation.

*(248 words — trim to ~200 if enforced.)*

**Keywords:** zebrafish, retina regeneration, Müller glia, single-cell RNA-seq,
*careg*, reproducibility

---

## 1. Introduction

### 1.1 Retinal structure and the regenerative divide

The vertebrate retina is a laminated neuroepithelium of three nuclear and two
plexiform layers, conserved between zebrafish and humans (Raymond and Barthel,
2004). Photoreceptors occupy the outer nuclear layer and synapse onto bipolar and
horizontal cells; Müller glia span the full retinal thickness, providing metabolic
and structural support.

Despite this conserved organisation, the consequences of photoreceptor loss
diverge. In mammals photoreceptor death is irreversible; zebrafish restore both
retinal structure and visual function (Wan and Goldman, 2016). The difference does
not arise from a dedicated stem cell population but from the intrinsic plasticity
of Müller glia.

### 1.2 Müller glia as facultative stem cells

After injury, zebrafish Müller glia dedifferentiate, re-enter the cell cycle and
divide asymmetrically to produce multipotent progenitors that differentiate into
replacement neurons (Bernardos *et al.*, 2007; Nagashima *et al.*, 2013; Goldman,
2014). Müller glia were first identified as a potential source of neural
regeneration in the postnatal chicken retina (Fischer and Reh, 2001). The
reprogramming step — quiescent glial support cell to neurogenic progenitor —
remains incompletely understood at the molecular level.

### 1.3 The *careg* reporter

Non-coding regulatory elements provide biosensors for the switch from quiescence
to regenerative mobilisation. Pfefferli and Jaźwińska (2017) identified *careg*, a
3.18 kb element from the *ctgfa* promoter, transiently upregulated in
regeneration-participating cells of fin and heart. Bise *et al.* (2023) extended
this to neural tissue: *careg:EGFP* is silent in uninjured retina, induced in a
subset of Müller glia from one day after MNU treatment, sustained throughout
regeneration, and extinguished by 90 dpMNU.

MNU is an alkylating agent inducing apoptotic photoreceptor death, predominantly
of rods, established as a non-surgical retinal degeneration model in zebrafish
(Tappeiner *et al.*, 2013; Maurer *et al.*, 2014). Because the transgene was added
to the alignment reference, reporter transcripts are countable per cell alongside
the endogenous transcriptome.

### 1.4 Rationale for independent reanalysis

Reproducibility has two distinct meanings in computational biology. Re-executing
published code establishes only that a pipeline is deterministic. Recovering the
same conclusion through independent methodological choices tests whether a finding
survives the decisions every analyst must make. The second is the stronger claim,
and the one attempted here.

This reanalysis uses Python/Scanpy (Wolf *et al.*, 2018) rather than Seurat,
log-normalisation rather than SCTransform, and Leiden community detection (Traag
*et al.*, 2019) at a different resolution, with entirely independent annotation.

### 1.5 Research questions

1. Which retinal populations are recoverable?
2. How does composition shift across control → 3 → 7 → 10 dpMNU?
3. Which populations express the reporter, when does it peak, does it subside?
4. Can Müller glia activation states be distinguished?
5. How do rod and cone transcriptomes respond to injury?

---

## 2. Materials and Methods

### 2.1 Dataset

Eight 10x Genomics Chromium 3′ filtered feature-barcode matrices from GEO
GSE202212: two uninjured controls (heat-inactivated MNU, controlling for chemical
exposure independent of injury) and two replicates each at 3, 7 and 10 dpMNU. The
reference was GRCz11 with the EGFP transgene added, giving 25,433 features; 22,813
genes were retained after filtering.

**Table 1. Cells per sample, before and after quality control.**

| Sample | Condition | Before | After | % lost | Median mito (%) |
|---|---|---|---|---|---|
| ctrl1 | Control | 1,155 | 906 | 21.6 | 5.96 |
| ctrl2 | Control | 2,116 | 1,609 | 24.0 | 5.48 |
| 3dp1 | 3 dpMNU | 2,873 | 2,317 | 19.4 | 8.74 |
| 3dp2 | 3 dpMNU | 2,388 | 1,844 | 22.8 | 6.20 |
| 7dp1 | 7 dpMNU | 3,444 | 2,562 | 25.6 | **10.30** |
| 7dp2 | 7 dpMNU | 1,480 | 1,165 | 21.3 | 6.31 |
| 10dp1 | 10 dpMNU | 4,124 | 3,333 | 19.2 | 8.09 |
| 10dp2 | 10 dpMNU | 2,517 | 2,157 | 14.3 | 6.65 |
| **Total** | | **20,097** | **15,893** | **20.9** | |

Recovery was unbalanced between conditions (controls 3,271 cells; 10 dpMNU 6,641)
and between replicates within a condition (7dp1 3,444 against 7dp2 1,480). A
systematic difference also exists between the two runs of each injured pair:
**replicate 1 consistently shows higher mitochondrial content than replicate 2**
(8.74 vs 6.20 at 3dp; 10.30 vs 6.31 at 7dp; 8.09 vs 6.65 at 10dp), indicating a
processing or dissociation difference between runs. This becomes important in §3.4.

All proportion analyses were computed **within each sample**, with replicate values
reported individually.

### 2.2 Software environment

Python 3.11.15, Scanpy 1.11.5 (Wolf *et al.*, 2018), AnnData 0.11.4, NumPy 2.4.6,
pandas 2.3.3, SciPy 1.17.1, matplotlib 3.11.1, leidenalg 0.12.0, gseapy 1.3.1.
Environments managed with `uv`; versions in `docs/session_info.json`, pinned in
`requirements.txt` and `uv.lock`. Random seed 0 throughout.

### 2.3 Quality control

Zebrafish gene symbols are lower case: mitochondrial genes were identified by the
`mt-` prefix, ribosomal by `rps`/`rpl`, haemoglobins by `hba`/`hbb`/`hb-`. Applying
the human `MT-` convention would match nothing and silently disable the
mitochondrial filter; this was verified explicitly.

Cells were retained with ≥200 detected genes and <15% mitochondrial counts; genes
detected in ≥3 cells. Threshold diagnostics were inspected before filtering.

### 2.4 Normalisation, feature selection and clustering

Counts were normalised to 10,000 per cell and log1p-transformed. Raw counts were
preserved in a layer and the complete log-normalised matrix in `.raw`.

**This distinction is not cosmetic.** Feature selection reduced the matrix to 2,000
highly variable genes for clustering, which is standard and appropriate — but all
differential expression and marker discovery were performed on the full 22,813-gene
matrix retained in `.raw`. Restricting statistical tests to the highly variable
subset would exclude 91% of the transcriptome, including genes central to the
comparisons attempted here (*meig1*, *guca1b*, *cx52.6*).

2,000 HVGs were selected with `sample` as batch key; technical covariates
(`total_counts`, `pct_counts_mt`) regressed out; matrix scaled with clipping at ±10.
PCA was followed by a 15-nearest-neighbour graph on 20 principal components, UMAP,
and Leiden clustering at resolution 0.6, giving 20 clusters.

### 2.5 Batch structure and the decision not to integrate

Bise *et al.* applied Seurat anchor-based integration. **Integration was not applied
here**, on evidence rather than by default. Variance in the first 20 principal
components was decomposed by factor:

| Factor | Mean variance explained |
|---|---|
| condition | 0.0201 |
| replicate | 0.0023 |
| sample | 0.0244 |

Replicate is the informative term, since replicates within a condition differ only
technically. **Condition explains nearly nine-fold more PC variance than replicate**,
so technical batch structure is small relative to biological signal. Correcting on
sample identity would risk removing the injury response, which is confounded with
sequencing run by design. That confound cannot be resolved computationally and is
stated as a limitation.

### 2.6 Cell type annotation

Cluster markers were identified by Wilcoxon rank-sum test on the full gene set.
Annotation was manual, requiring concordant evidence from multiple markers, using
the guide's canonical panel supplemented by markers from Ogawa and Corbo (2021)
and Hoang *et al.* (2020). Each assignment is recorded with supporting markers,
conflicting evidence, confidence and rationale in `tables/cluster_annotations.csv`.
Clusters without confident support were labelled **Unresolved**.

A methodological observation arose during annotation. `sc.tl.score_genes` panel
scores proved unreliable, assigning the highest rod score to clusters whose markers
were unambiguously oligodendrocyte (*mbpa*, *plp1b*, *cldnk*) and mesenchymal
(*col1a1b*, *col1a2*, *pdgfra*). This is attributable to ambient rod mRNA: rods
dominate the tissue and lyse readily during dissociation, so free *rho* transcripts
contaminate all droplets and inflate rod scores globally. Panel scores were
rejected and annotation performed from data-driven marker rankings.

### 2.7 Reporter positivity

A cell was EGFP-positive if its **raw transgene count was ≥1**. The observed count
distribution justifies this: it decays monotonically from 310 cells at one count
(76 at two, 34 at three, trailing to 14) with no bimodal separation, so any higher
threshold would be arbitrary. **64% of positive cells carry exactly one count**, and
requiring ≥2 would reduce the positive fraction from 3.03% to 1.08%. Continuous
log-normalised values are reported alongside every binary call.

### 2.8 Differences from the published analysis

| Step | Bise *et al.* (2023) | This analysis | Equivalent? |
|---|---|---|---|
| Cell QC | scater; ~12 ± 2% removed as doublets/dying | ≥200 genes, <15% mito | No |
| Doublet removal | explicit | **none** | No |
| Normalisation | SCTransform (Hafemeister and Satija, 2019) | `normalize_total` + `log1p` | No |
| Integration | Seurat anchors, 30 dims (Stuart *et al.*, 2019) | **none** (§2.5) | No |
| Clustering | SNN modularity, resolution 0.2 | Leiden, resolution 0.6 | No |
| Marker/DE test | Wilcoxon rank-sum | Wilcoxon rank-sum | **Yes** |
| Enrichment | topGO, Fisher exact | gseapy/Enrichr (`organism='zebrafish'`) | No |

SCTransform models expression as regularised negative binomial residuals — a
different estimator, not a variant. It weights highly expressed genes differently,
precisely the regime occupied by rod phototransduction transcripts. This is
relevant to §3.8.

**UMAP embeddings are not reproducible as images.** Coordinates depend on
normalisation, feature selection, PCA and random initialisation, and cluster
numbering is arbitrary between analyses. All comparisons are made on quantities.

---

## 3. Results

### 3.1 Quality control (Figure 1)

15,893 of 20,097 cells were retained (79.1%). Loss ranged 14.3–25.6% per sample.
Mitochondrial content was confirmed non-zero, verifying the `mt-` prefix matched
real genes.

Cell loss tracks the replicate rather than the condition (Table 1): the first run
of each injured pair consistently carries higher mitochondrial content. Injured
retinas dissociate into a population containing more stressed singlets, but the
magnitude of that effect differs between runs — a technical source of variation
that propagates into §3.4.

### 3.2 Thirteen retinal populations, including horizontal cells and RGCs (Figures 2, 3)

**Table 2. Annotated cell populations.**

| Cell type | Clusters | Cells | % | Key markers |
|---|---|---|---|---|
| Rods | 0, 2 | 4,314 | 27.1 | *rho*, *pde6a/b*, *gnb1b*, *guca1a* |
| Cones | 9, 10 | 2,614 | 16.4 | *gnat2*, *pde6c*, *arr3a*, *gnb3b* |
| Müller glia | 6, 7 | 1,807 | 11.4 | *rlbp1a*, *glula*, *slc1a2b*, *pleca* |
| Bipolar cells | 13, 17 | 1,758 | 11.1 | *cabp5a*, *grin1b*, *trpm1b* |
| Amacrine cells | 3 | 1,274 | 8.0 | *slc32a1*, *pax6a*, *snap25a* |
| Rods (low quality) | 4 | 1,062 | 6.7 | *rho* with *mt-co2* dominant |
| RPE | 11, 12 | 847 | 5.3 | *rpe65a*, *tyrp1b*, *pmela*, *stra6* |
| **Horizontal cells** | 8 | 712 | 4.5 | ***cx52.6***, *rem1*, *tnr*, *pcdh8* |
| Microglia | 16 | 661 | 4.2 | *apoc1*, *ccl34b.1*, *apoeb* |
| Unresolved | 1, 14 | 365 | 2.3 | — |
| Erythrocytes | 5 | 202 | 1.3 | *hbba1*, *hbaa1* |
| Oligodendrocytes | 19 | 129 | 0.8 | *mbpa*, *plp1b*, *cldnk* |
| **Retinal ganglion cells** | 15 | 97 | 0.6 | ***rbpms2b***, *gap43*, *islr2* |
| Pericytes | 18 | 51 | 0.3 | *col1a1b*, *col1a2*, *pdgfra* |

Confidence audit: **82.9% of cells lie in high-confidence assignments**, 14.4%
medium, 2.6% low. The three low-confidence clusters (1, 14, 18) support no
conclusion in this report.

**Horizontal cells** were identified by *cx52.6* (1.010) and *rem1* (1.342), each
five- to ten-fold above any other cluster. *cx52.6* encodes a connexin specific to
zebrafish horizontal cell gap junctions **[VERIFY]**. *cabp5a* was present (0.770)
but far below the bipolar clusters (3.004, 2.736), excluding bipolar identity.

**Retinal ganglion cells** were recovered as a discrete *rbpms2b*-positive cluster.
Bise *et al.* state their analysis captured very few ganglion cells and that no
cluster showed unique *rbpms2b*. This is a genuine divergence (§4.4).

Two clusters remain unresolved: cluster 1 (15 cells; 12,462 genes and 398,264
counts per cell, ~100-fold the dataset median) is multiplet droplets — the
population removed by the published doublet filtering; cluster 14 (350 cells)
expresses only cytoskeletal and heat-shock genes.

### 3.3 Composition shifts are dominated by capture bias (Figure 5)

**Table 3. Composition across the time course (%, mean of replicate values).**

| Cell type | Control | 3 dpMNU | 7 dpMNU | 10 dpMNU |
|---|---|---|---|---|
| Rods | 12.16 | **33.24** | 25.46 | 28.36 |
| Rods (low quality) | 0.66 | 5.00 | **9.82** | 6.09 |
| Cones | 26.05 | 16.80 | 14.30 | 15.02 |
| Bipolar cells | 25.76 | 10.64 | 9.01 | 6.77 |
| Müller glia | 8.39 | 11.67 | 10.89 | **13.52** |
| Horizontal cells | 9.29 | 3.84 | 4.47 | 2.85 |
| Microglia | 2.22 | 4.05 | **4.96** | 4.74 |
| RPE | 5.08 | 2.86 | 5.54 | 7.88 |
| Amacrine cells | 8.48 | 4.56 | 11.22 | 8.95 |

Rods rose 2.7-fold, from 12.16% to 33.24%. Bise *et al.* report the same direction
and attribute it to bias in cell survival and capture efficiency, noting comparable
biases in other vertebrate retinal injury datasets (Macosko *et al.*, 2015).

**Declines are not interpretable.** Cones fall to 15.02%, bipolar cells to 6.77%,
horizontal cells to 2.85%. None plausibly dies at that rate, and the published
immunofluorescence reports no such loss. Proportions are constrained to sum to
unity, so over-capture of rods mechanically depresses everything else.

**Increases against that dilution are meaningful.** Microglia rise from 2.22% to
4.96% *despite* the rod inflation, so the underlying increase must be larger than
the proportion suggests — consistent with an inflammatory response to
photoreceptor death. Müller glia (8.39 → 13.52%) and RPE (5.08 → 7.88%) rise
likewise. The Müller glia increase is corroborated independently by the reporter
data (§3.5–3.6).

### 3.4 A stressed rod population, with important replicate variability

Cluster 4 (1,062 cells) carried rod identity but was topped by *mt-co2*, with 12.6%
median mitochondrial content and 320 genes per cell against 2.6% and 496 for
healthy rods. It was annotated **Rods (low quality)**.

Its abundance rises with injury (0.66% of control cells to 9.82% at 7 dpMNU), but
**the replicate values disagree substantially**:

| | rep 1 | rep 2 | ratio |
|---|---|---|---|
| ctrl | 0.88 | 0.44 | 2.0× |
| 3dp | 8.76 | 1.25 | 7.0× |
| 7dp | **17.25** | **2.40** | **7.2×** |
| 10dp | 7.50 | 4.68 | 1.6× |

Replicates within a condition differ only technically, so a seven-fold gap requires
a technical explanation — and the QC data provide one. The mitochondrial content of
each sample (Table 1) tracks these values directly: 7dp1 has the highest
mitochondrial content in the dataset (10.30%) and the highest low-quality rod
fraction (17.25%), while 7dp2 has 6.31% and 2.40%.

**The conclusion this supports, and its limit.** A stressed rod population is
consistently present after injury and near-absent from controls — every injured
sample exceeds every control. Its abundance, however, tracks per-sample
dissociation quality, so the magnitude is not reliably estimated and **the apparent
peak at 7 dpMNU cannot be established from these data**. The population is real and
injury-associated; its time course is not determined.

### 3.5 The *careg:EGFP* time course reproduces closely (Figure 4)

**Table 4. EGFP-positive cells across the time course.**

| Condition | This analysis | Replicate values | Bise *et al.* | Difference |
|---|---|---|---|---|
| Control | **0.56%** | 0.55, 0.56 | 0.54% | +0.02 |
| 3 dpMNU | **5.40%** | 3.63, 7.16 | 5.40% | 0.00 |
| 7 dpMNU | **2.45%** | 3.43, 1.46 | 2.64% | −0.19 |
| 10 dpMNU | **2.80%** | 2.16, 3.43 | 2.66% | +0.14 |

Every timepoint agrees within 0.2 percentage points, with identical kinetics: a
near-silent baseline, a sharp peak at 3 dpMNU, and a partial decline that **does not
return to control levels by day 10** — approximately five-fold above control. This
answers research question 3 directly: expression declines from peak but does not
subside by day 10, consistent with the published immunofluorescence, in which
reporter protein persists to 30 dpMNU and is undetectable only at 90 dpMNU.

Three qualifications. The two controls agree to within 0.01 percentage points, so
the baseline is well determined. Injured conditions vary about two-fold between
replicates, and the direction is not systematic (replicate 2 higher at 3dp,
replicate 1 higher at 7dp), so this is per-sample variation rather than a batch
offset. The overlapping ranges at 7dp (1.46–3.43) and 10dp (2.16–3.43) mean **those
two timepoints are not distinguishable from one another**.

EGFP transcripts were detected at low frequency in every population, including
erythrocytes (0.50%) and RPE (0.24%), neither of which can transcribe a *careg*
reporter. This establishes an **ambient contamination floor of ~0.2–0.5%**, within
which the control value of 0.56% falls (see §3.7 for the resolution of this).

### 3.6 Müller glia are the reporter-enriched population

**Table 5. EGFP positivity by population.**

| Cell type | Cells | EGFP⁺ | % | Enrichment | Share of all EGFP⁺ |
|---|---|---|---|---|---|
| **Müller glia** | 1,807 | 205 | **11.34** | **3.75×** | **42.62%** |
| Pericytes | 51 | 3 | 5.88 | 1.94 | 0.62 |
| Retinal ganglion cells | 97 | 5 | 5.15 | 1.70 | 1.04 |
| Amacrine cells | 1,274 | 48 | 3.77 | 1.24 | 9.98 |
| Cones | 2,614 | 83 | 3.18 | 1.05 | 17.26 |
| Horizontal cells | 712 | 19 | 2.67 | 0.88 | 3.95 |
| Bipolar cells | 1,758 | 33 | 1.88 | 0.62 | 6.86 |
| Rods | 4,314 | 47 | 1.09 | 0.36 | 9.77 |
| Erythrocytes | 202 | 1 | 0.50 | 0.16 | 0.21 |
| RPE | 847 | 2 | 0.24 | 0.08 | 0.42 |

Müller glia comprise 11.4% of cells but hold 42.6% of reporter-positive cells. Bise
*et al.* report EGFP in nearly 10% of Müller glia, describing it as an outstanding
proportion among cell types; the present 11.34% agrees.

Nine of thirteen populations fall between 0.24% and 3.77%, forming a continuum
bracketed by erythrocytes and RPE — the ambient floor. Müller glia sit three-fold
above the next population. Pericytes (3 of 51 cells) and RGCs (5 of 97) show
nominally high percentages on denominators too small to interpret.

### 3.7 The reporter result survives correction for sequencing depth

EGFP-positive Müller glia proved to be sequenced far more deeply than negative
cells — **15,044 median counts against 1,360, an eleven-fold difference** — with
identical mitochondrial content (5.4% vs 5.3%). Deeper cells are more likely to
capture at least one transgene transcript by sampling alone, so positivity might
reflect detection rather than expression. This is the most serious potential
confound in the analysis and was tested directly.

**Table 6. EGFP-positive rate (%) within Müller glia depth quartiles.**

| Depth quartile | Median counts | n | ctrl | 3dp | 7dp | 10dp |
|---|---|---|---|---|---|---|
| Q1 | 627 | 455 | 3.57 | 5.88 | 2.22 | 2.89 |
| Q2 | 1,042 | 449 | 1.27 | 3.00 | 2.50 | 5.26 |
| Q3 | 2,640 | 451 | 5.26 | **13.40** | 5.17 | 4.97 |
| Q4 | 14,726 | 452 | **0.00** | **54.35** | 28.44 | 18.82 |

Both effects are real and separable. Depth strongly affects detection: positivity
rises from 3.74% in Q1 to 31.19% in Q4. **But the injury effect is independent of
depth.** In Q4, where every cell has comparable sequencing depth, **no control cell
is EGFP-positive while 54.35% of 3 dpMNU cells are**, and the full kinetics are
preserved (0 → 54.35 → 28.44 → 18.82). Q3 shows the same pattern.

Depth cannot explain a difference between groups of equal depth. The reporter
induction is therefore genuine, and the pooled figures **understate** it: at the
depth where detection is most sensitive, uninjured retina shows no reporter at all.
This agrees precisely with the published immunofluorescence, which finds no EGFP in
uninjured retina, and resolves the ambient-floor ambiguity of §3.5 — the 0.56%
control value is contamination, not baseline expression.

### 3.8 Müller glia resolve into resting and activated states (Figure 4F)

Subclustering the 1,807 Müller glia (resolution 0.3, over the full gene set) gave
eight subclusters, of which three are informative.

**Table 7. Müller glia subcluster evidence.**

| | Subcluster 0 | **Subcluster 2** | Subcluster 5 |
|---|---|---|---|
| Cells | 264 | 260 | 241 |
| % EGFP⁺ | 18.56 | **27.69** | 0.83 |
| % from control | 22.0 | **1.9** | 24.9 |
| *rlbp1a* | **2.921** | 1.148 | 1.278 |
| *glula* | **4.076** | 1.727 | 2.022 |
| *glulb* | **1.323** | 0.668 | 0.413 |
| *gfap* | **1.523** | 0.857 | 0.387 |
| *pcna* | 0.292 | **1.249** | 0.322 |
| *mki67* | 0.230 | **1.640** | 0.489 |
| *ascl1a* | 0.121 | **0.356** | 0.128 |

Subcluster 2 is proliferative (*mki67* seven-fold and *pcna* 4.3-fold above
subcluster 0), reporter-active, overwhelmingly injury-derived (1.9% from control
against ~16% expected), and has **downregulated glial identity markers**.
Critically, ***ascl1a* is 2.9-fold enriched** — the master regulator of Müller glia
reprogramming in zebrafish, not a generic proliferation marker. Its enrichment
specifically in the EGFP-high, identity-low, control-depleted population identifies
this as a **reprogramming state**, not merely dividing cells.

Subcluster 0 is the reciprocal — highest identity markers, lowest proliferation,
well represented in control: **resting Müller glia**. Subcluster 5 is an internal
control: comparable size and control representation to subcluster 0 but
twenty-fold lower reporter activity, showing EGFP does not simply track cell number
or condition.

**Two states are supported, not three.** Subclusters 3 and 6 are intermediate and
described as a gradient. Subcluster 7 shows the strongest activation signature
(31.58% EGFP⁺, *ascl1a* 0.480, *mdka* 3.064) on only 19 cells and is not
interpreted. **Research question 4:** quiescent and activated states are resolved;
a separately clustering progenitor population is not, at this depth.

### 3.9 The published EGFP⁺ Müller glia signature is fully recovered

Differential expression between EGFP-positive and negative Müller glia (205 vs
1,602 cells) recovered **all 23 genes tested from the published signature, every one
significantly upregulated**:

| Gene | log2FC | p-adj | Gene | log2FC | p-adj |
|---|---|---|---|---|---|
| *mmp9* | 2.009 | 2.2×10⁻¹⁴ | *icn* | 1.226 | 2.5×10⁻¹⁰ |
| *txn* | 1.771 | 3.2×10⁻²⁰ | *cabp5a* | 1.205 | 1.7×10⁻⁹ |
| *hbegfa* | 1.652 | 1.7×10⁻¹³ | *id1* | 1.170 | 1.2×10⁻¹¹ |
| *lgals2a* | 1.468 | 5.7×10⁻¹⁵ | *mdka* | 1.112 | 2.0×10⁻⁷ |
| *glulb* | 1.367 | 1.0×10⁻¹⁶ | *sncga* | 1.057 | 1.2×10⁻¹¹ |
| *stm* | 1.335 | 1.5×10⁻⁵ | *mt2* | 1.017 | 5.7×10⁻⁷ |
| *ddr1* | 1.327 | 3.1×10⁻¹⁵ | *col18a1a* | 1.009 | 4.5×10⁻¹² |
| *crabp1a* | 1.295 | 1.5×10⁻¹² | *selenop* | 0.935 | 8.7×10⁻⁸ |
| *gfap* | 1.271 | 6.1×10⁻¹⁴ | *col15a1b* | 0.934 | 4.3×10⁻¹³ |
| *apoeb* | 1.229 | 8.0×10⁻⁹ | *pleca* | 0.895 | 1.7×10⁻⁹ |
| | | | *six3b* | 0.823 | 5.0×10⁻⁹ |
| | | | *fxyd6l* | 0.797 | 4.0×10⁻⁶ |
| | | | *cahz* | 0.693 | 2.1×10⁻⁴ |

EGFP itself ranked first (log2FC 31.07), an internal control confirming the
comparison behaves as expected.

**The total DEG count is not interpretable.** 12,750 genes reached significance,
12,475 of them upregulated (98%) — a direct consequence of the eleven-fold depth
difference documented in §3.7. This comparison is therefore presented as **targeted
validation of the published signature**, not as an independent gene count.

### 3.10 The immature/mature rod distinction does not reproduce (Figure 6)

Bise *et al.* describe two rod populations distinguished by **inverse usage of
paralogous phototransduction genes**: immature rods high in *rho*, *pde6ga*,
*guca1a*; mature rods high in *rhol*, *pde6gb*, *guca1b*. Coordinated reciprocal
usage of three paralog pairs is specific and unlikely to arise by chance, making
this a decisive test.

**Table 8. Paralog expression across rod subclusters (mean log-normalised).**

| Subcluster | *rho* | *rhol* | *pde6ga* | *pde6gb* | *guca1a* | *guca1b* |
|---|---|---|---|---|---|---|
| 0 | 6.353 | 1.122 | 3.815 | 4.698 | 1.673 | 1.286 |
| 1 | 7.631 | 1.271 | 3.928 | 5.717 | 2.908 | 1.023 |
| 2 | 6.216 | 1.687 | 4.317 | 4.735 | 1.615 | 1.876 |
| 3 | 5.618 | 1.509 | 3.804 | 4.055 | 1.805 | 1.711 |

**Two of three pairs fail to split.** *rho* exceeds *rhol* in every subcluster;
*pde6gb* exceeds *pde6ga* in every subcluster. Only *guca1a*/*guca1b* divides, and
discordantly with the others. **No subcluster carries the *rhol*-high mature
signature**: *rhol* remains 1.1–1.7 throughout against *rho* at 5.6–7.6.

Rod subclusters instead separate by cellular stress. Median depth is comparable
across subclusters 0, 1 and 2 (659, 875, 888 counts), excluding depth as the
explanation:

| | Subcluster 0 | **Subcluster 1** | Subcluster 2 | Subcluster 3 |
|---|---|---|---|---|
| Cells | 418 | 3,149 | 1,786 | 23 |
| Median counts | 659 | 875 | 888 | 232,063 |
| Median genes | 362 | **317** | 511 | 12,154 |
| % mitochondrial | 4.8 | **10.1** | 2.6 | 5.4 |
| % from control | 8.4 | **1.1** | 13.7 | 0.0 |
| *meig1* | 0.315 | **1.078** | 0.253 | 0.798 |

Subcluster 1 detects **fewer genes at equal depth** (317 vs 511) with **four-fold
higher mitochondrial content**, drawing 1.1% of cells from controls. Reduced
transcriptional complexity at matched depth with high mitochondrial load indicates
stress, not immaturity. Subcluster 3 is a 23-cell multiplet artefact.

*meig1* is four-fold elevated in subcluster 1. But the published immature rods were
defined by **low** phototransduction expression alongside high *meig1*, whereas
subcluster 1 is simultaneously high in *rho*, *pde6gb*, *gnat1*, *gngt1* **and**
*meig1*. This is not the described signature. *prom1* was absent from the annotation.

**Rod injury response.** DE at 3 dpMNU vs control gave **332 genes (179 up, 153
down)** — a balanced result, in contrast to the skewed comparisons elsewhere. Rod
sequencing depth is flat across conditions (756, 818, 806, 982 counts), which
explains the balance. The count is comparable to the published 74–182 per rod
cluster per timepoint.

GO enrichment of the 179 upregulated genes recovered **three of the four categories
Bise *et al.* report** for their rod cluster A: cytoplasmic translation (GO:0002181,
p ≈ 0), ribosome assembly (GO:0042255, p ≈ 0), and energy-coupled proton transport
(GO:0015985, p = 3.2×10⁻¹⁴, driven by seven ATP synthase subunits) — their
"proton motive force-driven ATP synthesis".

**The same signature, a different interpretation.** Bise *et al.* read elevated
translation and oxidative phosphorylation as indicating growth and biogenesis in
immature rods. The present data are equally consistent with cellular stress, in
which a collapsed transcriptome concentrates reads into abundant housekeeping
transcripts — supported here by the mitochondrial content, the reduced gene
detection at matched depth, and the failed paralog test.

### 3.11 UV cones are resolved, reproducing the published subtype split (Figure 7)

Cone subclustering (resolution 0.3) gave seven subclusters. **Subcluster 4 (285
cells, 10.9% of cones) is the UV population**, identified by five concordant markers,
each outside the range spanned by all six other subclusters:

| Marker | Subcluster 4 | Range in others |
|---|---|---|
| *opn1sw1* (UV opsin) | **5.068** | 0.703 – 2.950 |
| *arr3b* | **3.344** | 0.783 – 2.075 |
| *arr3a* | **0.830** | 2.428 – 3.910 |
| *tbx2a* | **1.088** | 0.103 – 0.590 |
| *cngb3.2* | **1.770** | 0.295 – 1.103 |
| *guca1e* | **1.484** | 0.177 – 1.161 |

**The *arr3* paralog inversion reproduced exactly.** The *arr3a*/*arr3b* ratio is
1.7–4.7 in every other subcluster but **0.25** in subcluster 4 — an order of
magnitude away and in the opposite direction, with the two genes trading values
(*arr3a* 3.3 → 0.83; *arr3b* 1.6 → 3.34). *tbx2a* is a UV cone fate transcription
factor, so its enrichment reflects the identity programme rather than
phototransduction machinery.

**UV cones decline sharply after injury:**

| | ctrl | 3dp | 7dp | 10dp |
|---|---|---|---|---|
| UV cones (n) | 109 | 106 | 18 | 52 |
| % of cones | 16.7 | 15.7 | **3.8** | 6.4 |

A four-fold decline by 7 dpMNU with partial recovery at day 10. Bise *et al.* show
by immunofluorescence that MNU destroys UV cones, with UV cone outer segments
absent at 5 dpMNU and restored only at 40 dpMNU. **These transcriptomic data
independently support a conclusion the published work established by microscopy**,
though 18 cells at 7dp is a small denominator and the compositional caveat applies.

**Cone injury response.** 15 of 16 genes named by Bise *et al.* were recovered at
3 dpMNU: all nine downregulated (*neurod1*, *rorcb*, *cry3a*, *nfil3-5*, *rimkla*,
*ipmkb*, *opn6b*, *kera*, *grk1b*) and six of seven upregulated (*rbp4l*, *guca1c*,
*rcvrn3*, *ckbb*, *gpx4b*, *arr3a*). The exception, *arr3b*, is explained by
pooling: the paper reports *arr3a* up in non-UV and *arr3b* up in UV cones, and this
test pooled all cones, 89% of which are non-UV.

**Cone DEG counts are depth-confounded and not comparable.** 674 / 326 / 3,694 at
3 / 7 / 10 dpMNU, with 96% upregulated at 10dp. Median cone depth is 1,358 (ctrl),
1,483 (3dp), 1,116 (7dp) and 1,968 (10dp) — DEG count tracks depth almost exactly.

**No claim of cone recovery is available.** The series ends at 10 dpMNU; published
restoration occurs at 40 dpMNU by immunofluorescence.

---

## 4. Comparison with the published analysis

### 4.1 Figure correspondence

Of eleven published figures, four contain single-cell data. Figures 1–4
(immunofluorescence) and 9–11 (TOR signalling, phospho-rpS6, rapamycin) report
protein localisation and pharmacological intervention and **have no transcriptomic
counterpart by construction**.

**Table 9. All eighteen reproducible panels.**

| Panel | Content | Status |
|---|---|---|
| 5C | UMAP atlas | Content reproduced; coordinates not comparable |
| 5D | Canonical marker dot plot | **Reproduced** |
| 5E | Population % per timepoint | **Reproduced** |
| 6A | Rod clusters on UMAP | Reproduced (locator panel) |
| 6B | Rod paralog dot plot | **Not reproduced** — 2 of 3 pairs fail (§3.10) |
| 6C/D | GO terms, rods | **3 of 4 categories reproduced** |
| 6E | Rod gene heatmap | Reproduced as a panel |
| 6F | Rod volcano | **Reproduced** — 332 DEGs, balanced |
| 7A | Cone clusters on UMAP | Reproduced (locator panel) |
| 7B | Opsin heatmap | **Reproduced** — UV population identified |
| 7D | Cone subtype genes | **Reproduced** — *arr3* inversion exact |
| 7E | Cone volcano | Gene-level reproduced (15/16); counts depth-confounded |
| 7F | GO terms, cones | Reproduced |
| 8A | EGFP⁺ cells on UMAP | **Reproduced** |
| 8B | EGFP⁺ per timepoint | **Reproduced** — within 0.2 pp throughout |
| 8C | EGFP⁺ % per population | **Reproduced** — MG 11.34% vs ~10% |
| 8D | Müller glia on UMAP | Reproduced (locator panel) |
| 8E | EGFP⁺/⁻ MG heatmap | **Reproduced** — 23/23 signature genes |

Note that the project guide directs comparison against "Figures 1 & 2" for the
atlas and "Figures 3 & 4" for EGFP kinetics. In the published article those are
immunofluorescence figures; the single-cell atlas is **Figure 5** and the reporter
analysis **Figure 8**.

### 4.2 Why the reporter reproduced and the rods did not

The two pipelines differ at every major processing step, yet the reporter kinetics
agree within 0.2 percentage points. This is expected on mechanistic grounds: EGFP
positivity is a near-binary determination on a transgene with no endogenous
paralog, no ambiguity of gene identity, and a strong biological signal. Such a
measurement is insensitive to normalisation and clustering choices.

The rod distinction is the opposite case, resting on **graded differences in the
relative expression of paralogous genes** — quantities directly affected by
normalisation. SCTransform treats highly expressed genes differently from
log-normalisation, and rod phototransduction transcripts are among the most
abundant in the tissue.

**The cone result provides a decisive control.** The same code, the same
subclustering resolution and the same paralog-inversion logic reproduced the
*arr3a*/*arr3b* split exactly (§3.11). The pipeline is therefore capable of
resolving photoreceptor subtype structure, and the rod failure is specific rather
than a general limitation.

The distinction between the two cases is instructive. UV cone identity is
**categorical** — a cone expresses *arr3b* or *arr3a*, and the ratio differs
ten-fold between populations. Rod maturation is **graded** — both paralogs are
expressed in every rod, and the classification depends on relative levels.
Categorical identity survives a change of normalisation; graded ratios in highly
expressed genes need not.

A fourth factor is that stress appears to dominate rod variance in these data. The
rod subclusters differ in mitochondrial content and transcriptional complexity at
matched depth, and the stressed population is almost exclusively injury-derived. If
stress is the largest axis of variation among rods, it may partition the cells
before any maturation axis can.

### 4.3 Sequencing depth as a systematic confound

Three comparisons in this analysis differ in whether the compared groups had
matched sequencing depth, and the DEG results follow that difference exactly:

| Comparison | Depth difference | DEG direction | Interpretable? |
|---|---|---|---|
| Rods, 3dp vs ctrl | none (818 vs 756) | 179 up / 153 down | **Yes** |
| Cones, 10dp vs ctrl | 1,968 vs 1,358 | 96% up | No |
| MG, EGFP⁺ vs EGFP⁻ | 15,044 vs 1,360 | 98% up | No |

The directional skew tracks the depth imbalance precisely. This is an empirical
demonstration of why cell-level differential expression without depth correction is
unreliable, and it converts the pseudobulk limitation (§4.6.3) from an assertion
into an observation.

It also motivated the stratified analysis in §3.7, which showed the reporter
finding survives depth correction — and is in fact strengthened by it.

### 4.4 Retinal ganglion cells and horizontal cells

Recovery of a discrete *rbpms2b*-positive cluster where the published analysis
reported none is unexplained. At 97 cells (0.6%), stochastic capture variation is
plausible; alternatively, clustering at resolution 0.6 may resolve a population that
partitioning at 0.2 absorbed into a neighbour. This is consistent with the general
observation that this analysis produced 20 clusters from fewer cells than the
published 17. The observation is reported without a preferred explanation.

Horizontal cells were identified only after differential expression was extended to
the full gene set. *cx52.6* is not a highly variable gene, so a marker analysis
restricted to the 2,000-gene subset could not have detected it — a concrete
illustration of why feature selection must not propagate into statistical testing.

### 4.5 An observed correlate of the capture bias

Bise *et al.* attribute increased rod representation to differential survival and
capture, inferring the mechanism from composition. The present analysis identifies a
candidate population directly: mitochondrion-dominated, gene-poor rods, near-absent
from controls and present in every injured sample.

The strength of this observation is limited by replicate variability (§3.4). The
population's abundance varies seven-fold between replicates of the same condition
and tracks per-sample mitochondrial content, so it supports the capture-bias
mechanism **qualitatively** without establishing a time course.

A methodological corollary: removing these cells on standard quality-control
grounds would discard a biologically informative population. A mitochondrial
threshold that identifies technical artefact in healthy tissue may identify a
phenotype in damaged tissue.

### 4.6 Limitations

1. **No doublet detection.** Bise *et al.* removed ~12 ± 2% of cells using scater.
   Cluster 1 (15 cells) and rod subcluster 3 (23 cells) are evident multiplets
   caught by clustering rather than filtering; others may remain within larger
   clusters.
2. **Condition and batch are confounded by design.** Each timepoint is its own pair
   of 10x runs. The variance decomposition (§2.5) shows the technical component is
   small, but it cannot be eliminated.
3. **Differential expression is descriptive, not inferential.** Cell-level Wilcoxon
   tests treat cells as independent, but cells within a replicate are correlated and
   only two replicates exist per condition. Pseudobulk across the eight samples
   would be the correct design. §4.3 demonstrates the practical consequence.
4. **Unbalanced recovery.** Controls contributed fewest cells (3,271). Replicate 1
   of each injured pair consistently carries higher mitochondrial content,
   indicating a systematic processing difference between runs.
5. **No recovery claim is possible.** The series ends at 10 dpMNU; published
   restoration occurs at 30–40 dpMNU by immunofluorescence.
6. **Transcript is not protein.** EGFP mRNA need not track reporter protein, which
   persists on a different timescale.
7. **2.3% of cells remain unresolved**, and *prom1* was absent from the annotation,
   preventing one component of the rod comparison.
8. **GO enrichment is not methodologically equivalent.** Enrichr uses different
   background gene sets and versions from topGO, so term lists are not directly
   comparable even given identical input.

---

## 5. Conclusions

An independent Python/Scanpy reanalysis reproduces the central single-cell
conclusion of Bise *et al.* (2023). *careg:EGFP* is induced after MNU photoreceptor
ablation, peaks at 3 dpMNU, declines without returning to baseline by 10 dpMNU, and
marks Müller glia specifically (11.34% positive, 3.75-fold enriched, 42.6% of all
reporter-positive cells). Agreement within 0.2 percentage points at every timepoint,
obtained through different normalisation, clustering and annotation, establishes
these kinetics as robust to analytical choice.

The result also survives its most serious potential confound. EGFP-positive Müller
glia are sequenced eleven-fold more deeply than negative cells, but within matched
depth quartiles the induction is stronger, not weaker — rising from 0% of control
cells to 54.3% at 3 dpMNU.

All 23 tested genes from the published EGFP-positive Müller glia signature were
recovered as significantly upregulated. Subclustering resolved resting and activated
states, the latter reporter-active, *ascl1a*-positive, proliferative and depleted of
control cells — the expected dedifferentiation signature.

The published immature/mature rod distinction did not reproduce; the paralog
inversion failed for two of three pairs, and rod heterogeneity is instead organised
by cellular stress. That the equivalent test in cones reproduced exactly identifies
this as specific to the rod comparison rather than a limitation of the pipeline, and
suggests the distinction between categorical and graded transcriptional differences
determines which findings survive a change of method.

---

## 6. Bibliography

Bernardos RL, Barthel LK, Meyers JR, Raymond PA (2007). Late-stage neuronal
progenitors in the retina are radial Müller glia that function as retinal stem
cells. *Journal of Neuroscience* 27:7028–7040.

Bise T, Pfefferli C, Bonvin M, Taylor L, Lischer HEL, Bruggmann R, Jaźwińska A
(2023). The regeneration-responsive element *careg* monitors activation of Müller
glia after MNU-induced damage of photoreceptors in the zebrafish retina.
*Frontiers in Molecular Neuroscience* 16:1160707. doi:10.3389/fnmol.2023.1160707

Fischer AJ, Reh TA (2001). Müller glia are a potential source of neural
regeneration in the postnatal chicken retina. *Nature Neuroscience* 4:247–252.

Goldman D (2014). Müller glial cell reprogramming and retina regeneration. *Nature
Reviews Neuroscience* 15:431–442.

Hafemeister C, Satija R (2019). Normalization and variance stabilization of
single-cell RNA-seq data using regularized negative binomial regression. *Genome
Biology* 20:296.

Hoang T, Wang J, Boyd P, Wang F, Santiago C, Jiang L, *et al.* (2020). Cross-species
transcriptomic and epigenomic analysis reveals key regulators of injury response and
neuronal regeneration in vertebrate retinas. *bioRxiv* [Preprint].
doi:10.1101/717876 **[VERIFY — a peer-reviewed version was published; cite that if
you consulted it.]**

Macosko EZ, Basu A, Satija R, Nemesh J, Shekhar K, Goldman M, *et al.* (2015).
Highly parallel genome-wide expression profiling of individual cells using nanoliter
droplets. *Cell* 161:1202–1214.

Maurer E, Tschopp M, Tappeiner C, Sallin P, Jaźwińska A, Enzmann V (2014).
Methylnitrosourea (MNU)-induced retinal degeneration and regeneration in the
zebrafish: histological and functional characteristics. *Journal of Visualized
Experiments* 92:e51909.

Nagashima M, Barthel LK, Raymond PA (2013). A self-renewing division of zebrafish
Müller glial cells generates neuronal progenitors that require N-cadherin to
regenerate retinal neurons. *Development* 140:4510–4521.

Ogawa Y, Corbo JC (2021). Partitioning of gene expression among zebrafish
photoreceptor subtypes. *Scientific Reports* 11:17340.

Pfefferli C, Jaźwińska A (2017). The *careg* element reveals a common regulation of
regeneration in the zebrafish myocardium and fin. *Nature Communications* 8:15151.

Raymond PA, Barthel LK (2004). A moving wave patterns the cone photoreceptor mosaic
array in the zebrafish retina. *International Journal of Developmental Biology*
48:935–945.

Stuart T, Butler A, Hoffman P, Hafemeister C, Papalexi E, Mauck WM III, *et al.*
(2019). Comprehensive integration of single-cell data. *Cell* 177:1888–1902.

Tappeiner C, Balmer J, Iglicki M, Schuerch K, Jaźwińska A, Enzmann V, *et al.*
(2013). Characteristics of rod regeneration in a novel zebrafish retinal
degeneration model using N-methyl-N-nitrosourea (MNU). *PLoS ONE* 8:e71064.

Traag VA, Waltman L, van Eck NJ (2019). From Louvain to Leiden: guaranteeing
well-connected communities. *Scientific Reports* 9:5233.

Wan J, Goldman D (2016). Retina regeneration in zebrafish. *Current Opinion in
Genetics and Development* 40:41–47.

Wolf FA, Angerer P, Theis FJ (2018). SCANPY: large-scale single-cell gene expression
data analysis. *Genome Biology* 19:15.

Zhang Z, Shen X, Gude DR, Wilkinson BM, Justice MJ, Flickinger CJ, *et al.* (2009).
MEIG1 is essential for spermiogenesis in mice. *Proceedings of the National Academy
of Sciences USA* 106:17055–17060.

**Software citations to add if cited in Methods** **[VERIFY]** — confirm volume and
pages: NumPy (Harris *et al.*, 2020, *Nature* 585); SciPy (Virtanen *et al.*, 2020,
*Nature Methods* 17); Matplotlib (Hunter, 2007, *CiSE* 9); UMAP (McInnes *et al.*,
2018, arXiv:1802.03426); gseapy (Fang *et al.*, 2023, *Bioinformatics*); anndata
(Virshup *et al.*, *JOSS*).

*All references except those marked* **[VERIFY]** *appear in the reference list of
Bise* et al. *(2023) or are standard method citations. Confirm each against the
original source.*

---

## Appendix A — AI Usage Disclosure

**To be completed by the author.** Template: `docs/ai_usage_disclosure.md`.

### A.1 Tools and models

| Tool | Model | Applied to | Dates |
|---|---|---|---|
| *[fill in]* | | | |

### A.2 Prompt log

Record prompts verbatim; store long prompts under `docs/prompts/`.

| # | Prompt | Purpose | Outcome and modifications |
|---|---|---|---|
| 1 | | | |

### A.3 Division of authorship

The following are the author's own scientific judgement: all cluster annotations and
rationales in `tables/cluster_annotations.csv`; the decision to reject panel scores
in favour of marker-based annotation; the interpretation of every result.

### A.4 Errors identified and corrected

Documenting caught errors evidences critical rather than passive use. Errors
identified during this project:

- **Statistical tests confined to 2,000 of 22,813 genes.** Code selecting between
  the log-normalised layer and `.raw` tested only whether a layer existed; after
  feature selection that layer persisted in reduced form, so all differential
  expression and marker discovery ran on 8.8% of the transcriptome. Detected by
  investigating why GO enrichment returned no terms — only 10 genes had passed the
  significance filter. Correction changed a substantive conclusion: horizontal cells
  (712 cells, 4.5% of the dataset) were identified only after the fix, since
  *cx52.6* is not a highly variable gene.
- **Marker gene availability checked against the wrong matrix**, silently excluding
  *meig1*, *rom1b* and *guca1b* from the rod evidence table — the genes required for
  the paralog inversion test. Detected by comparing a dot plot against a table that
  should have contained the same genes.
- **Panel scores assigned rod identity to oligodendrocyte and mesenchymal clusters**,
  an artefact of ambient rod mRNA. Scores rejected in favour of marker rankings.
- **Invalid organism parameter** in the enrichment call (`'Fish'`), rejected by
  gseapy 1.3.1, which requires `'zebrafish'` or `'fish'`.
- **Misnumbered figure references in the project guide**, which directs comparison
  against Figures 1–4; those are immunofluorescence figures.
- **Cell-count arithmetic conflating pre- and post-QC totals** in an interim summary.

### A.5 Verification performed

- [ ] Confirmed every marker gene exists in the zebrafish annotation (ZFIN/Ensembl)
- [ ] Verified *cx52.6* as a horizontal cell marker in the primary literature
- [ ] Executed every notebook and inspected the outputs
- [ ] Confirmed `pct_counts_mt` is non-zero
- [ ] Verified the EGFP feature name in `features.tsv.gz`
- [ ] Checked all published values against the article, not a summary of it
- [ ] Re-derived at least one statistic manually from the raw matrix
- [ ] Can explain every line of the pipeline without assistance

### A.6 Statement

*[Signed statement that the analytical decisions, annotations and interpretations
are the author's own, that AI assistance is disclosed in full, and that the accuracy
of the code and biological claims has been verified.]*

**Mahabubur Rahman** — *[date]*

---

## Appendix B — Reproducibility

Repository: https://github.com/mr-mahabubur-rahman/sc-project1-zebrafish-retina
Tag: `v1.0-peerreview`

Eight notebooks execute in sequence, each loading a checkpoint rather than
inheriting kernel state. Versions in `docs/session_info.json`, pinned in
`requirements.txt` and `uv.lock`. Seed 0 throughout; UMAP coordinates may differ
across platforms while cluster membership is stable.

Data are not committed (GEO GSE202212); `data/README.md` documents the expected
layout. `docs/method_comparison.md` compares the published method with this
implementation step by step;
`tables/paper_figure_reproduction_status.csv` records the status of each panel.
