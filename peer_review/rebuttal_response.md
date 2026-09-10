# Response to peer review — Group 4

**Project 1 — Zebrafish retina regeneration**
Repository: https://github.com/mr-mahabubur-rahman/sc-project1-zebrafish-retina
Review received: 3 September 2026 · Release: `v2.0-final`

---

We thank Group 4 for a careful and constructive review. We are particularly
grateful that the reviewers cloned the repository, built the environment with
`uv`, and executed the full workflow rather than reading the code alone.
Independent confirmation that the pipeline runs cleanly on another machine is the
most useful thing a reviewer can provide, and the reported end-to-end run of 2
minutes 14 seconds gives us confidence the environment specification is correct.

The review raised no major methodological concerns. Every minor comment is
addressed below in the reviewers' own numbering. One comment we have answered with
a question rather than a change, and we say so plainly where that is the case.

---

## 1. GitHub repository, code quality and reproducibility

### 1a — Contributors, and README title

> *"Mentioning the contributors with their email IDs and affiliations will make it
> more professional. Moreover, the README title is not explanatory itself and does
> not match the actual report title."*

**Accepted on both counts.** The README now lists all nine group members with role,
affiliation and email address. Its title has been changed to the full report title,
so the repository and the manuscript are unambiguously the same piece of work.

### 1b — Overview to orient the reader

> *"An overview in the README section is needed to hook the reader's interest."*

**Accepted.** The README now opens with a short statement of the scientific
question and why a reanalysis was undertaken, followed by a headline results table
giving the reproduced *careg:EGFP* values against the published ones, the number of
cells and cell types recovered, and the panel reproduction score. A reader sees
what the project found within the first screen.

### 1c — Repository structure

> *"The repository structure can be cleaner. Please add the 'tables' in the
> 'results' folder. The docs folder contains too many files. Make a single final
> report... Otherwise, it is very confusing which report to review."*

**Accepted, and addressed fully.** This was the reviewers' only "No" in the
reproducibility check. The repository now follows the layout in the peer-review
guidelines:

- `results/tables/` — all CSV outputs, previously at `tables/`
- `results/figures/` — all figures, previously at `figures/`
- `report/` — a single final report as PDF and DOCX, with the AI disclosure
  appendix inside the same document
- `peer_review/` — this response and the review received
- `docs/` — supporting documentation only: method comparison, validation
  checklist, git workflow, session information and the prompt records

Three superseded drafts were deleted rather than archived, so there is no ambiguity
about which document is current. Paths in `scripts/config.py` were updated and the
pipeline re-run end to end to confirm the move broke nothing.

### 1d — Direct link to the data

> *"In the data README section, a link to the GEO data will help others to
> reproduce the project more easily."*

**Accepted.** `data/README.md` now links directly to the GEO accession, maps each
of the eight sample directories to its condition and replicate, and gives the
expected cell counts per sample so that a misplaced file is caught immediately
rather than surfacing later. A one-line pre-flight check confirms the directory
layout before any notebook is run.

---

## 2. Quality control and analytical integrity

### 2a — Whether regression was performed

> *"The report does not clearly specify if the regression is performed or not. If
> performed, there is no clear description of how many technical covariates are
> regressed out."*

**Accepted as a clarity problem.** Regression was performed, on **two** technical
covariates — `total_counts` and `pct_counts_mt` — using `sc.pp.regress_out` after
subsetting to highly variable genes. This was stated in §2.4 but embedded in a
longer sentence about the wider preprocessing sequence, which made it easy to miss.

It now has its own paragraph, names both covariates explicitly, and states why
regression runs after HVG subsetting rather than on the full matrix: the covariates
and the result on the retained genes are identical, but regressing 22,813 genes
across 15,893 cells costs hours of compute for genes discarded before PCA. The
deviation was already recorded in `cfg.REGRESS_ON_HVG_SUBSET`; it is now visible in
the manuscript as well.

---

## 3. Biological interpretation and literature context

### 3a — log2 fold-change threshold

> *"The log2F change threshold is not clearly mentioned, whether it is performed or
> not."*

**Accepted as a placement problem.** Differential expression used the Wilcoxon
rank-sum test with Benjamini-Hochberg correction, reporting genes at adjusted
p < 0.05 and |log2FC| ≥ 0.25.

The threshold appeared in §3.10 and in the volcano figure captions but not in
Methods, which is where a reader looks for it. It is now stated in §2.8 alongside
the test, and both cutoffs are defined as named constants in `scripts/config.py`
(`DE_PADJ_CUTOFF`, `DE_LOG2FC_CUTOFF`) so that no figure or table can apply a
different one. We also note there that the published analysis reports adjusted p
alone with no fold-change floor, which matters when DEG counts are compared.

---

## 4. Visualisation, figure quality and presentation

### 4a — Colour-blind accessibility

> *"Colour-blind accessibility is not well defined."*

**Accepted, and testing found a real defect.** This was not addressed anywhere in
the draft. Rather than assert accessibility, we measured it: under simulated
protanopia the previous palette rendered **rods and retinal ganglion cells 2.3
units apart in CAM02-UCS — effectively the same colour**, and both populations
appear on the same UMAP. Fourteen category pairs fell below the practical
separability threshold of 10 across the three deficiency forms.

The palette has been replaced with one drawn from colour-vision-deficiency-safe
schemes (Okabe-Ito; Tol, 2021). The assignment of colour to cell type was not
chosen by eye but optimised to maximise the smallest perceptual distance between
any two categories, evaluated simultaneously under normal vision, deuteranopia,
protanopia and tritanopia. The minimum distance is now **8.8**, no pair falls below
5, and pairs below 10 are reduced from 14 to 7. Every figure in the report was
regenerated with the new palette.

The verification is reproducible rather than asserted:

```
uv run python scripts/colour_check.py
```

It prints the before-and-after comparison and writes Supplementary Figure S3, which
shows the palette under each deficiency form. A paragraph in §2.2 records the
method and states honestly that thirteen categories cannot be separated by hue
alone under any palette — which is why every figure also identifies populations by
direct label or by a legend ordered to match the plot.

### 4b, 4c — Writing style and density

> *"It is filled with analysis without a clean transition from one section to
> another"* and *"it has become very complicated to understand the whole context."*

**Accepted in part.** We have added an orienting paragraph at the head of the
Results and of the Discussion, each stating what the section covers and in what
order, so a reader knows the shape of the argument before entering the detail.

We have not reduced the analytical content, and we would like to explain why. The
central claim of this report is a negative one — that the published
immature/mature rod distinction does not reproduce. A negative result is only
interpretable if the alternative explanations have been tested and excluded, which
is what §4.7 and §4.8 do. Removing that material would make the report shorter and
considerably weaker. We have instead tried to make its structure easier to follow.

### 4d — Biological narrative and flow chart

> *"The results are more of computational output rather than a connection with the
> biological story. A flow chart can provide a better understanding."*

**Accepted.** A pipeline schematic has been added to Methods showing the path from
the eight 10x matrices through quality control, normalisation, clustering and
annotation to the four analysis branches, with the notebook and the key quantity at
each stage. It is generated by `scripts/report_figures.py` rather than drawn by
hand, so it uses the same palette and settings as every other figure and can be
regenerated by anyone with the repository.

The Results section now opens with a paragraph framing the biological narrative:
photoreceptor ablation, the glial response, reporter induction as a marker of that
response, and the fate of the photoreceptors themselves.

### 4e — Graphical abstract

> *"Adding a graphical abstract will hook the reader."*

**Accepted.** A graphical abstract now follows the abstract, in four panels: the
injury model, the reporter time course against the published values, the thirteen
populations recovered, and a summary of which published conclusions survive an
independent analysis, which do not, and what this reanalysis adds. It is generated
by the same script as the pipeline schematic.

### 4f — Purpose of the reanalysis

> *"In the introduction part, the term descriptions are good, but there is a gap in
> the purpose of this reanalysis."*

**Accepted as a placement problem.** §1.4 sets out the rationale: re-executing
published code establishes only that a pipeline is deterministic, whereas
recovering the same conclusion through independent methodological choices tests
whether a finding survives the decisions every analyst must make. The second is the
stronger claim and the one attempted here.

We agree this arrived too late in the introduction. A sentence stating the purpose
has been added at the end of §1.1, so a reader knows why the reanalysis exists
before reading the biological background.

### 4g — Limitations

> *"The limitations are vaguely written, it should be clearly mentioned."*

**We would welcome more specific guidance here.** §4.6 lists eight limitations,
each naming a concrete methodological gap and the conclusion it bounds:

1. No doublet detection, against the ~12 ± 2% removed in the published analysis,
   with the two multiplet clusters we identified named
2. Condition and sequencing batch confounded by design, quantified in §2.5 and
   tested directly in §4.7
3. Cell-level differential expression overstating significance with n = 2
   replicates, with the practical consequence demonstrated in §4.3
4. Unbalanced recovery, with the control cell count and the replicate-1
   mitochondrial pattern given
5. No recovery claim possible, with the published restoration timepoints stated
6. Transcript abundance is not protein abundance
7. 2.3% of cells unresolved, and *prom1* absent from the annotation
8. GO enrichment not methodologically equivalent to the published topGO analysis

We have added a line at the head of the section stating which conclusions each
limitation bounds, which we hope addresses the concern. If a particular limitation
reads as vague to the reviewers, we would be glad to revise it — but we did not
want to dilute a section that already names specific gaps without knowing which one
was meant.

---

## Summary of changes

| Comment | Change |
|---|---|
| 1a | README contributors table with affiliations; title aligned with the report |
| 1b | README overview and headline results table |
| 1c | `results/tables/`, `results/figures/`, single report in `report/`, `peer_review/` added, `docs/` cleared of superseded drafts |
| 1d | Direct GEO link, per-sample mapping and expected cell counts in `data/README.md` |
| 2a | Regression stated in its own paragraph in §2.4 with both covariates named |
| 3a | padj < 0.05 and \|log2FC\| ≥ 0.25 added to Methods §2.8 and to `config.py` as named constants |
| 4a | Palette replaced and verified; minimum perceptual distance 2.3 → 8.8; all figures regenerated; Supplementary Figure S3 added; `scripts/colour_check.py` makes it reproducible |
| 4b, 4c | Orienting paragraphs added to Results and Discussion; analytical content retained, with reasons given |
| 4d | Pipeline schematic added to Methods; biological framing paragraph opens Results |
| 4e | Graphical abstract added after the abstract |
| 4f | Purpose of the reanalysis stated at the end of §1.1 |
| 4g | Summary line added at the head of §4.6; clarification requested |

The revised report is at `report/Rahman_Project1_Report.pdf` and the tagged release
is `v2.0-final`.

We thank the reviewers again for their time and for the evident care taken over the
review.

**Mahabubur Rahman**, on behalf of Group 1
