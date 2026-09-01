# Prompts submitted — selected excerpts

**This is a selection of 24 prompts from the 149 submitted during the working
session.** It is not the complete record. The prompts are reproduced verbatim and
in chronological order, and the number given for each is its position in the full
session.

- **Tool:** Claude Opus 5
- **Operator:** Mahabubur Rahman
- **Period:** 29–31 August 2026, to 06:10 UTC where the session log was captured
- **Opened with:** the prompt in `00_master_prompt.md` (reproduced in full in this directory)
- **Prompts submitted in the session:** 149
- **Prompts reproduced here:** 24

## Basis of selection

These were chosen as the prompts that most directly show how the analysis was
directed and checked: points where an output was questioned, an error reported, a
claim rejected, or a procedure confirmed before work was overwritten. They span
cell type annotation, figure verification, the paper-reproduction exercise, file
and kernel handling, and the GO enrichment failure.

The prompts not reproduced here are largely of three kinds: short procedural
questions ("where do I save this", "what next"), pasted console output supplied so
that an error could be diagnosed, and exchanges concerning report drafting. The
median prompt across the whole session was 190 characters and 65 of 149 were under
120, so the omitted material is mostly brief.

## Scope of the session log

Work after 31 August 06:10 UTC — the integration and resolution sensitivity
analyses (steps 08 and 09) and the report revisions incorporating them — continued
in the same session and is not covered by the log these prompts were taken from.
The record for that period is the notebooks, scripts and tables committed to the
repository.

## Format

Each prompt is reproduced between fenced code blocks so that whitespace and pasted
output are preserved as submitted. Nothing within a prompt has been edited or
shortened. A short note records what each achieved.

---

## Excerpt 1  ·  prompt 33 of 149

```
do I need to run the code : annotation_table = annotation.build_annotation_table(adata, annotations)
io_utils.save_table(annotation_table, "cluster_annotations.csv")
display(annotation_table)
adata = annotation.apply_annotations(adata, annotation_table) again or to run next code:
```

*Checked whether a cell needed re-running before proceeding rather than assuming state.*

## Excerpt 2  ·  prompt 38 of 149

```
I dont run the code for: Figure 3C — annotated UMAP
fig = plotting.figure_03_annotated_umap(adata)
4. Cell-type proportions (research question 2)
```

*Noticed a figure cell had not been executed.*

## Excerpt 3  ·  prompt 42 of 149

```
Did I get Figure 5: Quantification of Cell-type proportion shifts during regeneration?
```

*Verified that a required project figure had actually been produced.*

## Excerpt 4  ·  prompt 67 of 149

```
yes I need draft report. But before that Compare my tables and Figures side-by-side with the original paper and our reproduced figures
```

*Requested a side-by-side comparison with the published figures before any report was drafted.*

## Excerpt 5  ·  prompt 68 of 149

```
Did I reproduce any of figures of the main paper?
```

*Asked directly what had and had not been reproduced.*

## Excerpt 6  ·  prompt 69 of 149

```
So we can't reproduce Figures 1–4 and 9–11 are immunofluorescence, TOR signalling and rapamycin experiments. No scRNA-seq counterpart exists. Right. Do we reproduce others figures from our data?
```

*Established which published figures have no transcriptomic counterpart and why.*

## Excerpt 7  ·  prompt 70 of 149

```
You said, Figure 5 (design, UMAP atlas, marker dot plot, cluster percentages) — your Figures 2, 3 and 5 cover panels C, D and E. Done. But I can't match with them
```

*Rejected an over-stated claim of correspondence; the panels did not in fact match.*

## Excerpt 8  ·  prompt 71 of 149

```
Except Figures 1–4 and 9–11, Do we reproduce?
```

*Pressed for a complete accounting of the remaining panels.*

## Excerpt 9  ·  prompt 74 of 149

```
Before writing the report, we need to reproduce all the figures except Figures 1–4 and 9–11 are immunofluorescence, TOR signalling and rapamycin experiments. Do you rewrite the notebooks that I can run in jupyter lab?
```

*Directed that the reproducible panels be produced before writing began.*

## Excerpt 10  ·  prompt 75 of 149

```
one question: where are the notebooks for the figures 5a, 5B, 6A, 7A, 7C, 8D?
```

*Identified panels for which no notebook existed.*

## Excerpt 11  ·  prompt 77 of 149

```
how to run paper figures, Step 07 paper figure reproduction, Build notebooks?
```

*Asked how to execute the new reproduction step.*

## Excerpt 12  ·  prompt 99 of 149

```
In this figure, the position of figure name placed wrongly
```

*Caught a title overlapping the plot content.*

## Excerpt 13  ·  prompt 105 of 149

```
do I replace paper figures? or will both old and new exist
```

*Checked whether replacing files would leave stale outputs behind.*

## Excerpt 14  ·  prompt 107 of 149

```
I want the Fig 5D equivalent | canonical markers either top or below the figure without overlapping.
```

*Specified the required placement of a figure title.*

## Excerpt 15  ·  prompt 108 of 149

```
and same problems in step_05.  I wnat to fix this also
```

*Identified the same fault in a second notebook.*

## Excerpt 16  ·  prompt 112 of 149

```
Sorry I confused: here is step_05:
```

*Supplied the actual notebook after a mismatch was found.*

## Excerpt 17  ·  prompt 113 of 149

```
Do I need download it. Before download, Do I need to run in old step_05: Kernel → Restart Kernel and Clear All Outputs
```

*Checked whether a kernel restart was needed before replacing a notebook.*

## Excerpt 18  ·  prompt 114 of 149

```
Do I replace step_05
```

*Confirmed the replacement procedure before overwriting work.*

## Excerpt 19  ·  prompt 116 of 149

```
Before running 06 and 07, Do I need to follow: Kernel → Restart Kernel and Clear All Outputs then Run All:
```

*Confirmed the restart procedure before re-running downstream steps.*

## Excerpt 20  ·  prompt 118 of 149

```
where Do I insert: annotation_table = annotation.build_annotation_table(adata, annotations)
```

*Asked where a cell belonged rather than pasting it arbitrarily.*

## Excerpt 21  ·  prompt 120 of 149

```
We can fix NameError                                 Traceback (most recent call last) Cell In[14], line 3       1 # Confidence audit: how many cells rest on weak calls?       2 audit = ( ----> 3     annotation_table.groupby("confidence")["n_cells"].agg(["count", "sum"])       4     .rename(columns={"count": "n_clusters", "sum": "n_cells"})       5 )       6 audit["pct_of_dataset"] = (100 * audit["n_cells"] / adata.n_obs).round(1)  NameError: name 'annotation_table' is not defined
```

*Reported a NameError caused by cells being ordered wrongly in the notebook.*

## Excerpt 22  ·  prompt 121 of 149

```
here is this:     "6":  {"cell_type": "Muller glia", "positive_markers": ["rlbp1a", "glula", "mdka", "crabp1a", "fabp7a"],
           "conflicting_markers": [], "confidence": "high",
           "rationale": "MG identity (rlbp1a, glula). Present at 5.5% in controls, so this is the "
                        "main MG population, not an activated subset. mdka and crabp1a are elevated "
                        "and both appear in the paper's EGFP-positive MG signature, but activation "
                        "state is assessed in step 06 from EGFP and proliferation markers."},
```

*Supplied the annotation entry under discussion so the correction was made against the real text.*

## Excerpt 23  ·  prompt 122 of 149

```
without changing  NameError: name 'annotation_table' is not defined. can I continue to run next code?
```

*Asked whether the analysis could proceed with the error unresolved.*

## Excerpt 24  ·  prompt 132 of 149

```
last night night you said that GO enrichment fails gracefully. It uses gseapy's Enrichr with `organism='Fish'`, which needs network access. If it can't reach Enrichr it says so and continues. Note this is not the paper's method — they used topGO with a custom background, so term lists will differ even on identical input. Say that in Methods.. but in this case I think we did not GO enrichment, why?
```

*Questioned why GO enrichment had not run; this uncovered an invalid organism parameter that had silently disabled it.*
