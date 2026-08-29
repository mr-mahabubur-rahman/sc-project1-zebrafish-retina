# Appendix — AI Usage Disclosure

*Template. Fill in your own entries. Delete any row you did not actually do, and
add the ones you did. This appendix is graded on honesty and on evidence of
verification, not on how little AI you used.*

## 1. Tools used

| Tool | Model | Where used | Date range |
|---|---|---|---|
| Claude (Anthropic) | *[model name shown in the interface]* | Pipeline architecture, script authoring, notebook scaffolding | *[dates]* |
| *[add others]* | | | |

## 2. Prompt log

The course policy requires the **exact** prompts. Keep them verbatim, including
the long one; paraphrasing defeats the purpose.

| # | Prompt (verbatim) | Purpose | What I did with the output |
|---|---|---|---|
| 1 | *[paste the master prompt used to generate the workflow]* | Generate the reproducible Scanpy pipeline: repository structure, `scripts/` modules, seven notebooks, figure and table specifications | Reviewed every module; ran the pipeline on the real data; corrected *[list what you changed]* |
| 2 | | | |
| 3 | | | |

Suggested storage: keep the full prompts in `docs/prompts/` as plain text files
and reference them here by filename if the table becomes unwieldy.

## 3. What the AI produced, and what I changed

| Component | AI-generated | My modifications |
|---|---|---|
| `scripts/config.py` | yes | *[parameters you changed and why]* |
| `scripts/io_utils.py` | yes | |
| `scripts/qc.py` | yes | |
| `scripts/preprocessing.py` | yes | |
| `scripts/clustering.py` | yes | |
| `scripts/annotation.py` | yes | |
| `scripts/egfp_analysis.py` | yes | |
| `scripts/plotting.py` | yes | |
| Notebooks (structure and markdown) | yes | |
| **Cluster annotations** | **no — my own scientific judgement** | the annotation dictionary in Step 05 and every rationale string |
| **Biological interpretation in the report** | **no** | |
| Report text | *[state honestly]* | |

## 4. Verification I personally performed

Tick only what you actually did.

- [ ] Confirmed every marker gene named in `config.py` exists in the zebrafish
      annotation, using ZFIN / Ensembl rather than trusting the model.
- [ ] Confirmed the guide's marker panel genes against the primary literature
      (Hoang et al. 2020; Ogawa & Corbo 2021; Bise et al. 2023).
- [ ] Ran every notebook end to end and inspected the outputs, rather than
      assuming the code was correct because it did not error.
- [ ] Checked that `pct_counts_mt` is non-zero, i.e. the `mt-` prefix matched real
      genes in this annotation.
- [ ] Verified the EGFP feature name against `features.tsv.gz` directly.
- [ ] Cross-checked the paper's reported figures against the article itself, not
      against the AI's summary of it. *(One specific reason: the project guide's
      own figure references are misnumbered — the scRNA-seq atlas is Figure 5 and
      the EGFP analysis is Figure 8.)*
- [ ] Confirmed each cluster annotation against multiple concordant markers.
- [ ] Re-derived at least one statistic by hand to check the pipeline
      (e.g. EGFP⁺ percentage in one sample from the raw matrix).
- [ ] Can explain every line of the scripts and notebooks without assistance.

## 5. Errors or questionable output I caught

*Recording these strengthens the appendix — it is evidence of critical use rather
than copy-pasting.*

| What the AI produced | Why it was wrong or doubtful | What I did |
|---|---|---|
| | | |

## 6. Statement

The analysis decisions, cluster annotations, and biological interpretations in
this report are my own. AI assistance was used for the components listed above and
is disclosed in full. I have verified the accuracy of the code and the biological
claims, and I can explain and defend every step of the workflow.

Signed: **Mahabubur Rahman** · Date: *[date]*
