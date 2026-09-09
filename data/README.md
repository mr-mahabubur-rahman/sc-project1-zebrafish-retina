# data/

The eight 10x sample folders are **not committed** to this repository. Count
matrices are large, and GitHub is not a data archive. This file explains how to
obtain them.

## Source

Bise T, Pfefferli C, Bonvin M, Taylor L, Lischer HEL, Bruggmann R, Jaźwińska A
(2023). *Front. Mol. Neurosci.* 16:1160707.

**GEO accession GSE202212:**
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202212

The reference used for alignment was GRCz11 with the *careg:EGFP* transgene added,
so the feature table contains an `EGFP` entry alongside the endogenous
transcriptome. `scripts/io_utils.find_egfp_feature()` confirms its presence before
any analysis proceeds rather than assuming it.

## Getting the data

Open the accession page above and download the supplementary files from the
**Supplementary file** section at the foot of the page. GEO offers them either
individually or as a single archive via the *(http)* and *(custom)* links.

Each sample provides the three files of a standard 10x filtered feature-barcode
matrix. Unpack them so that each sample sits in its own directory named as below.
The sample names matter: `scripts/io_utils.discover_samples()` matches on them, and
`cfg.EXPECTED_SAMPLES` asserts that all eight are present.

| Directory | Condition | Replicate |
|---|---|---|
| `ctrl1/` | Uninjured control (heat-inactivated MNU) | 1 |
| `ctrl2/` | Uninjured control (heat-inactivated MNU) | 2 |
| `3dp1/` | 3 days post-MNU | 1 |
| `3dp2/` | 3 days post-MNU | 2 |
| `7dp1/` | 7 days post-MNU | 1 |
| `7dp2/` | 7 days post-MNU | 2 |
| `10dp1/` | 10 days post-MNU | 1 |
| `10dp2/` | 10 days post-MNU | 2 |

GEO's own file names differ from these directory names. Map each downloaded file
to the correct directory using the sample titles shown on the accession page.

## Expected layout

```
data/
  ctrl1/
    filtered_feature_bc_matrix/
      barcodes.tsv.gz
      features.tsv.gz
      matrix.mtx.gz
    web_summary.html          optional, not read by the pipeline
  ctrl2/
    filtered_feature_bc_matrix/ ...
  3dp1/    3dp2/
  7dp1/    7dp2/
  10dp1/   10dp2/
```

The three `.gz` files inside `filtered_feature_bc_matrix/` are required and are
checked by name; `cfg.REQUIRED_10X_FILES` lists them. Leave them compressed —
Scanpy reads them as they are.

## Checking the layout before running

Step 01 validates the directory structure and reports what it finds, but you can
check first:

```bash
uv run python -c "from scripts import io_utils; print(io_utils.discover_samples())"
```

That should list all eight sample directories. If a sample is missing or a file is
misplaced, the error names the specific path rather than failing later inside the
loading code.

## What you should see after loading

Step 01 reports these cell counts, which are a useful check that the right files
landed in the right directories:

| Sample | Cells loaded | Cells after QC |
|---|---|---|
| ctrl1 | 1,155 | 906 |
| ctrl2 | 2,116 | 1,609 |
| 3dp1 | 2,873 | 2,317 |
| 3dp2 | 2,388 | 1,844 |
| 7dp1 | 3,444 | 2,562 |
| 7dp2 | 1,480 | 1,165 |
| 10dp1 | 4,124 | 3,333 |
| 10dp2 | 2,517 | 2,157 |
| **Total** | **20,097** | **15,893** |

25,433 features are present before filtering, including `EGFP`; 22,813 genes are
retained after removing those detected in fewer than three cells.

## Why this directory is git-ignored

`.gitignore` excludes `data/*` with the exception of this file. The matrices are
several hundred megabytes and are redistributable from GEO, so committing them
would bloat the repository without making the analysis any more reproducible. The
`.h5ad` checkpoints in `results/` are excluded for the same reason: they are
regenerable by re-running the notebooks.
