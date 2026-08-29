# data/

This directory is **git-ignored**. The eight 10x sample folders are not committed.

Expected layout:

```
data/
├── ctrl1/filtered_feature_bc_matrix/{barcodes,features,matrix}.{tsv,mtx}.gz
├── ctrl1/web_summary.html
├── ctrl2/ ...
├── 3dp1/ ...   3dp2/ ...
├── 7dp1/ ...   7dp2/ ...
└── 10dp1/ ...  10dp2/ ...
```

Source: the project data folder supplied with the assignment. The underlying
sequencing data are deposited at GEO accession **GSE202212** (Bise et al. 2023).

`scripts/io_utils.discover_samples()` finds these folders automatically; no path
is hard-coded anywhere in the workflow.
