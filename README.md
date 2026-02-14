# AZFc CNV Annotation and CNVkit Post-processing Pipeline

Utilities for processing CNVkit `.cnr` outputs into integer copy-number
estimates across AZFc amplicons and annotating known AZFc structural
variants, including gr/gr deletion signatures.

This pipeline was developed to support the analysis described in:

**Association of Y-chromosomal gr/gr deletions with testicular germ cell
tumour: whole-genome analysis of 198,306 individuals**

------------------------------------------------------------------------

## Overview

The AZFc region of the Y chromosome contains highly repetitive
ampliconic sequences that are prone to structural variation, including
partial deletions such as the gr/gr deletion. Accurate inference of
integer copy numbers across these amplicons is required to identify and
classify these variants.

This repository provides two Python utilities to:

1.  Convert CNVkit `.cnr` log2 ratio outputs into per-sample integer
    copy-number estimates across AZFc amplicons\
2.  Annotate known AZFc CNV signatures using predefined structural
    patterns

------------------------------------------------------------------------

## Repository structure

    AZFc_CNV_annotation/
    │
    ├── CNR2CN.py
    ├── AZFc_CNV_annotation.py
    ├── README.md
    ├── requirements.txt
    └── .gitignore

------------------------------------------------------------------------

## Requirements

Python 3.8 or higher

Install dependencies:

    pip install pandas

------------------------------------------------------------------------

## Input data

Input files must be CNVkit `.cnr` files containing at minimum:

  Column   Description
  -------- -------------------------------
  gene     Amplicon or marker identifier
  log2     Log2 copy-number ratio

Each `.cnr` file corresponds to one sample.

------------------------------------------------------------------------

## Step 1: Convert `.cnr` files to integer copy numbers

Script: `CNR2CN.py`

Example usage:

    python CNR2CN.py --cnr_file_path path/to/cnr_directory

Output: per-sample integer copy numbers across AZFc amplicons.

------------------------------------------------------------------------

## Step 2: Annotate AZFc CNV signatures

Script: `AZFc_CNV_annotation.py`

Example:

``` python
import pandas as pd
from AZFc_CNV_annotation import AZFcCNV_annotation

df = pd.read_csv("cnr_processed.tsv", sep="\t")

annotated_df = AZFcCNV_annotation(df)

annotated_df.to_csv("cnr_processed_annotated.tsv", sep="\t", index=False)
```

------------------------------------------------------------------------

## Method

Copy number inference:

1.  Convert log2 ratio to linear scale:

    copy_number ∝ 2\^(log2)

2.  Scale relative to reference copy number

3.  Average across bins

4.  Round to nearest integer

------------------------------------------------------------------------

## Intended use

Designed for:

-   AZFc structural variant analysis\
-   Y chromosome CNV analysis\
-   Whole genome sequencing CNV pipelines\
-   TGCT genetic association studies

------------------------------------------------------------------------

## Author

Subin Choi\
Institute of Cancer Research, London

------------------------------------------------------------------------

## License

MIT License
