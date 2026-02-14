# AZFc CNV Annotation and CNVkit Post‑processing Pipeline

Utilities for converting CNVkit `.cnr` outputs into integer copy‑number
estimates across AZFc amplicons and annotating known AZFc CNV
signatures, including gr/gr deletions.

This pipeline was developed to support the analysis described in:

**Association of Y-chromosomal gr/gr deletions with testicular germ cell tumour: whole-genome analysis of 198,306 individuals**


------------------------------------------------------------------------

## Repository structure

    AZFcCNV_TGCT/
    │
    ├── src/
    │   ├── CNR2CN.py
    │   └── AZFc_CNV_annotation.py
    │
    ├── example/
    │   ├── HG38_Y_HG00096.cnr
    │   ├── HG38_Y_NA18970.cnr
    │   └── HG38_Y_NA18504.cnr
    │
    ├── requirements.txt
    └── README.md

------------------------------------------------------------------------

## Requirements

Python ≥ 3.8

Install dependencies:

``` bash
pip install -r requirements.txt
```

------------------------------------------------------------------------

## Input data

Input files must be CNVkit `.cnr` files containing at minimum:

  | Column  | Description                     |
  | ------- | :------------------------------:|
  | gene    |  Amplicon or marker identifier  |
  | log2    |  Log2 copy‑number ratio         |

Each `.cnr` file represents one sample.

Example files included:

    example/HG38_Y_HG00096.cnr
    example/HG38_Y_NA18970.cnr
    example/HG38_Y_NA18504.cnr

------------------------------------------------------------------------

## Step 1: Generate copy‑number summary from CNVkit output

Run:

``` bash
python src/CNR2CN.py --cnr_file_path example --output_name Example
```

This will automatically create:

    example/output/

and generate:

    example/output/Example_cnr_processed.tsv

Example output:

  |   SAMPLE_ID    | IR1 | IR5 | Blue | Teal | Green | Red | Gray | Yellow |
  | -------------- | --- | --- | ---- | ---- | ----- | --- |  ---- | ----- |
  | HG38_Y_HG00096 |  2  |  4  |   4  |  2   |   3   |  2  |   2   |   1   |

------------------------------------------------------------------------

## Step 2: Annotate AZFc CNV signatures

Run:

``` bash
python src/AZFc_CNV_annotation.py \
    --processed_cnr_file example/output/Example_cnr_processed.tsv
```

Output will be saved to:

    example/output/Example_cnr_processed_output.tsv

Additional columns added:

  | Column              | Description
  | --------------------| -------------------------------------
  | HighConfidence_CNV  | Exact CNV signature match
  | LowConfidence_CNV   | Near match (±1 amplicon difference)

------------------------------------------------------------------------

## Output directory structure

After running the pipeline:

    example/
    │
    ├── HG38_Y_HG00096.cnr
    ├── HG38_Y_NA18970.cnr
    ├── HG38_Y_NA18504.cnr
    │
    └── output/
        ├── Example_cnr_processed.tsv
        └── Example_cnr_processed_output.tsv

------------------------------------------------------------------------

## Method overview

Copy number inference:

1.  Convert log2 ratio to linear scale\
2.  Scale relative to reference copy number\
3.  Average across bins for each amplicon\
4.  Round to nearest integer

CNV annotation:

1.  Compare inferred copy numbers to known AZFc signatures\
2.  Assign high‑confidence matches\
3.  Assign low‑confidence matches for near signatures

------------------------------------------------------------------------

## Example: full pipeline

From repository root:

``` bash
pip install -r requirements.txt

python src/CNR2CN.py --cnr_file_path example --output_name Example

python src/AZFc_CNV_annotation.py \
    --processed_cnr_file example/output/Example_cnr_processed.tsv
```

------------------------------------------------------------------------

## Author

Subin Choi\
Institute of Cancer Research, London

------------------------------------------------------------------------

## License

MIT License
