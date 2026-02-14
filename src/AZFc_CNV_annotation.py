import pandas as pd
from pathlib import Path
import argparse

# -------------------------------
# Reference definitions
# -------------------------------

AZFc = ['IR1', 'IR5', 'Blue', 'Teal', 'Green', 'Red', 'Gray', 'Yellow']
cn_2 = ['IR3', 'IR1', 'P8', 'P7', 'P6', 'P5', 'P4', 'IR2', 'Teal', 'Gray', 'Yellow']
cn_3 = ['Green']
cn_4 = ['IR5', 'Blue', 'Red']

#Create a dictionary for AZFc Amplicon : referenece copy number
copy_numbers = {}
for amp in cn_2:
    copy_numbers[amp] = 2
for amp in cn_3:
    copy_numbers[amp] = 3
for amp in cn_4:
    copy_numbers[amp] = 4

df = pd.DataFrame(columns=AZFc, index=['v' + str(x) for x in range(1, 14)], )
ANNOTATION = pd.DataFrame(data=df)

# KNOWN NAHR PHENOTYPE
ANNOTATION_DICT = dict(
    v1={'IR5': -1, 'Blue': -1, 'Green': -1, 'Red': -2, 'Gray': -1, 'Yellow': -1},
    v2={'IR5': 1, 'Blue': 1, 'Green': 1, 'Red': 2, 'Gray': 1, 'Yellow': 1},
    v3={'IR1': -1, 'IR5': -1, 'Blue': -1, 'Green': -2, 'Red': -2, 'Yellow': -1},
    v4={'IR1': 1, 'IR5': 1, 'Blue': 1, 'Green': 2, 'Red': 2, 'Yellow': 1},
    v5={'IR1': 1, 'Green': 1, 'Gray': -1},
    v6={'IR1': -1, 'Green': -1, 'Gray': 1},
    v7={'IR1': 1, 'IR5': 2, 'Blue': 2, 'Green': 3, 'Red': 4, 'Gray': 1, 'Yellow': 2},
    v8={'IR1': 1, 'Blue': 2, 'Teal': 2, 'Green': 1, 'Red': 2, 'Gray': 1},
    v9={'IR1': -1, 'Blue': -2, 'Teal': -2, 'Green': -1, 'Red': -2, 'Gray': -1},
    v10={'IR5': 2, 'Blue': 2, 'Green': 2, 'Red': 4, 'Gray': 2, 'Yellow': 2},
    v11={'IR1': 2, 'IR5': 2, 'Blue': 2, 'Green': 4, 'Red': 4, 'Yellow': 2},
    v12={'IR5': 4, 'Blue': 4, 'Green': 4, 'Red': 8, 'Gray': 4, 'Yellow': 4},
    v13={'IR1': 3, 'IR5': 3, 'Blue': 3, 'Green': 6, 'Red': 6, 'Yellow': 3},
    v14={'IR1': 2, 'IR5': 1, 'Blue': 1, 'Green': 3, 'Red': 2, 'Gray': -1, 'Yellow': 1},
    whole={'IR1': -1, 'IR5': -2, 'Blue': -2, 'Green': -3, 'Red': -4, 'Gray': -1, 'Yellow': -2},
)
ANNOTATION = pd.DataFrame(ANNOTATION_DICT).fillna(0).astype(int).T
ANNOTATION = ANNOTATION[AZFc]

def prep_output_column(d):

    d['HighConfidence_CNV'] = str(0)
    d['LowConfidence_CNV'] = str(0)
    return d
def calculate_cn_change(d):
    for key in copy_numbers.keys():
        d[key] = d[key] - copy_numbers[key]
    return d
def known_pheno_annotation(d):

    for i in ANNOTATION.index:
        d.loc[
            ((d[ANNOTATION.loc[i].keys()]) == ANNOTATION.loc[i]).all(axis=1), 'HighConfidence_CNV'] = i
        d.loc[
            (abs(d[AZFc] - ANNOTATION.loc[i][AZFc]).sum(axis=1) == 1)
            & ((d[ANNOTATION_DICT[i].keys()] == ANNOTATION.loc[i][ANNOTATION_DICT[i].keys()]).all(axis=1))
            | ((d[ANNOTATION.loc[i].keys()]) == ANNOTATION.loc[i]).all(axis=1)
            , 'LowConfidence_CNV'] = i
    return d

def AZFcCNV_annotation(d: pd.DataFrame) -> pd.DataFrame:
    """Run full AZFc CNV annotation pipeline on input dataframe."""
    return (
        d.pipe(prep_output_column)
          .pipe(calculate_cn_change)
          .pipe(known_pheno_annotation)
    )
if __name__ == "__main__":

   parser = argparse.ArgumentParser(description="Annotate AZFc CNVs")
   parser.add_argument("--processed_cnr_file", help="Path to the processed cnr file (.tsv)",
                        required=True)
   args = parser.parse_args()
   cnr_path = Path(args.processed_cnr_file)
   d=pd.read_csv(cnr_path, sep="\t")
   annotated_df = AZFcCNV_annotation(d)
   annotated_df.to_csv(cnr_path.parent / f"{cnr_path.stem}_output.tsv", sep="\t", index=False)
   print("Output saved to output.tsv")