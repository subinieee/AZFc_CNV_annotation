import pandas as pd


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


d=process_cnr_files(cnr_files)
def prep_output_column(d):

    d['HighConfidence_CNV'] = str(0)
    d['LowConfidence_CNV'] = str(0)
    return d
def calculate_cn_change(d):
    for key in copy_numbers.keys():
        d[key] = d[key] - copy_numbers[key]
    return d
def known_pheno_annotation(d):

    for i in AZFcCNV_ANNOTATION.index:
        d.loc[
            ((d[AZFcCNV_ANNOTATION.loc[i].keys()]) == AZFcCNV_ANNOTATION.loc[i]).all(1), 'HighConfidence_CNV'] = i
        d.loc[
            (abs(d[AZFc] - AZFcCNV_ANNOTATION.loc[i][AZFc]).sum(axis=1) == 1)
            & ((d[ANNOTATION_DICT[i].keys()] == AZFcCNV_ANNOTATION.loc[i][ANNOTATION_DICT[i].keys()]).all(1))
            | ((d[AZFcCNV_ANNOTATION.loc[i].keys()]) == AZFcCNV_ANNOTATION.loc[i]).all(1)
            , 'LowConfidence_CNV'] = i
    return d

def AZFcCNV_annotation(d: pd.DataFrame) -> pd.DataFrame:
    """Run full AZFc CNV annotation pipeline on input dataframe."""
    return (
        d.pipe(prep_output_column)
          .pipe(calculate_cn_change)
          .pipe(known_pheno_annotation)
    )