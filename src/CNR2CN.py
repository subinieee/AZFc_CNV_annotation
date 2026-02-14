import pandas as pd
from pathlib import Path
import argparse
import math

# -------------------------------
# Locate .cnr files, calculate mean read depth grouped by STS and AZFc amplicons and round up to nearest integer
# -------------------------------

def main():

    parser = argparse.ArgumentParser(description="Process CNR files")
    parser.add_argument("--cnr_file_path", help="Path to the folder directory containing .cnr files",
                        required=True
    )
    parser.add_argument("--output_name", default='Example', help="Output file name for the processed cnr file",
                        required=False
    )


    args = parser.parse_args()
    cnr_path = Path(args.cnr_file_path)
    output_dir = cnr_path / "output"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_name = args.output_name
    cnr_files = list(cnr_path.glob('*.cnr'))

    
    #Create a dictionary for AZFc Amplicon : referenece copy number
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


    # Amplicon regions grouping
    amplicon_regions = {'IR3': ['IR3-1', 'IR3-2', 'IR3-3'],
                        'IR1': ['IR1-1', 'IR1-2', 'IR1-3'], 'P8': ['P8'],
                        'P7': ['P7'], 'P6': ['P6'], 'P5': ['P5'],
                        'IR5': ['IR5'], 'P4': ['P4'],
                        'IR2': ['IR2'], 'Blue': ['Blue'],
                        'Teal': ['Teal-1', 'Teal-2'], 'Green': ['Green'],
                        'Red': ['Red-1', 'Red-2'],
                      #  'DAZ': ['DAZ'],
                        'Gray': ['Gray'],
                        'Yellow': ['Yellow-1', 'P1.1/2', 'P1.1/2-spacer',
                                   'P1-chr15-1', 'P1.3/4', 'P1.3/4-spacer', 'P1-ch15-2']}

    # Create a DataFrame with processed cnr data
    cnr_processed = pd.DataFrame([f.stem for f in cnr_files])
    cnr_processed.columns = ['SAMPLE_ID']

    # Calculate log mean read depth of each STS and AZFc amplicon regions
    for i in cnr_files:
        print(i)
        temp = pd.read_table(i)
        cnr_processed.loc[cnr_files.index(i), 'sY1191'] = math.log2((2 ** (temp.loc[temp.gene.str.contains('sY1191'), 'log2'])).dropna().mean())
        cnr_processed.loc[cnr_files.index(i), 'sY1291'] = math.log2((2 ** (temp.loc[temp.gene.str.contains('sY1291'), 'log2'])).dropna().mean())
        cnr_processed.loc[cnr_files.index(i), 'gr/gr'] = math.log2((2 ** (temp.loc[temp.gene.str.contains('gr/gr'), 'log2'])).dropna().mean())
        for k in amplicon_regions.keys():
            cnr_processed.loc[i, k] = round(copy_numbers[k] * (
                            2 ** (temp.loc[temp.gene.str.contains('|'.join(amplicon_regions[k])), 'log2']).mean())
                            )
    cnr_processed.to_csv(cnr_path / f"output/{output_name}_cnr_processed.tsv", sep="\t", index=False)
    



if __name__ == "__main__":
    main()