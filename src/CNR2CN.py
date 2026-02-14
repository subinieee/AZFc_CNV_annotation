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

    args = parser.parse_args()
    cnr_path = Path(args.cnr_path)
    cnr_files = list(cnr_path.glob('*cnr'))

    # Amplicon regions grouping
    amplicon_regions = {'IR3': ['IR3-1', 'IR3-2', 'IR3-3'],
                        'IR1': ['IR1-1', 'IR1-2', 'IR1-3'], 'P8': ['P8'],
                        'P7': ['P7'], 'P6': ['P6'], 'P5': ['P5'],
                        'IR5': ['IR5'], 'P4': ['P4'],
                        'IR2': ['IR2'], 'Blue': ['Blue'],
                        'Teal': ['Teal-1', 'Teal-2'], 'Green': ['Green'],
                        'Red': ['Red-1', 'Red-2'],
                        'DAZ': ['DAZ'],
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
            for j in range(0, len(amplicon_regions[k])):
                cnr_processed.loc[i, k] = round(copy_numbers[k] * (
                            2 ** (temp.loc[temp.gene.str.contains('|'.join(amplicon_regions[k])), 'log2']).mean()))
    return cnr_processed


if __name__ == "__main__":
    main()