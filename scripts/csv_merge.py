import argparse
import pandas as pd

from pathlib import Path

def shuffle_df(csv_df : pd.DataFrame) -> pd.DataFrame:
    return csv_df.sample(frac=1)

def merge_csv(fst : Path, sec : Path) -> pd.DataFrame:
    fst_df = pd.read_csv(fst, low_memory=False)
    sec_df = pd.read_csv(sec, low_memory=False)

    fst_df = fst_df.fillna(0)
    sec_df = sec_df.fillna(0)
    return fst_df.append(sec_df)
    
parser = argparse.ArgumentParser(description="Shuffle rows in csv file")
parser.add_argument('--fst',    '-f', required=True, type=Path, help="Path to 1st csv file")
parser.add_argument('--sec',    '-s', required=True, type=Path, help="Path to 2nd csv file")
parser.add_argument('--output', '-o', required=True, type=Path)

args = parser.parse_args()

FIRST  = args.fst
SECOND = args.sec
OUTPUT = args.output

shuffle_df(merge_csv(FIRST, SECOND)).to_csv(OUTPUT, index=False)