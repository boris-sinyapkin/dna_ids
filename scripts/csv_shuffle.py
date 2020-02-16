import argparse
import pandas as pd

from pathlib import Path

parser = argparse.ArgumentParser(description="Shuffle rows in csv file")
parser.add_argument('--input',  '-i', required=True, type=Path)
parser.add_argument('--output', '-o', required=True, type=Path)

args = parser.parse_args()

INPUT  = args.input
OUTPUT = args.output

df = pd.read_csv(INPUT)
shuffled_df = df.sample(frac=1)
shuffled_df.to_csv(OUTPUT, index=False)