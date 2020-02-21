"""
    This script is using to train & test this IDS on different datasets.
    Note: For now, only *.csv format of datasets are supported.
    
    Parameters:
        @dir: directory in which train.csv and test.csv files are located.
    
    Artifacts:
        This script generate metrics.csv file.
"""
import argparse
from pathlib import Path

#------------------
# Argument parsing
#------------------
parser = argparse.ArgumentParser(description="This script is using to train & test this IDS on different datasets.")
parser.add_argument("--dir", "-d", type=Path, required=True)
parser.add_argument("--output" , "-o", type=Path, required=True) 
args = parser.parse()

DIR     = args.dir
OUTPUT  = args.output

#-----------------
# Entry point
#-----------------
def main():
    pass
    

if __name__ == "__main__":
    main()