"""
    This script is using to train & test this IDS on different datasets.
    Note: For now, only *.csv format of datasets are supported.
    
    Parameters:
        @dir: directory in which train & test dataset files will be searched
    
    Artifacts:
        This script generate metrics.csv file.

    Warning: This script must be located in scripts folder to correct import of IDS modules.
"""
import argparse

from pathlib import Path
from os      import walk, path
from sys     import path as syspath
from Bio     import Align

#------------------
# Argument parsing
#------------------
parser = argparse.ArgumentParser(description="This script is using to train & test this IDS on different datasets.")
parser.add_argument("--dir",     "-d", type=Path, required=True)
parser.add_argument("--output" , "-o", type=Path, required=True) 
args = parser.parse_args()

DIR         = args.dir
OUTPUT      = args.output
SCRIPT_DIR  = Path(path.dirname(path.abspath(__file__)))

# Hardcoded dataset filenames
TEST_DS_NAME  = Path("test.csv")
TRAIN_DS_NAME = Path("train.csv")
CODETABLE     = Path("codetable.json")

# For this fractions of train dataset sizes metrics will be obtained
TRAIN_SIZES = [100, 90, 80, 70]

# Import IDS modules
syspath.append(path.join(SCRIPT_DIR, ".."))
from src.datasets.interfaces    import JSON_Codetable
from src.datasets.csv_ds        import CSV_Dataset  
from src.ids                    import Metrics, IDS

def run_test(TEST_DIR : Path) -> Metrics:

    TEST_DS    = CSV_Dataset.from_file(TEST_DIR / TEST_DS_NAME)
    TRAIN_DS   = CSV_Dataset.from_file(TEST_DIR / TRAIN_DS_NAME) 
    JSON_CODES = JSON_Codetable(CODETABLE)
    
    # By default, a global pairwise alignment is performed
    ALIGNER = Align.PairwiseAligner()

    # Match/mismatch scores of 5/-4 and gap penalties (open/extend) of 2/0.5):
    ALIGNER.match            =  5.0
    ALIGNER.mismatch         = -4.0
    ALIGNER.open_gap_score   = -2.0
    ALIGNER.extend_gap_score = -0.5

    # Create IDS instance with Codetable & Aligner and calculate testing & training metrics
    return IDS(JSON_CODES, ALIGNER).analyze(TRAIN_DS, TEST_DS, sizes=TRAIN_SIZES)
    

#-----------------
# Entry point
#-----------------
def main():

    TEST_DIRS = []
    
    # Find of test directories
    for dirpath, dirnames, filenames in walk(DIR):
        if not dirnames:
            if TEST_DS_NAME and TRAIN_DS_NAME and CODETABLE in filenames:
                TEST_DIRS.append(Path(dirpath))
    
    if not TEST_DIRS:
        raise Exception("No test directories was found")
    else:    
        # Execute IDS testing for each directory
        for test_dir in TEST_DIRS:
            run_test(test_dir).to_csv(test_dir / Path("metrics.csv"),  index=False)
        
if __name__ == "__main__":
    main()