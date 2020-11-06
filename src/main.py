from .datasets.csv_ds      import CSV_Dataset
from .datasets.interfaces  import JSON_Codetable
from .utils                import create_shuffled_test_df
from .ids                  import IDS, Align
from pathlib               import Path

def run( train_ds_path: Path, test_ds_path: Path, codetable_path : Path):

    CODETABLE = JSON_Codetable(codetable_path)          if codetable_path   else None
    TRAIN_DS  = CSV_Dataset.from_file(train_ds_path)    if train_ds_path    else None
    TEST_DS   = CSV_Dataset.from_file(test_ds_path)     if test_ds_path     else None

    # By default, a Smith-Waterman alignment if performed
    ALIGNER = IDS.Aligner()

    # Create IDS instance with Codetable & Aligner
    ids = IDS(CODETABLE, ALIGNER)
    
    mixed_test_ds = CSV_Dataset(create_shuffled_test_df(TEST_DS, TRAIN_DS))
    
    ids.analyze(TRAIN_DS, mixed_test_ds, sizes=[1]).to_excel("Metrics.xlsx")