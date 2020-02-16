
from .datasets.ieee        import IEEE_Dataset
from .datasets.interfaces  import JSON_Codetable

from .ids                  import IDS, Align, IdealSequence
from pathlib               import Path

def run( train_ds_path: Path, test_ds_path: Path, codetable_path : Path, normal_act_path: Path):

    CODETABLE = JSON_Codetable(codetable_path) if codetable_path   else None
    TRAIN_DS  = IEEE_Dataset(train_ds_path)    if train_ds_path    else None
    TEST_DS   = IEEE_Dataset(test_ds_path)     if test_ds_path     else None

    # By default, a global pairwise alignment is performed
    ALIGNER = Align.PairwiseAligner()

    # Match/mismatch scores of 5/-4 and gap penalties (open/extend) of 2/0.5):
    ALIGNER.match            =  5.0
    ALIGNER.mismatch         = -4.0
    ALIGNER.open_gap_score   = -2.0
    ALIGNER.extend_gap_score = -0.5

    # Create IDS instance with Codetable & Aligner
    ids = IDS(CODETABLE, ALIGNER)

    # Get Ideal Sequence
    # ideal_sequence = ids.train(TRAIN_DS)
    # ideal_sequence.dump(Path("ideal1"))

    ideal_sequence = IdealSequence.load(Path("ideal1"))

    # Get IDS metrics
    ids.test(TEST_DS, ideal_sequence)