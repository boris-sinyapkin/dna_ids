import ids

from kdd         import JSON_Codetable, KDD_Dataset 
from pathlib     import Path

def run(dataset_path: Path, codetable_path: Path, normal_path: Path):

    # Prepare dataset and codetable
    codetable     = JSON_Codetable(codetable_path) 
    train_dataset = KDD_Dataset(dataset_path)

    normal_act = []

    # Prepare file with normal activities (if not exist)
    if normal_path == None or normal_path.exists() == False:
        # Select normal activities from dataset
        normal_act = train_dataset.get_normal_seq(codetable)

        # Put it in FASTA fromat file
        ids.SeqIO.write(normal_act, "normal.faa", "fasta")
    else:
        # Obtain records from file
        normal_act = [seq_record for seq_record in ids.SeqIO.parse(normal_path, "fasta")]

    assert(normal_act != [])

    # By default, a global pairwise alignment is performed
    aligner = ids.Align.PairwiseAligner()

    # Match/mismatch scores of 5/-4 and gap penalties (open/extend) of 2/0.5):
    aligner.match = 5.0
    aligner.mismatch = -4.0
    aligner.open_gap_score = -2.0
    aligner.extend_gap_score = -0.5

    ideal_seq = ids.get_best_signature(aligner, normal_act)
