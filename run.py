
import csv
import ids

from pathlib import Path

def run(dataset_path : Path, codetable_path : Path, normal_path : Path):

    # Prepare dataset
    dataset = ids.kdd_dataset(dataset_path)

    # Set codetable
    dataset.set_codetable(codetable_path)

    normal_act = []

    # Prepare file with normal activities (if not exist)
    if normal_path.exists() == False:
        # Select normal activities from dataset
        normal_act = dataset.findall("normal",
                                     lambda val, seq_rec: val == seq_rec.name)
        # Put it in FASTA fromat file
        ids.SeqIO.write(normal_act, "normal.faa", "fasta")
    else:
        # Obtain records from file
        normal_act = [seq_record for seq_record in ids.SeqIO.parse(normal_path, "fasta")]

    assert(normal_act != [])

    # By default, a global pairwise alignment is performed
    aligner = ids.Align.PairwiseAligner()

    # Match/mismatch scores of 5/-4 and gap penalties (open/extend) of 2/0.5):
    aligner.match            =  5.0
    aligner.mismatch         = -4.0
    aligner.open_gap_score   = -2.0
    aligner.extend_gap_score = -0.5

    ideal_seq = ids.get_normal_sequence(aligner, normal_act)
        


        

