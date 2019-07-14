
import csv
import ids

from pathlib import Path

from Bio     import Align
from Bio     import SeqIO

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
        SeqIO.write(normal_act, "normal.faa", "fasta")
    else:
        # Obtain records from file
        normal_act = [seq_record for seq_record in SeqIO.parse(normal_path, "fasta")]

    raise "List of normal activities is empty." if normal_act == [] else {}

    # By default, a global pairwise alignment is performed
    aligner = Align.PairwiseAligner()

    # Match/mismatch scores of 5/-4 and gap penalties (open/extend) of 2/0.5):
    aligner.match            =  5.0
    aligner.mismatch         = -4.0
    aligner.open_gap_score   = -2.0
    aligner.extend_gap_score = -0.5

        


        

