
import csv
import json
import regex


from numpy                   import array as numpy_array
from tqdm                    import tqdm
from pathlib                 import Path
from .datasets.interfaces    import Codetable, Dataset

from Bio           import SeqIO, Align
from Bio.Seq       import Seq
from Bio.SeqRecord import SeqRecord

from os            import path

class AlignInfo:

    _info = {}

    sum       = property()
    seqrec    = property()
    threshold = property()

    def __init__(self, align_info : dict):
        self._info = align_info

    def to_csv(self) -> list:
        csv_line = [self._info["seq"].id]
        csv_line += [item["score"] for item in self._info["score_list"]]
        csv_line += [self._info["threshold"], self._info["sum"]]

        return csv_line

    @sum.getter
    def sum(self):
        return self._info["sum"]
    @seqrec.getter
    def seqrec(self):
        return self._info["seq"]
    @threshold.getter
    def threshold(self):
        return self._info["threshold"]

class IdealSequence(SeqRecord):

    threshold = property()

    def __init__(self, seqrec : SeqRecord, threshold):
        super().__init__(seqrec.seq, 
                            id = str(seqrec.id), 
                            name = seqrec.name, 
                            description = seqrec.description)

        self._threshold = threshold

    @threshold.getter
    def threshold(self):
        return self._threshold

    # True - attack, False - normal activity
    def test(self, aligner : Align.PairwiseAlignment, other_seqrec : SeqRecord) -> bool:
        return True \
                if aligner.score(other_seqrec.seq, self.seq) <= self._threshold \
                    else False

    def dump(self, dest_dir : Path):
        """
            Puts ideal sequence in FASTA format in specified directory
        """
        dest_dir.mkdir()

        SeqIO.write(self, (dest_dir / "sequence.faa"), "fasta")

        with (dest_dir / Path("info.json")).open("w") as thold_file:
            json.dump({ "threshold" : self._threshold }, thold_file)
    
    @staticmethod
    def load(src_dir : Path):
        """
            Reads ideal sequence in FASTA format from specified directory

            src_dir/sequence.faa - SeqRecord in FASTA format
            src_dir/info.json    - Threshold
        """
        if path.isdir(src_dir):

            ideal_seqrec = SeqIO.read((src_dir / "sequence.faa"), "fasta")

            with (src_dir / Path("info.json")).open("r") as thold_file:
                threshold = json.load(thold_file)["threshold"]

            return IdealSequence(ideal_seqrec, threshold)

        else:
            raise Exception(f"Directory with ideal sequence was not found: {src_dir}")

            

class IDS:
    _codetable  = None
    _aligner    = None

    codetable   = property()
    aligner     = property()

    def __init__(self, codetable : Codetable, aligner : Align.PairwiseAlignment):
        self.codetable = codetable
        self.aligner   = aligner

    @codetable.setter
    def codetable(self, codetable : Codetable):
        self._codetable = codetable
    @aligner.setter
    def aligner(self, aligner : Align.PairwiseAlignment):
        self._aligner = aligner
    
    @codetable.getter
    def codetable(self):
        return self._codetable
    @aligner.getter
    def aligner(self):
        return self._aligner

    # This function search ideal sequence in train dataset
    def train(self, train_dataset : Dataset) -> IdealSequence:

        # Ideal sequence attributes
        ideal_seq = None 
        threshold = None

        train_ds_dna = train_dataset.as_DNA_records(self.codetable)

        max_sum = 0
        
        for ds_record in tqdm(train_ds_dna, desc="Training process"):

            align_info = calculate_vector_alignment(self.aligner, ds_record, train_ds_dna)

            if align_info.sum > max_sum:
                ideal_seq   = align_info.seqrec
                threshold   = align_info.threshold
                max_sum     = align_info.sum

        return IdealSequence(ideal_seq, threshold)   

    def test(self, test_dataset : Dataset, ideal_seq : IdealSequence) -> None:

        # Metrics
        false_pos, false_negative = 0, 0
        true_pos,  true_negative  = 0, 0
        
        for i in tqdm(range(0, len(test_dataset)), desc="Testing process"):

            # Obtain DatasetRecord instance
            ds_rec = test_dataset[i]
            # Encode it in DNA
            test_dna_seq = ds_rec.encode_into_DNA(self.codetable)
            # Test it with IdealSequence
            RESULT = ideal_seq.test(self.aligner, test_dna_seq)
            LABEL  = test_dna_seq.name

            # Calculate metrics
            if RESULT is True:
                if LABEL == "normal":
                    false_pos +=1
                else:
                    false_negative +=1

            elif RESULT is False:
                if LABEL == "normal":
                    true_pos +=1
                else:
                    true_negative +=1

        print(f"False Pos: {false_pos}")
        print(f"False Neg: {false_negative}")
        print(f"True  Pos: {true_pos}")
        print(f"True  Neg: {true_negative}")
           

# Align 'seq_rec' with each sequence in 'seq_vector'
def calculate_vector_alignment(aligner : Align.PairwiseAligner, seq_rec : SeqRecord, seq_vector : numpy_array) -> AlignInfo:
    
    score_list = []

    # Sum of align scores
    score_sum  = 0

    for list_item in seq_vector:

            current_score = aligner.score(seq_rec.seq, list_item.seq)  
            score_list.append({
                "id"    : list_item.id,
                "score" : current_score
            })
            score_sum += current_score

    threshold  = score_sum / len(score_list)

    return AlignInfo({
        "seq"         : seq_rec,
        "score_list"  : score_list,
        "threshold"   : threshold,
        "sum"         : score_sum
    })