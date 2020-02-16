
import csv
import json
import regex

from pandas                  import DataFrame
from numpy                   import array as numpy_array
from tqdm                    import tqdm
from pathlib                 import Path
from .datasets.interfaces    import Codetable, Dataset

from Bio           import SeqIO, Align
from Bio.Seq       import Seq
from Bio.SeqRecord import SeqRecord

from os            import path

class Metrics:
    false_pos, false_negative = 0, 0
    true_pos,  true_negative  = 0, 0
    
    accuracy    = property()
    precision   = property()
    recall      = property()
    specificity = property()

    @accuracy.getter
    def accuracy(self):
        return (self.true_pos + self.true_negative) / \
            (self.false_pos + self.false_negative + self.true_pos + self.true_negative)

    @precision.getter
    def precision(self):
        return self.true_pos / \
            (self.true_pos + self.false_pos)

    @recall.getter
    def recall(self):
        return self.true_pos / \
            (self.true_pos + self.false_negative)

    @specificity.getter
    def specificity(self):
        return self.true_negative / \
            (self.true_negative + self.false_pos)

    def update(self, test_result : bool, condition : bool) -> None:
        if test_result is True:
            if condition is True:
                self.true_pos += 1
            else:
                self.true_negative += 1
        else:
            if condition is True:
                self.false_pos += 1
            else:
                self.false_negative += 1

    def show(self):
        errors = DataFrame({"True" : [self.true_pos, self.true_negative], \
                        "False" : [self.false_pos, self.false_negative]}, index=["Positive", "Negative"])

        metrics = DataFrame({"Value" : [self.accuracy, self.precision, self.recall, self.specificity]}, \
                            index=["Accuracy", "Precision", "Recall", "Specificity"])
        print(errors)
        print(metrics)


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
    _ideal_seq  = None

    codetable       = property()
    aligner         = property()
    ideal_sequence  = property()

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
    @ideal_sequence.getter
    def ideal_sequence(self):
        return self._ideal_seq

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

        self._ideal_seq = IdealSequence(ideal_seq, threshold)
        return self._ideal_seq

    def classify(self, test_dna_seq : SeqRecord) -> bool:
        """
            Align DNA sequence with ideal and do prediction
        """
        return self.ideal_sequence.test(self.aligner, test_dna_seq)

    def test(self, test_dataset : Dataset) -> Metrics:
        
        if self.ideal_sequence is None:
            raise Exception("Ideal sequence is None")

        metrics = Metrics()
        
        for i in tqdm(range(0, len(test_dataset)), desc="Testing process"):

            # Obtain DatasetRecord instance
            ds_rec = test_dataset[i]
            # Encode it in DNA
            test_dna_seq = ds_rec.encode_into_DNA(self.codetable)
            # Test it with IdealSequence
            RESULT = self.classify(test_dna_seq)
            LABEL  = test_dna_seq.name

            metrics.update(RESULT, LABEL == "attack")

        return metrics

    def analyze(self, train_ds : Dataset, test_ds : Dataset, sizes=[10, 8, 6, 4, 2, 1]):
        """
            Test IDS metrics on different sample size of train dataset.
            @sizes - the size of parts of test dataset is using.

            For example:
                sizes=[10, 5] means that (test_ds_size/10) and (test_ds_size/5) will be taken. 
        """
        TRAIN_DS_SIZE = len(train_ds)
        METRICS_DICT  = {}     
        for s in sizes:
            SAMPLE_SIZE = int(TRAIN_DS_SIZE / s)
            # Get sample of train dataset
            current_train_ds = train_ds.random_sample(SAMPLE_SIZE)
            # Obtain ideal sequence
            self.train(current_train_ds)
            # Get metrics
            metrics = self.test(test_ds) 

            METRICS_DICT.update({
                "size"    : s,
                "metrics" : metrics
            })
        
        return METRICS_DICT

           

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