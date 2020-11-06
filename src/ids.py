
import json
import multiprocessing as mp  
import pandas as pd

from multiprocessing.dummy   import Pool as ThreadPool
from math                    import floor
from numpy                   import array as numpy_array
from tqdm                    import tqdm
from pathlib                 import Path
from .datasets.interfaces    import Codetable, Dataset
from Bio                     import SeqIO, Align
from Bio.Seq                 import Seq
from Bio.SeqRecord           import SeqRecord
from os                      import path

class Metrics:
    false_pos, false_negative = 0, 0
    true_pos,  true_negative  = 0, 0
    
    accuracy    = property()
    precision   = property()
    recall      = property()
    specificity = property()

    @accuracy.getter
    def accuracy(self):
        try:
            value = (self.true_pos + self.true_negative) / \
                (self.false_pos + self.false_negative + self.true_pos + self.true_negative)
        except ZeroDivisionError:
            value = None
        return value

    @precision.getter
    def precision(self):
        try:
            value = self.true_pos / (self.true_pos + self.false_pos)
        except ZeroDivisionError:
            value = None
        return value

    @recall.getter
    def recall(self):
        try:
            value = self.true_pos / (self.true_pos + self.false_negative)
        except ZeroDivisionError:
            value = None
        return value

    @specificity.getter
    def specificity(self):
        try:
            value = self.true_negative / (self.true_negative + self.false_pos)
        except ZeroDivisionError:
            value = None
        return value
    
    def __add__(self, other):
        result = self
        result.true_pos         += other.true_pos
        result.true_negative    += other.true_negative
        result.false_pos        += other.false_pos
        result.false_negative   += other.false_negative
        return result
        
    def update(self, test_result : bool, condition : bool) -> None:
        if test_result is True:
            if condition is True:
                self.true_pos += 1
            else:
                self.false_pos += 1
        else:
            if condition is True:
                self.false_negative += 1
            else:
                self.true_negative += 1
                
    def as_dataframe(self, train_size, test_size) -> pd.DataFrame:
        return pd.DataFrame({
                    "True Positive"     : [self.true_pos],
                    "True Negative"     : [self.true_negative],     
                    "False Positive"    : [self.false_pos],
                    "False Negative"    : [self.false_negative],
                    "Accuracy"          : [self.accuracy],
                    "Precision"         : [self.precision],
                    "Recall"            : [self.recall],
                    "Specificity"       : [self.specificity],
                    "Test subset size"  : [test_size],
                    "Train subset size" : [train_size]})

    def show(self):
        errors = pd.DataFrame({"True" : [self.true_pos, self.true_negative], \
                        "False" : [self.false_pos, self.false_negative]}, index=["Positive", "Negative"])

        metrics = pd.DataFrame({"Value" : [self.accuracy, self.precision, self.recall, self.specificity]}, \
                            index=["Accuracy", "Precision", "Recall", "Specificity"])
        print(errors)
        print(metrics)

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
                if aligner.score(self.seq, other_seqrec.seq) < self._threshold \
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
    class Aligner(Align.PairwiseAligner):
        def __init__(self, algo='Smith-Waterman') -> None:
            super().__init__()
            if algo == 'Smith-Waterman':
                self.match            = 1.0
                self.mismatch         = 0
                self.open_gap_score   = 0
                self.extend_gap_score = 0
                self.mode             = 'local'
            elif algo == 'Gotoh':
                self.match            =  5.0
                self.mismatch         = -4.0
                self.open_gap_score   = -2.0
                self.extend_gap_score = -0.5
                self.mode             = 'global'
            print(f"Alignment algo: {self.algorithm}")
    
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

    def get_multiple_align_score(self, seq : Seq, dna_sequences : numpy_array, tasks : list, pool : ThreadPool) -> float:
        """
            Align seq with each record in dna_sequences list and return sum of alignment scores.
                @dna_sequences - list of encoded DNA sequences
                @tasks - list of (start_index, finish_index) tuples for each thread in @pool
        """                
        def _calc_multiple_alignment_score(wrapped_data : tuple) -> int:
            """return value is a sum of all alignment scores in input rangerange."""   
            (start, finish) = wrapped_data       
            score_sum = 0.
            for dna_record in tqdm(dna_sequences[start : finish + 1], total=(finish + 1 - start), desc="Training process"):
                score_sum += self.aligner.score(seq, dna_record.seq)
            return score_sum
        return sum(s for s in pool.map(_calc_multiple_alignment_score, tasks))

    # This function search ideal sequence in train dataset
    def train(self, train_dataset : Dataset, proc_num = mp.cpu_count()) -> IdealSequence:
        
        mean_row = train_dataset.get_median()
                
        mean_row_dna = mean_row.encode_into_DNA(self.codetable, id='')
        train_ds_dna = train_dataset.as_DNA_records(self.codetable)
        
        # Prepare environment for parallel execution
        SIZE, THREADS = len(train_ds_dna), proc_num
        TASKS = [ (start, finish) for start, finish in self._intervals(THREADS, SIZE) ]

        with ThreadPool(THREADS) as pool:
            score_sum = self.get_multiple_align_score(mean_row_dna.seq, train_ds_dna, TASKS, pool) 

        self._ideal_seq = IdealSequence(mean_row_dna, score_sum / SIZE)
        return self._ideal_seq

    def classify(self, test_dna_seq : SeqRecord) -> bool:
        """Align DNA sequence with ideal and do prediction"""
        return self.ideal_sequence.test(self.aligner, test_dna_seq)
    
    @staticmethod
    def _intervals(parts, duration) -> list:
        """ Private function. Split duration into ranges."""
        part_duration = duration / parts
        return [(floor(i * part_duration), floor((i + 1) * part_duration)) for i in range(parts)]
    
    def test(self, test_dataset : Dataset, proc_num = mp.cpu_count()) -> Metrics:
        
        if self.ideal_sequence is None:
            raise Exception("Ideal sequence is None")
        
        def test_worker(wrapped_data : tuple) -> Metrics:
            (start, finish), metrics = wrapped_data, Metrics()
            for i in tqdm(range(start, finish), desc="Testing process"):  
                # Obtain DatasetRecord instance
                ds_rec = test_dataset.raw_index(i)
                # Encode it in DNA
                test_dna_seq = ds_rec.encode_into_DNA(self.codetable, id=str(i))            
                # Test it with IdealSequence
                metrics.update(self.classify(test_dna_seq), test_dna_seq.name == "attack")
                
            return metrics
                        
        SIZE, PROCS = len(test_dataset), proc_num 
        TASKS   = [ (start, finish) for start, finish in self._intervals(PROCS, SIZE)]
        METRICS = Metrics()
        
        # Make the Pool of workers
        with ThreadPool(PROCS) as pool:
            for met in pool.map(test_worker, TASKS):   
                METRICS = METRICS + met

        return METRICS

    def analyze(self, train_ds : Dataset, test_ds : Dataset, sizes=[10, 8, 6, 4, 2, 1]) -> pd.DataFrame:
        """
            Test IDS metrics on different sample size of train dataset.
            @sizes - the size of parts of test dataset is using.

            For example:
                sizes=[10, 5] means that (test_ds_size/10) and (test_ds_size/5) will be taken. 
        """
        TRAIN_DS_SIZE = len(train_ds)
        METRICS       = pd.DataFrame(None)    

        # Test IDS instance
        ids = IDS(self.codetable, self.aligner)

        for s in sizes:
            SAMPLE_SIZE = int(TRAIN_DS_SIZE / s)
            print(f"Train Dataset size : {SAMPLE_SIZE}")  
            # Get sample of train dataset
            current_train_ds = train_ds.random_sample(SAMPLE_SIZE)            
            # Obtain ideal sequence
            ids.train(current_train_ds)
            # Get metrics
            metrics = ids.test(test_ds).as_dataframe(SAMPLE_SIZE, len(test_ds))
            
            METRICS = METRICS.append(metrics)
                    
        return METRICS