
import csv
import json
import regex

from pathlib import Path

from Bio           import Align
from Bio           import SeqIO
from Bio.Seq       import Seq
from Bio.SeqRecord import SeqRecord

from sys           import stderr


class kdd_dataset:

    _path = ""
    _dataset = []
    _codetable = {}
    
    def __init__(self, path : str):
        self._path = path

        print("Loading the dataset from : %s" % (path))
        with open(path, 'r') as raw_dataset:
            self._dataset = [row for row in csv.reader(raw_dataset)]

    def __getitem__(self, index):
        return self._dataset[index]

    #load codetable from json file
    def set_codetable(self, path : str):

        print("Loading the codetable from : %s" % (path))
        with open(path, 'r') as json_codetable:
            self._codetable = json.load(json_codetable)

    # Encoding KDD record into DNA sequence
    def encode_kdd_record(self, index : int) -> SeqRecord:
        if not self._codetable:
            raise Exception("Codetable is empty")
        else:
            raw_record = self._dataset[index]
            
            payload = ""
            label   = ""
            extra   = []
            
            for idx, item in enumerate(raw_record):

                # in kdd record only 41 attribute
                if idx > 40:
                    break

                if regex.match(r'\d{1,}', item) or regex.match(r'\d{1,}\.\d{1,}', item):
                    # encoding numeric value
                    for literal in list(item):
                        if regex.match(r'\d{1,}', literal): 
                            # digit
                            payload += self._codetable["digits"][literal]
                        else:
                            # dot
                            payload += self._codetable["dot"]

                else:
                    try:
                        # hardcoded value
                        ident = self._codetable["prefixes"][f"{idx}"][0]
                        prefix = self._codetable["prefixes"][f"{idx}"][1]

                        payload += prefix + self._codetable["hardcoded"][ident][item]
                    except KeyError:
                        print("Key Error: \n\
                               \r- idx    = %d \n\
                               \r- ident  = \'%s\' \n\
                               \r- item   = \'%s\'" % (idx, ident, item), file=stderr)

            label = raw_record[41]

            for i in range(42, len(raw_record)):
                extra.append(raw_record[i])
                
        return SeqRecord \
        (
            Seq(payload), 
            id=str(index),
            name=label, 
            features=extra, 
            description=f"({index}) KDD record encoded into DNA. [{label}]"
        )

    # Search all sequences, which applies function in parameter list.
    def findall(self, value : str, func) -> list:
        retlist = []
        dataset_len = len(self._dataset)

        print("Processing sequences: ")
        print("                      ]\r[", end='', flush=True)

        for index in range(0, dataset_len):
            
            seq_rec = self.encode_kdd_record(index)

            if func(value, seq_rec):
                retlist.append(seq_rec)

            print("#", end='', flush=True) if index % int(dataset_len/20) == 0 else {}

        print("")
        return retlist

class AlignInfo:

    _info = {}

    def __init__(self, align_info : dict):
        self._info = align_info

    def to_csv(self) -> list:
        csv_line = [self._info["seq"].id]
        csv_line += [item["score"] for item in self._info["score_list"]]
        csv_line += [self._info["treshold"], self._info["sum"]]

        return csv_line

    def get_sum(self):
        return self._info["sum"]
    def get_seqrec(self):
        return self._info["seq"]
    def get_treshold(self):
        return self._info["treshold"]

        

# This function search normal sequence in list of normal activities             
def get_normal_sequence(aligner : Align.PairwiseAligner, normal_act : list) -> dict:
    
    print("Processing normal sequences: ")

    with open(Path("matrix.csv"), 'w') as csv_file:

        csv_writer = csv.writer(csv_file, dialect='excel')

        # Create csv headers
        csv_writer.writerow(["id"] + [seq.id for seq in normal_act] + ["treshold", "sum"])

        # Ideal sequence attributes
        seq      = None 
        treshold = None

        max_sum = 0
        
        # Calculate alignment with each normal sequence
        for seq_rec in normal_act:       

            align_info = calculate_vector_alignment(aligner, seq_rec, normal_act)
            # Dump align list in csv file
            csv_writer.writerow(align_info.to_csv())

            if align_info.get_sum() > max_sum:
                seq      = align_info.get_seqrec()
                treshold = align_info.get_treshold()   

    return {
        "seq"      : seq,
        "treshold" : treshold
    }

# Align 'seq_rec' with each sequence in 'seq_vector'
def calculate_vector_alignment(aligner : Align.PairwiseAligner, seq_rec : SeqRecord, seq_vector : list) -> AlignInfo:

    # Only for printing loading line
    print(f"{' '*(20 + 10)}]", end='', flush=True)
    print(f"\r{seq_rec.id}\t[", end='', flush=True)

    score_list = []
    for index, list_item in enumerate(seq_vector):
            score_list.append({
                "id"    : list_item.id,
                "score" : aligner.score(seq_rec.seq, list_item.seq)
            })
            # Only for printing loading line
            print("#", end='', flush=True) if index % int(len(seq_vector)/20) == 0 else {}

    print("")

    # Calculate sum of align scores
    score_sum = sum(item["score"] for item in score_list)

    treshold  = score_sum / len(score_list)

    return AlignInfo({
        "seq"        : seq_rec,
        "score_list" : score_list,
        "treshold"   : treshold,
        "sum"        : score_sum
    })





