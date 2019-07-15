
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

# This function search normal sequence in list of normal activities             
def get_normal_sequence(aligner : Align.PairwiseAligner, normal_act : list, dump_matrix_in_csv=False) -> dict:
    
    if dump_matrix_in_csv:
        csv_file   = open(Path("matrix.csv"), 'w')
        csv_writer = csv.writer(csv_file, dialect='excel')

    # The selection criterion is the minimum of align scores sum
    sum_info = {
        "sum_list"  : [],
        "tresholds" : []
    }

    print("Processing normal sequences: ")

    # For each normal seq
    for seq_rec in normal_act:       
        # Calculate alignment with each normal sequence
        align_info = calculate_vector_alignment(aligner, seq_rec, normal_act)

        # Dump align list in csv file
        csv_writer.writerow(align_info["score_list"] + \
                            ["treshold", align_info["treshold"]] + \
                            ["sum",      align_info["sum"]]) if dump_matrix_in_csv else {}

        sum_info["sum_list"].append(align_info["sum"])
        sum_info["tresholds"].append(align_info["treshold"])
        print("")

    min_idx = sum_info["sum_list"].index(min(sum_info["sum_list"]))
    
    csv_file.close() if dump_matrix_in_csv else {}
    return {
        "seq_rec"  : normal_act[min_idx],
        "treshold" : sum_info["tresholds"][min_idx]
    }

# Align 'seq_rec' with each sequence in 'seq_vector'
def calculate_vector_alignment(aligner : Align.PairwiseAligner, seq_rec : SeqRecord, seq_vector : list) -> dict:

    # Only for printing loading line
    print(f"{' '*(20 + 10)}]", end='', flush=True)
    print(f"\r{seq_rec.id}\t[", end='', flush=True)

    score_list = []
    for index, list_item in enumerate(seq_vector):
            score_list.append(aligner.score(seq_rec.seq, list_item.seq))
            # Only for printing loading line
            print("#", end='', flush=True) if index % int(len(seq_vector)/20) == 0 else {}

    score_sum = sum(score_list)
    treshold  = score_sum / len(score_list)

    return {
        "score_list" : score_list,
        "treshold"   : treshold,
        "sum"        : score_sum
    }





