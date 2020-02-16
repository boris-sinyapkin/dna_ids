"""
    This file represents implementation of interfaces for KDD dataset.
"""
import regex
import json
import csv

from interfaces         import JSON_Codetable, Dataset, DatasetRecord
from pathlib            import Path

from Bio.Seq            import Seq
from Bio.SeqRecord      import SeqRecord
from sys                import stderr

class KDD_DatasetRecord(DatasetRecord):
    def create_record(self, data):
        # data represents csv line
        return data

    # Encoding KDD record into DNA sequence
    def encode_into_DNA(self, codetable : JSON_Codetable):
        if not codetable.data:
            raise Exception("Codetable is empty")
        else:
            raw_record = self.record

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
                            payload += codetable["digits"][literal]
                        else:
                            # dot
                            payload += codetable["dot"]

                else:
                    try:
                        # hardcoded value
                        ident = codetable["prefixes"][f"{idx}"][0]
                        prefix = codetable["prefixes"][f"{idx}"][1]

                        payload += prefix + \
                            codetable["hardcoded"][ident][item]
                    except KeyError:
                        print("Key Error: \n\
                               \r- idx    = %d \n\
                               \r- ident  = \'%s\' \n\
                               \r- item   = \'%s\'" % (idx, ident, item), file=stderr)

            label = raw_record[41]

            for i in range(42, len(raw_record)):
                extra.append(raw_record[i])

        return SeqRecord(
            Seq(payload),
            name=label,
            features=extra,
            description=f"KDD record encoded into DNA. [{label}]"
        )
        

class KDD_Dataset(Dataset):

    def get_normal_seq(self, codetable : Codetable) -> list:
        retlist = []
        dataset_len = len(self)

        print("Getting normal sequences: ")
        print("                      ]\r[", end='', flush=True)

        for index in range(0, dataset_len):
            seq_rec    = self[index].encode_into_DNA(codetable)
            seq_rec.id = str(index)

            if "normal" == seq_rec.name:
                retlist.append(seq_rec)

            print("#", end='', flush=True) if index % int(
                dataset_len/20) == 0 else {}

        print("")
        return retlist

    def parse_dataset_file(self, path : Path) -> list:
        retlist = []
        with path.open() as dataset_file:
            retlist = [KDD_DatasetRecord(line) for line in csv.reader(dataset_file)]
        return retlist