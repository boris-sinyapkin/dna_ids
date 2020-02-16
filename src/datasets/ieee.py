"""
    This file represents implementation of interfaces for IoT IEEE dataset.

    Link: https://ieee-dataport.org/open-access/iot-network-intrusion-dataset
"""
import json
import csv
import regex

from .interfaces     import JSON_Codetable, Dataset, DatasetRecord
from pathlib         import Path
from pandas          import DataFrame
from numpy           import array as numpy_array

from Bio.Seq         import Seq
from Bio.SeqRecord   import SeqRecord

class IEEE_Dataset(Dataset):
    def __getitem__(self, index):
        return IEEE_DatasetRecord(self.loc[index, :].tolist())

    def as_DNA_records(self, codetable : JSON_Codetable) -> numpy_array:
        result = []
        for id, row in self.iterrows():
            seq        = IEEE_DatasetRecord(row.tolist())
            dna_seq    = seq.encode_into_DNA(codetable)
            dna_seq.id = id

            result.append(dna_seq)

        return numpy_array(result)

class IEEE_DatasetRecord(DatasetRecord):
   # Encoding KDD record into DNA sequence
    def encode_into_DNA(self, codetable : JSON_Codetable):
        if not codetable.data:
            raise Exception("Codetable is empty")
        else:

            raw_record  = self.record
            payload     = ""

            for field in raw_record[:-1]:
                for literal in list(str(field)):
                    if regex.match(r'\d{1,}', literal):
                        payload += codetable["digits"][literal]
                    else:
                        payload += codetable[literal]

        return SeqRecord(Seq(payload), name=raw_record[-1], description=f"IEEE record encoded into DNA.")