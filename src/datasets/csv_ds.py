"""
    This file represents implementation of interfaces for IoT IEEE dataset.

    Link: https://ieee-dataport.org/open-access/iot-network-intrusion-dataset
"""
import json
import csv
import regex
import pandas

from .interfaces     import JSON_Codetable, Dataset, DatasetRecord
from pathlib         import Path
from pandas          import DataFrame
from numpy           import array as numpy_array

from Bio.Seq         import Seq
from Bio.SeqRecord   import SeqRecord

class CSV_Dataset(Dataset):

    @staticmethod
    def from_file(path : Path):
        df = pandas.read_csv(path, low_memory=False)
        return CSV_Dataset(df.fillna(0))

    def __getitem__(self, index):
        relative_indexes = self.index
        absolute_index   = relative_indexes[index]
        return CSV_DatasetRecord(numpy_array(self.loc[absolute_index, :].tolist()))

    def as_DNA_records(self, codetable : JSON_Codetable) -> numpy_array:
        result = []
        for id, row in self.iterrows():
            seq        = CSV_DatasetRecord(row.tolist())
            dna_seq    = seq.encode_into_DNA(codetable)
            dna_seq.id = id
            result.append(dna_seq)

        return numpy_array(result)
    
    def random_sample(self, size : int):
        return CSV_Dataset(self.sample(size))
        
class CSV_DatasetRecord(DatasetRecord):

    from math import ceil

   # Encoding IEEE dataset record into DNA sequence
    def encode_into_DNA(self, codetable : JSON_Codetable):
        if not codetable.data:
            raise Exception("Codetable is empty")
        else:
            raw_record  = self.record
            payload     = ""

            for field in raw_record[:-1]:
                field = str(round(field, 10)) if type(field) is float else str(field)

                if field.startswith("0x"):
                    # Convert Hex value to int
                    field = str(int(field, base=16))

                for literal in list(field):
                    payload += codetable[literal]

        return SeqRecord(Seq(payload), name=raw_record[-1], description=f"IEEE record encoded into DNA.")