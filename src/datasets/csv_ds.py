"""
    This file represents implementation of interfaces for CSV dataset.
"""

import pandas as pd
import numpy  as np

from tqdm import tqdm

from .interfaces     import JSON_Codetable, Dataset, DatasetRecord
from ..utils         import normalize_df
from pathlib         import Path
from numpy           import array as numpy_array

from Bio.Seq         import Seq
from Bio.SeqRecord   import SeqRecord

class CSV_Dataset(Dataset):
    FLOAT_FIELDS = [ 'tcp.time_delta' ]

    @staticmethod
    def from_file(path : Path):
        df = pd.read_csv(path, low_memory=False).fillna(0).astype({'tcp.len' : int, 
                                                'tcp.time_delta' : 'float64',
                                                'tcp.seq_raw' : 'int64',
                                                'tcp.ack_raw' : 'int64',
                                                'tcp.hdr_len' : 'int64',
                                                'tcp.flags'   : str,
                                                'tcp.window_size_value' : int,
                                                'tcp.checksum.status' : int,
                                                'tcp.urgent_pointer' : int,
                                                'tcp.options' : str})
        
        df['tcp.flags']   = df['tcp.flags'].apply(lambda x : int(x, 16))
        df['tcp.options'] = df['tcp.options'].apply(lambda x : int(x, 16) if x != '<MISSING>' else 0)
                
        return CSV_Dataset(df)

    def raw_index(self, index):
        relative_indexes = self.index
        absolute_index   = relative_indexes[index]
        return CSV_DatasetRecord(numpy_array(self.loc[absolute_index, :].tolist()))
        
    def get_median(self):
        def value_to_str(value, column):
            if column in CSV_Dataset.FLOAT_FIELDS:
                return str(round(value, 6))
            else:
                return str(int(value))
            
        columns    = [ column for column in self if column != "label" ]
        median_row = [ value_to_str(self[column].median(), column) for column in columns ]
        median_row += [ 'median' ]
        return CSV_DatasetRecord(median_row)

    def as_DNA_records(self, codetable : JSON_Codetable) -> pd.Series:
        result = []
        for id, row in tqdm(self.iterrows(), total=self.shape[0], desc='Encoding Dataset into DNA records'):
            seq        = CSV_DatasetRecord(row.astype(str).tolist())
            dna_seq    = seq.encode_into_DNA(codetable, str(id))
            dna_seq.id = id
            result.append(dna_seq)
        return pd.Series(result)
    
    def random_sample(self, size : int):
        return CSV_Dataset(self.sample(size))
        
class CSV_DatasetRecord(DatasetRecord):

   # Encoding IEEE dataset record into DNA sequence
    def encode_into_DNA(self, codetable : JSON_Codetable, id='', description=f"CSV dataset record encoded into DNA.") -> SeqRecord:
        if not codetable.data:
            raise Exception("Codetable is empty")
        else:
            raw_record  = self.record
            payload     = ""
            for field in raw_record:
                field = str(round(field, 6)) if type(field) is np.float64 else str(field)
                for literal in list(field):
                    payload += codetable[literal]
        return SeqRecord(Seq(payload), id=id, name=self.label, description=description)