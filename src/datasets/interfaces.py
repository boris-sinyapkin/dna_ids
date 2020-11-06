import json
import csv
import pandas

from numpy         import array as numpy_array
from pandas        import DataFrame
from pathlib       import Path
from abc           import ABC, abstractmethod

from Bio.SeqRecord import SeqRecord
from Bio.Seq       import Seq
from Bio           import Align

class AlignInfo:
    _info = {}

    sum       = property()
    seqrec    = property()
    threshold = property()

    def __init__(self, align_info : dict):
        self._info = align_info

    @sum.getter
    def sum(self):
        return self._info["sum"]
    @seqrec.getter
    def seqrec(self):
        return self._info["seq"]
    @threshold.getter
    def threshold(self):
        return self._info["threshold"]

class Codetable(ABC):
    """
        This class represents the DNA codes table, 
        which using when DatasetRecord's encode into DNA sequences.
    """
    _codetable = None
    data       = property()

    def __init__(self, path : Path):
        self.data = path
    
    @data.setter
    def data(self, path : Path):
        print("Loading the codetable from : %s" % (path))
        self._codetable = self.parse_codetable_file(path)

    @data.getter
    def data(self):
        return self._codetable

    @abstractmethod
    def parse_codetable_file(self, path : Path):
        pass

    @abstractmethod
    def __getitem__(self, value):
        pass

# Codetable loaded from *.json file. 
class JSON_Codetable(Codetable):
    def parse_codetable_file(self, path : Path):
        retval = {}
        with open(path, 'r') as json_codetable:
            retval = json.load(json_codetable)
        return retval

    def __getitem__(self, value):
        return self._codetable[value]

class DatasetRecord:
    """
        This class represents the 'Dataset' record. 
        Reimplement concrete methods to create your own logic of this class.
    """
    _record = None
    _label  = None
    record  = property()
    label   = property()

    def __init__(self, data):
        self.record = data[:-1]
        self.label  = data[-1]

    @record.setter
    def record(self, data):
        self._record = self.create_record(data)
        
    @label.setter
    def label(self, data):
        self._label = data

    @record.getter
    def record(self):
        return self._record
    
    @label.getter
    def label(self):
        return self._label

    def create_record(self, data):
        """
            This method should parse raw data and create record in concrete format.
        """
        return data

    @abstractmethod
    def encode_into_DNA(self, codetable : Codetable) -> SeqRecord:
        """
            This method should encode record from the Dataset into 
            DNA sequence using concrete codetable.
        """
        raise NotImplementedError
    
class Dataset(ABC, DataFrame):
    """
        This is the interface of 'Dataset'. 
        Derive this class and reimplement some of methods to create your own dataset. 
    """
    width  = property()
    height = property()
    
    def __init__(self, df : DataFrame):
        super().__init__(df)

    @staticmethod
    def from_file(path : Path):
        """
            This method should load raw dataset records from file, 
            parse and store there into the Dataset instance.

            Note: Basically, use this steps:
                df = pandas.read_csv(path)
                return Dataset(pandas.read_csv(path))
        """
        raise NotImplementedError
    
    @width.getter
    def width(self):
        return self.shape[1]
    
    @height.getter
    def height(self):
        return self.shape[0]
        
    @abstractmethod
    def raw_index(self, index) -> DatasetRecord:
        """
            Basically, should return DatasetRecord instance
        """
        raise NotImplementedError

    @abstractmethod
    def as_DNA_records(self, codetable : JSON_Codetable) -> numpy_array:
        """
            This method should return list of encoded activities.
        """
        raise NotImplementedError

    @abstractmethod
    def random_sample(self, size : int):
        """
            This method should return random sample of DatasetRecords

            Note: Basically, DataFrame.sample(...) method is used and result
                  is wrapped into Dataset class
        """
        raise NotImplementedError
    
    @abstractmethod
    def get_median(self) -> DatasetRecord:
        """
            This method should return a DatasetRecord instance. 
            Each element in vector (except "label") should be a median.
        """
        raise NotImplementedError