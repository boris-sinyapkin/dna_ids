import json
import csv
import pandas

from numpy         import array as numpy_array
from pandas        import DataFrame
from pathlib       import Path
from abc           import ABC, abstractmethod
from Bio.SeqRecord import SeqRecord

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
    record  = property()

    def __init__(self, data):
        self.record = data

    @record.setter
    def record(self, data):
        self._record = self.create_record(data)

    @record.getter
    def record(self):
        return self._record

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
    def __init__(self, dataset_path : Path):
        print("Loading the dataset from : %s" % (dataset_path))
        super().__init__(self.parse_dataset_file(dataset_path))

    @abstractmethod
    def __getitem__(self, index) -> DatasetRecord:
        """
            Basically, should return DatasetRecord instance
        """
        raise NotImplementedError

    def parse_dataset_file(self, path : Path) -> DataFrame:
        """
            This method should load raw dataset records from file, 
            parse and store there into the pandas DataFrame instance.
        """           
        return pandas.read_csv(path)

    @abstractmethod
    def as_DNA_records(self, codetable : JSON_Codetable) -> numpy_array:
        """
            This method should return list of encoded activities.
        """
        raise NotImplementedError