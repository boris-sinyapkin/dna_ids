from pathlib import Path
from abc     import ABC, abstractmethod

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

class Dataset(ABC):
    """
        This is the interface of 'Dataset'. 
        Derive this class and reimplement some of methods to create your own dataset. 
    """

    # List of DatasetRecords's
    _dataset   = []
    dataset    = property()

    def __init__(self, dataset_path : Path):
        self.dataset = dataset_path

    def __getitem__(self, index):
        return self._dataset[index]

    def __len__(obj):
        return len(obj._dataset)

    @dataset.setter
    def dataset(self, path : Path):
        print("Loading the dataset from : %s" % (path))
        self._dataset = self.parse_dataset_file(path)

    @abstractmethod
    def parse_dataset_file(self, path : Path) -> list:
        """
            This method should load raw dataset records from file, 
            parse and store there into the list as DatasetRecord's.
        """
        pass

    @abstractmethod
    def get_normal_seq(self, codetable : Codetable) -> list:
        """
            This method should return list of encoded normal activities.
        """
        pass

class DatasetRecord(ABC):
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

    @abstractmethod
    def create_record(self, data):
        """
            This method should parse raw data and create record in concrete format.
        """
        pass

    @abstractmethod
    def encode_into_DNA(self, codetable : Codetable):
        """
            This method should encode record from the Dataset into 
            DNA sequence using concrete codetable.
        """
        pass