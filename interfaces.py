import json
import csv

from pathlib import Path
from abc     import ABC, abstractmethod

class Dataset:
    _dataset   = []
    _codetable = {}

    codetable  = property()
    dataset    = property()

    def __init__(self, dataset_path : Path, codetable_path : Path):
        if dataset_path.exists():
            self.codetable = codetable_path

        if codetable_path.exists():
            self.dataset   = dataset_path

    def __getitem__(self, index):
        return self._dataset[index]

    def __len__(obj):
        return len(obj._dataset)

    @codetable.setter
    def codetable(self, path : Path):
        print("Loading the codetable from : %s" % (path))
        with open(path, 'r') as json_codetable:
            self._codetable = json.load(json_codetable)

    @dataset.setter
    def dataset(self, path : Path):
        print("Loading the dataset from : %s" % (path))
        with open(path, 'r') as dataset_file:
            self._dataset = [line for line in csv.reader(dataset_file)]

    @codetable.getter
    def codetable(self):
        return self._codetable

    @dataset.getter
    def dataset(self):
        return self._dataset

    @abstractmethod
    def encode_record(self, index : int):
        pass

