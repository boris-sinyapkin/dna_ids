
import csv
import json
import regex

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
    def encode_kdd_record(self, index : int) -> dict:
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
                    # hardcoded value
                    ident   = self._codetable["prefixes"][f"{idx}"][0]
                    prefix  = self._codetable["prefixes"][f"{idx}"][1]

                    payload += prefix + self._codetable["hardcoded"][ident][item]     

            label = raw_record[41]

            for i in range(42, len(raw_record)):
                extra.append(raw_record[i])
                
        return {
                "payload" : payload,
                "label"   : label,
                "extra"   : extra
        }