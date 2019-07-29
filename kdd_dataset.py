import regex

from interfaces         import Dataset
from Bio.Seq            import Seq
from Bio.SeqRecord      import SeqRecord
from sys                import stderr

class KDD_Dataset(Dataset):
    
    # Encoding KDD record into DNA sequence
    def encode_record(self, index: int):
        if not self._codetable:
            raise Exception("Codetable is empty")
        else:
            raw_record = self._dataset[index]

            payload = ""
            label = ""
            extra = []

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

                        payload += prefix + \
                            self._codetable["hardcoded"][ident][item]
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
            id=str(index),
            name=label,
            features=extra,
            description=f"({index}) KDD record encoded into DNA. [{label}]"
        )

    # Search all sequences, which applies function in parameter list.
    def findall(self, value: str, func) -> list:
        retlist = []
        dataset_len = len(self._dataset)

        print("Processing sequences: ")
        print("                      ]\r[", end='', flush=True)

        for index in range(0, dataset_len):

            seq_rec = self.encode_kdd_record(index)

            if func(value, seq_rec):
                retlist.append(seq_rec)

            print("#", end='', flush=True) if index % int(
                dataset_len/20) == 0 else {}

        print("")
        return retlist
