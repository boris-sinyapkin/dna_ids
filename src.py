
import os
import csv

import ids



def main():
    dataset_dir = os.path.abspath("NSL_KDD")
    dataset_filename = "KDDTrain+.csv"

    codetable_dir = os.path.abspath(".")
    codetable_filename = "dna_codes.json"

    dataset = ids.kdd_dataset(os.path.join(dataset_dir, dataset_filename))

    dataset.set_codetable(os.path.join(codetable_dir, codetable_filename))
    
    dataset.encode_kdd_record(0)



if __name__ == '__main__':
    main()