
import argparse

from src.main        import run
from sys             import platform
from pathlib         import Path

def main():

    parser = argparse.ArgumentParser(description="DNA Intrusion Detection System")

    parser.add_argument('--train_dataset',      type=Path,  required=True, help="Path to test dataset. [*.csv]",  default=None)
    parser.add_argument('--test_dataset',       type=Path,  required=False, help="Path to train dataset. [*.csv]", default=None)     
    parser.add_argument('--codetable',  '-c',   type=Path,  required=True,
                        help="Path to codetable. (Using for encoding dataset records in DNA sequences. [*.json])")

    args = parser.parse_args()

    TRAIN_DS   = args.train_dataset
    TEST_DS    = args.test_dataset
    CODETABLE  = args.codetable

    # Execute the main IDS function
    run(TRAIN_DS, TEST_DS, CODETABLE)

if __name__ == '__main__':
    main()