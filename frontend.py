
import argparse

from run     import run
from sys     import platform
from pathlib import Path

def main():

    # Temporary
    assert ('linux' in platform)

    parser = argparse.ArgumentParser(
        description="DNA Intrusion Detection System")

    parser.add_argument('--dataset',    '-d',   type=Path,  required=True, 
                        help="Path to KDD dataset. [*.csv]")
    parser.add_argument('--codetable',  '-c',   type=Path,  required=True,
                        help="Path to codetable. (Using for encoding dataset records in DNA sequences. [*.json])")
    parser.add_argument('--normal',     '-n',   type=Path,  required=False, default=None, 
                        help="Path to FASTA format file with normal activities.[*.faa]")

    args = parser.parse_args()

    # Execute the main IDS function
    run(args.dataset, args.codetable, args.normal)

if __name__ == '__main__':
    main()
