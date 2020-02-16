"""
    This script is using to format IEEE dataset *.csv file, obtained from *.pcap file
    
    1. *.pcap -> *.csv: 
            tshark -r <*.pcap> -q -z conv,tcp > <*.csv>
    
    2. *.csv -> dataset *.csv:
            python prepare.py -i <*.csv> -o <*.csv>
"""
import csv
import argparse
import re

from pathlib import Path

parser = argparse.ArgumentParser(description="DNA Intrusion Detection System")

parser.add_argument('--input',    '-i',   type=Path,  required=True, help="Path to dataset. [*.csv]")
parser.add_argument('--label',    '-l',   type=Path,  required=True, help="Traffic label")
parser.add_argument('--output',   '-o',   type=Path,  required=True)
args = parser.parse_args()

INPUT_DS  = args.input
OUTPUT_DS = args.output 
LABEL     = args.label

content = []
with INPUT_DS.open("r") as input_ds:
    content = input_ds.readlines()

CSV_HEADER = [ "ip_src", "ip_dst", "in_frames", "in_bytes", "out_frames", "out_bytes", "total_frames", "total_bytes", "rel_start", "duration", "label" ]

with OUTPUT_DS.open("w") as output_ds: 

    csv_writer = csv.writer(output_ds) 

    csv_writer.writerow(CSV_HEADER)

    for line in content[5:-1]:
        format_line = re.split(r' +', line.strip())
        format_line = [format_line[0]] + format_line[2:] + [ LABEL ]

        csv_writer.writerow(format_line)