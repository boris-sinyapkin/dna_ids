#!/bin/bash

echo "Warning: You need to have JRE 1.6 or higher installed."

echo "Download and unpack NSL KDD Dataset"
wget https://iscxdownloads.cs.unb.ca/iscxdownloads/NSL-KDD/NSL-KDD.zip

mkdir NSL-KDD

unzip NSL-KDD.zip -d NSL-KDD/

rm NSL-KDD.zip

echo "Install python modules"
pip install regex
pip install biopython

