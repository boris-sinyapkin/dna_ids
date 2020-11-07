# DNA-based Intrusion Detection System

The main idea of this project is using bioinformatic features in detecting masquerading attacks. 

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Installing

```bash
git clone https://github.com/Dr-Livsey/dna_ids.git
cd dna_ids
pip3 install -r requirements.txt
python run.py -h
```

## Examples

### Run
Let you dataset located in *datasets/CSV/IEEE-IoT/dos-syn-flooding-1* folder. So, running options of this IDS suit will be the next:
```
py -3 run.py \
--train_dataset datasets/CSV/IEEE-IoT/dos-syn-flooding-1/train.csv \
--test_dataset  datasets/CSV/IEEE-IoT/dos-syn-flooding-1/test.csv \
--codetable     datasets/CSV/IEEE-IoT/dos-syn-flooding-1/codetable.json
```

If you want to do test analysis of the IDS use the "--test_dataset" option.

### Obtain machine learning metrics of IDS
To test this IDS prepare you datasets folders in the following way. 
For example, let we have 2 datasets to test. Firstly, put them into different folders, each of which should be contains 3 files:
- **train.csv** - Train subset
- **test.csv** - Test subset
- **codetable.json** - Will be used in process of encoding traffic into DNA sequences.

So, we have 2 folders. For example:
- **datasets/CSV/IEEE-IoT/dos-syn-flooding-1**
- **datasets/CSV/IoT/mirai-httpflooding-1**

To run IDS test suit for above folder structure, use:
```
py -3 scripts/ids_test.py --dir datasets
```
This run will produce metrics.csv file with IDS characteristics which will be located in each above folder.

### TODO
- [x] Implement interfaces for Dataset and IDS classes
- [x] Implement multithreaded searching of the best signature
- [x] Use "Builder" pattern when implementing the main IDS algo
- [ ] To implement Smith-Waterman local alignment
- [ ] Smith-Waterman & Needleman–Wunsch comparasion
- [ ] Add ROC-curve metric

### Datasets
- [IEEEDataPort IoT Network Intrusion Dataset](https://ieee-dataport.org/open-access/iot-network-intrusion-dataset)
- [Stratosphere Lab Malware on IoT Dataset](https://www.stratosphereips.org/datasets-iot23)
- [NSW Canberra The BoT-IoT Dataset](https://www.unsw.adfa.edu.au/unsw-canberra-cyber/cybersecurity/ADFA-NB15-Datasets/bot_iot.php)

## Built With

* [Biopython](http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc32) - Main python module to work with DNA sequences
