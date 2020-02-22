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
Let be you dataset located in Datasets/IoT/IEEE folder. So, running options of this IDS suit will be the next:
```python
py -3 run.py \
--train_dataset Datasets/IoT/IEEE/train.csv \
--codetable     iot_dna.json
```

If you want to do test analysis of the IDS use the "--test_dataset" option.

### Obtain machine learning metrics of IDS
To test this IDS prepare you datasets folders in the following way. 
For example, let we have 2 datasets to test. Firstly, put it on different folders, which must contain 3 files:
- **train.csv** - Train subset
- **test.csv** - Test subset
- **codetable.json** - Will be use in process of encoding traffic into DNA sequences.

So, we have 2 folders. For example:
- **Datasets/IEEE/dos_synflooding-1**
- **Datasets/IoT/mirai-httpflooding-1**

To run IDS test suit for above folder structure, use:
```python
py -3 scripts/ids_test.py --dir Datasets
```
This run will produce metrics.csv file with IDS characteristics which will be located in each above folder.

### TODO
- [x] Implement interfaces for Dataset and IDS classes
- [x] Implement multithreaded searching of the best signature
- [x] Use "Builder" pattern when implementing the main IDS algo
- [ ] Enhance traffic encoding method
- [ ] Parse input traffic and learn the IDS "on the fly"


## Built With

* [Biopython](http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc32) - Main python module to work with DNA sequences
