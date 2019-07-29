# DNA-based Intrusion Detection System

The main idea of this project is using bioinformatic features in detecting masquerading attacks. 

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Installing

```
git clone https://github.com/Dr-Livsey/dna_ids.git
cd dna_ids
./install_packages.sh
python frontend.py -h
```

### TODO
- [ ] Implement interfaces for Dataset and IDS classes
- [ ] Use "Builder" pattern when implementing the main IDS algo
- [ ] Create separate object for matrix.csv records

## Built With

* [Biopython](http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc32) - Main python module to work with DNA sequences
* [NSL-KDD](https://www.unb.ca/cic/datasets/nsl.html)
