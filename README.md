# ARMS-PCR Primer Design Tool

## Overview
The ARMS-PCR Primer Design Tool is a Python-based bioinformatics application designed to generate highly specific primers for Amplification Refractory Mutation System (ARMS) PCR. 

## Installation Requirements:
Python 3.9+
Primer3
Required Python packages: sys, os, subprocess, re, json, argparse, datetime

## Usage:

python ARMS_complete.py [--help] [--version] [--fasta FASTA] [--length LENGTH] [--nomask] [--param PARAM] [--output OUTPUT]

**example:**
```
python ARMS_complete.py --fasta test.fasta --output result.html
```

**optional arguments:**
| Arguments      | Note |
| -------------  | -----------------------------------------------------------------------------------------------|
| -h, --help     | #show this help message and exit                                                               |
|  --version     |#show program's version number and exit                                                         |
| --fasta FASTA  | #Input FASTA file (required)                                                                   |
|--length LENGTH |#Length of PCR template (default: 1001)                                                         |
|--nomask        |#Do not mask common SNPs in template                                                            |
|--param PARA    |#Custom parameters from JSON file (default: parameterDefault.json)                              |
|--output OUTPUT |#Output HTML file                                                                               |


## Input Requirements
The input FASTA file must contain exactly one SNP site formatted as [WildType/Mutation], e.g., [A/G]. Example:
```
>test
TAACTATACTAAGGCGTGACGTCATCTAAACTATGTGCCTCCGGGGTCCGTGCGCGTTGTTAGCTCCTTTACCACCACTATTACGAGTCGCCTAAAGAGGACAACCCGCTCAAAATGA
CTGGGTACGCTATTTCTACTTAGGCCCCACAGGATCTCCAAGCATATCAGTGTGGGGTATCACACTAGACTGGGGTGATGCCGCCGCCGCGTCGCGGGTGTTGATATAAATAGCCATT
GGCATTTTTGGGTTAGAGGCCGGCCAGCTTCAGGAGTCCTATGCATAACACTTTAAGTTTGGCTCGTTGTCACCCCCCCACGTCTTGCCTGGCAGAATTGATTCTACTCCCTTAATAC
GTAGAGTGCTCAGTCCTATCTCTCTTCAAGCCCGTGTCTAAACAGAAGATGCGGGACGTGCCTACCCGTCATCTCATTTGCGCCCCAAACGTTTCTGAGAGACGTCTTCTAGTTCTTC
ACAGAGTAAAGTAGTTAAGTTCTTGTAC[G/C]CGCCATTGTATTCATGGTACGCAGGCTCAAGACGCGCGCCATCCCCAAATTCTCCGAAGGGTAATCGGTCAATAGCGCTCGGCAC
CCCTCTCAGTTGACTATAGATTTAAGAGAAGTCACCATCCCTTCGGGCAGCCGGATTCCAGCTAGTGCGAGACAATTTTAAGAGAATACACGATGGTAGAGGTCGGGCCTAACGCTAT
AGGTTCGTAGAGTAAACATAGTGTAGGTTGTAGGATCCCCCAGACCCTTGGGTCCATATCAACCAACATGCTTAACATTAAGCGCCGAAGGCGGAAACCGCCTGATATCTGTTTTTGG
TGGAAGGCGTGAAAAAGAGATGTTACGTCAGTCTCGCGCACCCGTATATCGCAGTAAATCCGTGAGTTTACCCGCGTTTGCGTTTCATTAGATTACCACGTCGGGTACACAACCTACT
CAACATTTGGTACCTGTACAGGAGCTTCACTATCATTATAGTGAACGAGATGTTCCAAAGG
```


## Cite
The main content of this work was not independently completed. Here, [ABM-PCR](https://www.degruyterbrill.com/document/doi/10.1515/cclm-2024-0962/html)and [ARMSprimer3](https://pubmed.ncbi.nlm.nih.gov/40421167/) are integrated. 
