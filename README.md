ARMS-PCR Primer Design Tool
Overview
The ARMS-PCR Primer Design Tool is a Python-based bioinformatics application designed to generate highly specific primers for Amplification Refractory Mutation System (ARMS) PCR. 

Installation Requirements:
Python 3.9+
Primer3
Required Python packages: sys, os, subprocess, re, json, argparse, datetime

Usage:
Command Line Arguments
plaintext
usage: ABM-PCR Primer Design Tool (Supports HTML Output) [-h] [--version] -f FASTA [--length LENGTH] [--nomask] [--param PARAM] [-o OUTPUT]

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  -f FASTA, --fasta FASTA
                        Input FASTA file (required)
  --length LENGTH       Length of PCR template (default: 1001)
  --nomask              Do not mask common SNPs in template
  --param PARAM         Custom parameters from JSON file (default: parameterDefault.json)
  -o OUTPUT, --output OUTPUT
                        Output HTML file path
Input Requirements
The input FASTA file must contain exactly one SNP site formatted as [WildType/Mutation], e.g., [A/G]. Example:

fasta
> test
TAACTATACTAAGGCGTGACGTCATCTAAACTATGTGCCTCCGGGGTCCGTGCGCGTTGTTAGCTCCTTTACCACCACTATTACGAGTCGCCTAAAGAGGACAACCCGCTCAAAATGACTGGGTACGCTATTTCTACTTAGGCCCCACAGGATCTCCAAGCATATCAGTGTGGGGTATCACACTAGACTGGGGTGATGCCGCCGCCGCGTCGCGGGTGTTGATATAAATAGCCATTGGCATTTTTGGGTTAGAGGCCGGCCAGCTTCAGGAGTCCTATGCATAACACTTTAAGTTTGGCTCGTTGTCACCCCCCCACGTCTTGCCTGGCAGAATTGATTCTACTCCCTTAATACGTAGAGTGCTCAGTCCTATCTCTCTTCAAGCCCGTGTCTAAACAGAAGATGCGGGACGTGCCTACCCGTCATCTCATTTGCGCCCCAAACGTTTCTGAGAGACGTCTTCTAGTTCTTCACAGAGTAAAGTAGTTAAGTTCTTGTAC[G/C]CGCCATTGTATTCATGGTACGCAGGCTCAAGACGCGCGCCATCCCCAAATTCTCCGAAGGGTAATCGGTCAATAGCGCTCGGCACCCCTCTCAGTTGACTATAGATTTAAGAGAAGTCACCATCCCTTCGGGCAGCCGGATTCCAGCTAGTGCGAGACAATTTTAAGAGAATACACGATGGTAGAGGTCGGGCCTAACGCTATAGGTTCGTAGAGTAAACATAGTGTAGGTTGTAGGATCCCCCAGACCCTTGGGTCCATATCAACCAACATGCTTAACATTAAGCGCCGAAGGCGGAAACCGCCTGATATCTGTTTTTGGTGGAAGGCGTGAAAAAGAGATGTTACGTCAGTCTCGCGCACCCGTATATCGCAGTAAATCCGTGAGTTTACCCGCGTTTGCGTTTCATTAGATTACCACGTCGGGTACACAACCTACTCAACATTTGGTACCTGTACAGGAGCTTCACTATCATTATAGTGAACGAGATGTTCCAAAGG
The tool generates:


The main content of this work was not independently completed. Here, ABM-PCR and ARMSprimer3 are integrated. Please cite：
Lau, Cia-Hin et al. “Artificial base mismatches-mediated PCR (ABM-PCR) for detecting clinically relevant single-base mutations.” Clinical chemistry and laboratory medicine vol. 63,7 1301-1314. 17 Mar. 2025
Guo, Jingwen et al. “ARMSprimer3: An open-source primer design Python program for amplification refractory mutation system PCR (ARMS-PCR).” Journal of pathology informatics vol. 17 100442. 15 Apr. 2025
