# Lyme-In-Silico-Typing-Tool (LISTT)
a computational pipeline for typing the *Borrelia* OspA protein

LISTT assigns *in-silico* types based on sequence of the OspA protein. Illumina NGS reads of a *Borrelia* isolate or its genome sequence can be provided as input.
23 *in-silico* types (ISTs) can currently be assigned, with ISTs 1-8 corresponding to the canonical OspA serotypes 1-8: 

| Genospecies  | OspA *in-silico* type |
| ------------- | ------------- |
| *B. burgdorferi s.s.* | 1  |
| *B. afzelii* | 2  |
| *B. garinii* | 3, 5, 6, 7, 8, 11, 12 |
| *B. bavariensis* | 4, 9, 10 |
| *B. spielmanii* | 13|
| *B. mayonii* | 14 |
| *B. valaisiana* | 15, 25 |
| *B. turdi* | 16, 24 |
| *B. yangtzensis* | 17 |
| *B. americana* | 18 |
| *B. bissettiae* | 19 |
| *B. carolinensis* | 20 |
| *B. finlandensis* | 21 |
| *B. japonica* | 22 |
| *B. lusitaniae* | 23 |


The full method for IST assignment is described in the publication: [https://pubmed.ncbi.nlm.nih.gov/38787376/](https://pubmed.ncbi.nlm.nih.gov/38787376/)

## Installation
Clone the repository and enter the directory:

    git clone https://github.com/Pfizer-opensource/LISTT.git
    cd LISTT

Build and activate the environment:

    conda env create -f environment.yml
    conda activate ospa_ist

## Usage

LISTT takes an assembled _Borrelia_ genome or Illumina NGS reads of the _ospA_ gene as input:

    python ospa_typing.py -r1 [read1 fastq] -r2 [read2 fastq] -m r

    python ospa_typing.py -g [genome fasta] -m a

Users may also specify the input file(s) location [-p] if an absolute path is not provided.

Analysis will fail if there is insufficient coverage of _ospA_ or if _ospA_ is truncated in the genome assembly.



## Additional
Test genomes are sourced from the PubMLST genome sequence database

For more information, contact jonathan.lee@pfizer.com
