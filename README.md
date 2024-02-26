# Lyme-In-Silico-Typing-Tool (LISTT)
a computational pipeline for typing the Borrelia OspA protein

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
