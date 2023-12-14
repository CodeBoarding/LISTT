# Lyme-In-Silico-Typing-Tool (LISTT)
a computational pipeline for typing the Borrelia OspA protein

## Installation
Clone the repository and enter the directory:

    git clone https://github.com/Pfizer-rd/Lyme-In-Silico-Typing-Tool.git
    cd Lyme-In-Silico-Typing-Tool

Build and activate the environment:

    conda env create -f environment.yml
    conda activate ospa_ist

## Usage

LISTT takes an assembled _Borrelia_ genome or Illumina NGS reads of the _ospA_ gene as input:

    python ospa_typing.py -r1 [read1 fastq] -r2 [read2 fastq] -p [path to fastq files] -m r

    python ospa_typing.py -g [genome fasta] -p [path to genome] -m a

Analysis will fail if there is insufficient coverage of _ospA_ or if _ospA_ is truncated in the genome assembly.



## Additional
Test genomes are sourced from the PubMLST genome sequence database

For more information, contact jonathan.lee@pfizer.com
