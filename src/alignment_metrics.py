# -*- coding: utf-8 -*-
'''
Python 3.11.6
4 December 2023; jonathan.lee@pfizer.com

'''

import csv
import os
import re
from subprocess import call, DEVNULL, STDOUT
import sys

import numpy as np
import pandas as pd



class AlignmentQuality:

    # needs sample name and path to alignment files
    def __init__(self, sample, dataPath, cov_threshold=5):
        self.sample = sample
        self.dataPath = dataPath #path to original fastq files
        self.ref = None
        self.cov_threshold = cov_threshold
        self.outputs = None
        self.run()

    def get_reference(self):
        inpath = '{0}/aln/{1}/thresholds.csv'.format(self.dataPath, self.sample)
        df = pd.read_csv(inpath)
        if max(df.threshold != 0):
            self.ref = df['variant'][0]

    # counts total reads based on batches in bwa mem log
    def parse_log(self):
        logpath = '{0}/aln/{1}/{1}_log.txt'.format(self.dataPath, self.sample)
        total_reads = 0
        print('alignment_metrics:parse_log() - logpath: ' + logpath)
        with open(logpath, 'r') as inf:
            for line in inf:
                if line.startswith('[M::mem_process_seqs] Processed'):
                    reads = int(re.findall(r'\d+', line)[0])
                    total_reads += reads
                if line.startswith('[mpileup]'):
                    break
        return total_reads

    def parse_flagstat(self):
        infile = open('{0}/aln/{1}/{1}_stats.txt'.format(self.dataPath, self.sample))
        content = infile.readlines()
        # total_reads = int(content[0].split(' ')[0])
        mapped_line = [i for i, line in enumerate(content) if line.split(' ')[3] == 'mapped'][0]
        aligned_reads = int(content[mapped_line].split(' ')[0])
        # aligned_reads = int(content[4].split(' ')[0])
        # total_reads = self.count_lines_fastq()

        total_reads = self.parse_log()

        try:
            perc_aligned = (aligned_reads / total_reads) * 100
        except:
            perc_aligned = 0
        infile.close()
        return [total_reads, aligned_reads, perc_aligned]

    def count_reads(self):
        # pull statistics with flagstat command
        cmd = 'samtools flagstat {0}/aln/{1}/{1}_ospA_{2}.srt.bam \
        > {0}/aln/{1}/{1}_stats.txt'.format(self.dataPath, self.sample, self.ref)
        print('alignment_metrics:count_reads() - cmd: ' + cmd)
        call(cmd, stdout=DEVNULL, stderr=STDOUT, shell=True)
        read_counts = self.parse_flagstat()
        return read_counts

    def pull_metrics(self):
        print('alignment_metrics:pull_metrics() - hardcoded alleles')  # hardcoded
        sero3_refs = {'allele13', 'allele15', 'allele17', 'allele22'}
        if self.ref in sero3_refs:
            max_pos = 806
        else:
            max_pos = 819
        infile = '{0}/aln/{1}/{1}_ospA_{2}_depth.csv'.format(self.dataPath, self.sample, self.ref)
        df = pd.read_csv(infile)
        # get median coverage and percent above 5 read depth
        depth_tbl = pd.Series(df.depth.values, index=df.pos).to_dict()
        for pos in range(1, 819):
            depth_tbl[pos] = depth_tbl.get(pos, 0)
        depth_array = np.asarray([v for k, v in depth_tbl.items() if k > 32 and k < max_pos])
        median_cov = np.median(depth_array)
        x5 = round(min(np.sum(depth_array >= 5) / len(depth_array), 1), 3)
        return [median_cov, x5]

    def run(self):
        self.get_reference()
        if self.ref:
            counts = self.count_reads()
            depths = self.pull_metrics()
        else:  # if alignments failed
            self.ref = 'allele1'  # check allele1 .bam to pull total reads
            counts = self.count_reads()
            depths = [0, 0]
        self.outputs = counts + depths
        # print(self.outputs)

class CoverageCheck:

    def __init__(self, sample, dataPath):
        self.sample = sample
        self.dataPath = dataPath
        self.min_med = None
        self.min_breadth = None
        self.min_reads = None
        self.metrics = None
        print('>>> alignment_metrics:__init__() - CoverageCheck: self')
        self.run()

    def set_minimums(self):
        print('alignment_metrics:set_minimums()- loading')
        fields = {}
        with open('variants/min_cov_metrics.csv', 'r') as inf:
            reader = csv.reader(inf)
            for row in reader:
                fields[row[0]] = float(row[1])
        self.min_med = fields['min_med']
        self.min_breadth = fields['min_breadth']
        self.min_reads = fields['min_reads']

    def run_checks(self):
        '''
        metrics: list - obtained from AlignmentQuality.outputs
            metrics[0] = total reads
            metrics[1] = aligned reads
            metrics[2] = percent reads aligned
            metrics[3] = median cov
            metrics[4] = 5x coverage breadth
        '''
        print(self.metrics)
        check1 = self.metrics[1] >= self.min_reads
        check2 = self.metrics[3] >= self.min_med
        check3 = self.metrics[4] >= self.min_breadth
        # if check1 + check2 + check3 < 3 and self.metrics[2] > 0:
        print("alignment_metrics: run_checks() - metrics min_reads: {0} >= {1}\n".format(self.metrics[1], self.min_reads))
        print("alignment_metrics: run_checks() - metrics min_med: {0} >= {1}\n".format(self.metrics[3], self.min_med))
        print("alignment_metrics: run_checks() - metrics min_breadth: {0} >= {1}\n".format(self.metrics[4],
                                                                                         self.min_breadth))
        allChecks = int(check1) + int(check2) + int(check3)
        print("alignment_metrics: run_checks() - allChecks: {0}\n".format(allChecks))
        if allChecks < 3: return False
        else: return True

    def run(self):
        self.set_minimums()
        self.metrics = AlignmentQuality(self.sample, self.dataPath).outputs
        self.pass_qc = self.run_checks()
