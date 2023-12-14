# -*- coding: utf-8 -*-
'''
Python 3.10.2
Nov 2023; jonathan.lee@pfizer.com
@author: LEEJ322

Compares quality of each alignment to references to select best
sequence for consensus building

'''

import csv
import os
import sys

import numpy as np
import pandas as pd

class BestAlignment:

    def __init__(self, sample, path):
        self.sample = sample
        self.qual_files = None
        self.thresholds = None
        self.best = None
        self.lengths = self.load_lengths()
        self.path = path
        self.run()

    # load reference lengths to check for sites with no coverage
    def load_lengths(self):
        length_tbl = {}
        file = 'variants/reference_lengths.csv'
        with open(file, 'r') as infile:
            reader = csv.reader(infile)
            next(reader)
            for row in reader:
                length_tbl[int(row[0])] = int(row[1])
        return length_tbl

    def get_files(self):
        indir = '{0}/aln/{1}'.format(self.path, self.sample)
        self.qual_files = [x for x in os.listdir(indir) if x.endswith('_depth.csv')]

    # compute mean and sd coverage - accounting for sites with no coverage
    def threshold_helper(self, cov_file):
        inf = '{0}/aln/{1}/{2}'.format(self.path, self.sample, cov_file)
        print('qc_alignments:threshold_helper() - inf: ' + inf)
        if not os.path.isfile(inf):
            raise ValueError('qc_alignments:threshold_helper() - Input file does not exist:\n' + inf)
        df = pd.read_csv(inf)
        seq_depths = list(df['depth'])
        reference = int(cov_file.split('_')[-2][6:])
        no_cov_sites = self.lengths[reference] - len(seq_depths)
        if no_cov_sites > 0:
            seq_depths = seq_depths + [0]*no_cov_sites
        mean = np.mean(seq_depths)
        sd = np.std(seq_depths)
        diff = mean - sd
        return diff

    # compute mean - sd coverage for each alignment
    def pull_thresholds(self):
        thresholds = {}
        print('qc_alignments:pull_thresholds() - qual_files: ', self.qual_files)
        for cov_file in self.qual_files:
            allele = cov_file.split('_')[-2]
            thresholds[allele] = self.threshold_helper(cov_file)
        # highest threshold used for building consensus
        if max(thresholds.values()) > 0: # accounting for cases of no coverage
            self.best = max(thresholds, key = thresholds.get)
        dft = pd.DataFrame.from_dict(thresholds, orient = 'index')
        dft.reset_index(inplace = True)
        dft.columns = ['variant', 'threshold']
        dft.sort_values('threshold', inplace = True, ascending = False)
        self.thresholds = dft

    def run(self):
        self.get_files()
        self.pull_thresholds()
        outfile = '{0}/aln/{1}/thresholds.csv'.format(self.path, self.sample)
        self.thresholds.to_csv(outfile, index = False)
