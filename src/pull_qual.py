# -*- coding: utf-8 -*-
"""
Python 3.6.13
11 May 2021; jonathan.lee@pfizer.com
@author: LEEJ322

exports read depth and quality at each position from variant vcf file

python pull_qual.py [isolate]

"""
import os
import sys

import pandas as pd
import vcf


def extract_quality(sample, vcf_file, dataPath):
    muts = []
    infile = '{0}/aln/{1}/{2}'.format(dataPath, sample, vcf_file)
    print('pull_qual:extract_quality() - vcf.Reader:\n' + infile)

    if not os.path.isfile(infile):
        raise ValueError('pull_qual:extract_quality() - vcf.Reader: file does not exist\n' + infile)

    reader = vcf.Reader(open(infile, 'r'))
    for record in reader:
        dp = record.INFO['DP4']
        ref_count = dp[0] + dp[1]
        alt_count = dp[2] + dp[3]
        high_qual_depth = ref_count + alt_count
        if high_qual_depth > 0:
            purity = max([ref_count / high_qual_depth, alt_count / high_qual_depth])
        else:
            purity = 0
        outputs = (record.POS, record.REF, str(record.ALT[0]), record.INFO['DP'],
                   round(purity, 3))
        muts.append(outputs)

    df = pd.DataFrame(muts, columns=['pos', 'ref', 'alt', 'depth', 'purity'])
    outputFile = '{0}/aln/{1}/{2}_depth.csv'.format(dataPath, sample, vcf_file.split('_var.')[0])
    df.to_csv(outputFile, index=False)
