# -*- coding: utf-8 -*-
"""
Python 3.6.13
10 June 2021; jonathan.lee@pfizer.com
@author: LEEJ322

aligns one isolate to each reference

arguments: [isolate]
"""

import os
import sys
from subprocess import call, DEVNULL, STDOUT

class ospA_Aligner:

    def __init__(self, read1, read2, path):
        self.read1 = read1
        self.read2 = read2
        self.path = path
        self.sample = read1.split('_')[0]
        self.run()

    def run_command(self, cmd):
        log_file = '{0}/aln/{1}/{1}_log.txt'.format(self.path, self.sample)
        log = open(log_file, 'a+')
        print(cmd)
        call(cmd, stdout = log, stderr = STDOUT, shell = True)
        log.close()

    def get_out_prefix(self, ref):
        return '{0}_{1}'.format(self.sample, ref.split('.fasta')[0])

    def align_reads(self, prefix, ref):
        if self.read2 == None: # assumes interleaved fastq file if no read2 fastq is provided
            cmd = 'bwa mem -p ref/{0} {1}/{2} | samtools view -b -F 4 | samtools sort -o {1}/aln/{3}/{4}.srt.bam'.format(
                ref, self.path, self.read1, self.sample, prefix)
        else: # assumes paired end fastq files
            cmd = 'bwa mem ref/{0} {1}/{2} {1}/{3}| samtools view -b -F 4 | samtools sort -o {1}/aln/{4}/{5}.srt.bam'.format(
                ref, self.path, self.read1, self.read2, self.sample, prefix)
        print('aligning {0} to {1}'.format(self.sample, ref))
        print(cmd)
        if os.path.exists('{0}/aln/{1}.srt.bam'.format(self.path, self.sample)): # skip alignment if bam file exists already
            print('skipping alignment to {} because sorted bam exists'.format(ref))
        else:
            self.run_command(cmd)

    def build_bcf(self, prefix, ref):
        bam = prefix + '.srt.bam'
        outf = prefix + '_align.bcf'

        cmd = 'bcftools mpileup -o {0}/aln/{1}/{2} -A --indel-bias 0.75 --max-depth 2147483647 -L 2147483647 -f ref/{3} {0}/aln/{1}/{4}'.format(
            self.path, self.sample, outf, ref, bam)
        print('building bcf file for ' + bam)

        self.run_command(cmd)

    def filter_variants(self, outp):
        inf = outp + '_align.bcf'
        outf = outp + '_var.vcf'

        cmd = 'bcftools call -o {0}/aln/{1}/{2} -O v --ploidy 1 -c {0}/aln/{1}/{3}'.format(self.path, self.sample, outf, inf)
        print('calling variants for ' + inf)

        self.run_command(cmd)

    def run(self):
        refs = [f for f in os.listdir('ref') if f.endswith('.fasta')]
        print('aligning '+ self.sample)
        for ref in refs:
            out_prefix = self.get_out_prefix(ref)
            self.align_reads(out_prefix, ref)
            self.build_bcf(out_prefix, ref)
            self.filter_variants(out_prefix)
