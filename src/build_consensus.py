# -*- coding: utf-8 -*-
"""
Python 3.6.13
22 April 2021; jonathan.lee@pfizer.com
@author: LEEJ322

python build_consensus.py [vcf variant file] [variant]

"""
import os
import sys
from subprocess import call, DEVNULL, STDOUT

from Bio import SeqIO
from Bio.Seq import Seq

class ConsensusBuilder:

    def __init__(self, vcf, ref, path):
        self.vcf = vcf
        self.ref = ref
        self.sample = vcf.split('_')[0]
        self.path = path
        self.run()

    def run_command(self, cmd):
        log_file = '{0}/aln/{1}/{1}_log.txt'.format(self.path, self.sample)
        log = open(log_file, 'a+')
        call(cmd, stdout = log, stderr = STDOUT, shell = True)
        log.close()

    def call_consensus(self, q = 30):
        outf = self.sample + '_consensus.fasta'
        # filter variants by quality (can add depth of coverage)
        cmd1 = "bcftools view -o {0}/aln/{1}/{1}_calls.vcf -i 'QUAL>={2}' {0}/aln/{1}/{3}".format(self.path, self.sample, q, self.vcf)
        # filter adjacent indels
        cmd2 = "bcftools filter --IndelGap 5 {0}/aln/{1}/{1}_calls.vcf -o {0}/aln/{1}/{1}_calls_indel-filt.vcf".format(self.path, self.sample)
        # compress and index vcf
        cmd3 = "bgzip {0}/aln/{1}/{1}_calls_indel-filt.vcf".format(self.path, self.sample)
        cmd4 = "bcftools index {0}/aln/{1}/{1}_calls_indel-filt.vcf.gz".format(self.path, self.sample)
        # build consensus sequence
        cmd5 = "bcftools consensus -f ref/{0} {1}/aln/{2}/{2}_calls_indel-filt.vcf.gz > {1}/aln/{2}/{3}".format(self.ref, self.path, self.sample, outf)
        cmds = [cmd1, cmd2, cmd3, cmd4, cmd5]
        for cmd in cmds:
            self.run_command(cmd)

    def translate_consensus(self):
        inf = '{0}/aln/{1}/{1}_consensus.fasta'.format(self.path, self.sample)
        for record in SeqIO.parse(inf, 'fasta'):
            seq = str(record.seq.translate())
        outfile = '{0}/aln/{1}/{1}_prot.fasta'.format(self.path, self.sample)
        with open(outfile, 'w') as handle:
            handle.write('>consensus\n')
            handle.write(seq + '\n')

    def run(self):
        print('calling consensus for ' + self.vcf)
        self.call_consensus()
        print('translating to protein sequence')
        self.translate_consensus()
