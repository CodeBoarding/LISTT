# -*- coding: utf-8 -*-
'''
Python 3.11.6
30 March 2023
'''

import os
import re
from subprocess import call, DEVNULL, STDOUT
import xml.etree.ElementTree as ET

from Bio.Seq import Seq
import screed

class OspABLASTer:

    def __init__(self, fasta, path):
        self.fasta = fasta
        self.path = path
        self.sample = fasta.split('.f')[0]
        self.blst = 'dummy.xml'
        self.flag = ''
        self.seq = None
        self.prot = None
        self.run()

    def blast(self):
        ref = 'alleles/ref_alleles.fasta'
        self.blst = '{0}/aln/{1}/{1}_blast.xml'.format(self.path, self.sample)
        cmd = 'blastn -db {0} -query {1}/{2} -out {3} -outfmt 5'.format(ref, self.path, self.fasta, self.blst)
        print(cmd)
        call(cmd, stdout = DEVNULL, stderr = STDOUT, shell = True)

    def pull_results(self):
        print('pulling results from BLAST search')
        tree = ET.parse(self.blst)
        root = tree.getroot()
        # path to blast results
        hp = './BlastOutput_iterations/Iteration/Iteration_hits/Hit/Hit_hsps/Hsp'
        # find best scoring alignment
        scorepath = hp + '/Hsp_score'
        scores = [ET.tostring(x) for x in root.findall(scorepath)]
        num_scores = [float(re.search('>(.+?)<', str(x)).group(1)) for x in scores]
        if num_scores:
            best = num_scores.index(max(num_scores))
            attributes = ('Hsp_hit-from', 'Hsp_hit-to', 'Hsp_qseq')
            outputs = []
            for attribute in attributes:
                path = '{0}/{1}'.format(hp, attribute)
                element = ET.tostring(root.findall(path)[best])
                output = re.search('>(.+?)<', str(element)).group(1)
                outputs.append(output)
            blast_result = {'start':int(outputs[0]),
                            'end':int(outputs[1]),
                            'seq':outputs[2]}
            match_length = blast_result['end'] - blast_result['start']
            print('ospA BLAST aligned {} bases'.format(match_length))
        else:
            print('no significant alignment')
            blast_result = {'start':0,
                            'end':1,
                            'seq':'AAA'}
        return blast_result

    # reverse complement sequence if needed and remove dashes
    def curate_sequence(self, blast_result):
        seq = blast_result['seq'].replace('-', '')
        if blast_result['start'] > blast_result['end']:
            self.seq = screed.rc(seq)
        else: self.seq = seq
        if len(seq) < 813: #require 99% sequence alignment
            self.flag = 'blast_align:OspABlaster() - Unable to obtain complete ospA sequence'
            self.seq = None

    def translate_and_write(self):
        if self.seq:
            self.prot = str(Seq(self.seq).translate())
        else: self.prot = ''
        outfile = '{0}/aln/{1}/{1}_prot.fasta'.format(self.path, self.sample)
        with open(outfile, 'w') as handle:
            handle.write('>consensus\n')
            handle.write(self.prot + '\n')

    # remove blast output files
    def cleanup(self):
        if os.path.exists(self.blst):
            os.system('rm ' + self.blst)

    def run(self):
        try:
            self.blast()
            blast_result = self.pull_results()
            corrected = self.curate_sequence(blast_result)
            self.translate_and_write()
        except:
            self.flag = 'BLAST alignment failed'
        finally:
            #self.cleanup()
            pass