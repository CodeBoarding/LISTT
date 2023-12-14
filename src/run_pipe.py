# -*- coding: utf-8 -*-
"""
Python 3.6.13
12 April 2023; jonathan.lee@pfizer.com
@author: LEEJ322

"""

import os
import sys
from subprocess import call, DEVNULL, STDOUT

from src.ospA_aligner import ospA_Aligner
import src.pull_qual as pull_qual
from src.qc_alignments import BestAlignment
from src.build_consensus import ConsensusBuilder
from src.align_to_refs import VariantSimilarity
from src.blast_align import OspABLASTer
from src.alignment_metrics import CoverageCheck

def ngs_mode(read1, read2, dataPath):

    # make aln output folder and sample subdirectory
    sample = read1.split('_')[0]
    if not os.path.exists('{}/aln'.format(dataPath)):
        call('mkdir {}/aln'.format(dataPath), stdout = DEVNULL, stderr = STDOUT, shell = True)
        if not os.path.exists('{}/aln'.format(dataPath)):
            raise ValueError('run:run_pipe() - ERROR: '+'{}/aln'.format(dataPath)+' does not exist')
    if not os.path.exists('{0}/aln/{1}'.format(dataPath, sample)):
        call('mkdir {0}/aln/{1}'.format(dataPath, sample), stdout = DEVNULL, stderr = STDOUT, shell = True)
        if not os.path.exists('{0}/aln/{1}'.format(dataPath, sample)):
            raise ValueError('run:run_pipe() - ERROR: '+'{0}/aln/{1}'.format(dataPath, sample)+' does not exist')
    
    # run alignments to each reference sequence
    print('run:run_pipe() - Calling ospA_Aligner')
    ospA_Aligner(read1, read2, dataPath)
    
    # pull quality of alignments
    indir = '{0}/aln/{1}'.format(dataPath, sample)

    vcf_files = [x for x in os.listdir(indir) if x.endswith('_var.vcf')]
    for vcf_file in vcf_files:
        print('run:run_pipe() - Extracting alignment quality for' + vcf_file)
        pull_qual.extract_quality(sample, vcf_file, dataPath)
    
    # determine best alignment
    print('run:run_pipe() - Checking Best')
    best = BestAlignment(sample, dataPath).best

    if best: # only continue if alignment was successful
        # check alignment qc metrics
        pass_alignment_qc = CoverageCheck(sample, dataPath).pass_qc
    else: pass_alignment_qc = False

    if pass_alignment_qc: 
        #build consensus sequence
        best_vcf = '{0}_ospA_{1}_var.vcf'.format(sample, best)
        best_ref = 'ospA_{}.fasta'.format(best)
        print('run:run_pipe() - Best alignment to ' + best)
        ConsensusBuilder(best_vcf, best_ref, dataPath)
        consensus = sample + '_prot.fasta'

        # call serotype and species and write out result
        print('run:run_pipe() - VariantSimilarity 1')
        print('  Consensus: '+consensus)
        print('  dataPath: '+dataPath)
        result = VariantSimilarity(consensus, dataPath)
        print('  Result: ', result)
    else:
        print('run:run_pipe() - VariantSimilarity 2')
        print('  Sample: '+ sample)
        print('  dataPath: '+ dataPath)
        result = VariantSimilarity(sample, dataPath, False)
        print('  Result: ', result)

    
def assembly_mode(genome, dataPath):
    # make aln output folder and sample subdirectory
    sample = genome.split('.f')[0]
    if not os.path.exists('{}/aln'.format(dataPath)):
        call('mkdir {}/aln'.format(dataPath), stdout = DEVNULL, stderr = STDOUT, shell = True)
        if not os.path.exists('{}/aln'.format(dataPath)):
            raise ValueError('run:run_pipe() - ERROR: '+'{}/aln'.format(dataPath)+' does not exist')
    if not os.path.exists('{0}/aln/{1}'.format(dataPath, sample)):
        call('mkdir {0}/aln/{1}'.format(dataPath, sample), stdout = DEVNULL, stderr = STDOUT, shell = True)
        if not os.path.exists('{0}/aln/{1}'.format(dataPath, sample)):
            raise ValueError('run:run_pipe() - ERROR: '+'{0}/aln/{1}'.format(dataPath, sample)+' does not exist')

    alignment = OspABLASTer(genome, dataPath)
    consensus = sample + '_prot.fasta'
    if alignment.prot:
        print('ospA sequence aligned. Calling Serotype')
        result = VariantSimilarity(consensus, dataPath)
    else: 
        print(alignment.flag)
        result = VariantSimilarity(consensus, dataPath, False)
    

if __name__ == '__main__':
    test()
