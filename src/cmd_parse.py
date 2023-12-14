# -*- coding: utf-8 -*-
'''
Python 3.6.13
15 November 2022
'''
import argparse

class CommandLine:
    def __init__(self):
        self.args = {'read1':None, 'read2':None, 'genome':None, 'path':None, 'mode':None}
        parser = argparse.ArgumentParser(description = "Lyme OspA in-silico typing")
        parser.add_argument("-r1", "--read1", help = "Read 1 or interleaved FASTQ file", required = False, default = "")
        parser.add_argument("-r2", "--read2", help = "Read 2 FASTQ file (paired-end only)", required = False, default = "")
        parser.add_argument("-g", "--genome", help = "Assembled genome FASTA file", required = False, default = "")
        parser.add_argument("-p", "--path", help = "path to FASTQ files directory", required = False, default = "")
        parser.add_argument("-m", "--mode", help = "reads [r] or genome assembly [a]", required = False, default = "")

        argument = parser.parse_args()
        if argument.read1: self.args['read1'] = argument.read1
        if argument.read2: self.args['read2'] = argument.read2
        if argument.genome: self.args['genome'] = argument.genome
        if argument.path: self.args['path'] = argument.path
        if argument.mode: 
            self.args['mode'] = argument.mode
        else: 
            print('no mode selected - defaulting to reads')
            self.args['mode'] = 'r'
