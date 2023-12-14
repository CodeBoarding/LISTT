# -*- coding: utf-8 -*-
"""
Python 3.6.13
12 October 2023; jonathan.lee@pfizer.com
@author: LEEJ322

"""
import os

from src import cmd_parse, run_pipe

# verify either ngs or assembly mode or print an error if missing inputs
def verify_args(args):
    valid_args = True
    # check that a valid mode is provided
    if args['mode'] == 'r':
        if not args['read1']:
            print('Error: Fastq file(s) required for reads mode. Use -r1 and -r2')
            valid_args = False
        print('reads ')
    elif args['mode'] == 'a':
        if not args['genome']:
            print('Error: Fasta file required for assembly mode. Use -g')
            valid_args = False
    else:
        print('Error: invalid run mode. see --help for more details.')
        valid_args = False
    return valid_args

# choose which version of the pipeline to run
def choose_pipe(args):
    if args['mode'] == 'r': run_pipe.ngs_mode(args['read1'], args['read2'], args['path'])
    else: run_pipe.assembly_mode(args['genome'], args['path'])


def test():
    indir = 'test/pubmlst'
    fastas = [f for f in os.listdir(indir) if f.endswith('.fas')]
    for fasta in fastas:
            run_pipe.assembly_mode(fasta, indir)

def main():
    args = cmd_parse.CommandLine().args
    for arg in args: print(arg, args[arg])
    if verify_args(args): choose_pipe(args)


if __name__ == '__main__':
    #test()
    main()
