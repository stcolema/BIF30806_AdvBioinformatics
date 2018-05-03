#!/usr/bin/env python3
"""
Salmon indexing and quantifying of RNA-seq reads

Author: Stephen Coleman
"""

from subprocess import check_call
from os.path import exists
from argparse import ArgumentParser
from pathlib import Path
import tempfile
import csv

def argument_details():
    # argument parsing
    AP = ArgumentParser(description='This script is used to first index, and \
    then quantify a csv file of RNA-seq reads library by pseudomapping using Salmon.')
    AP.add_argument('filename', help='file to remove all whitespace from')
    AP.add_argument('output', help='file to write to')
    args = AP.parse_args()
    return args

#def check_arguments(args):
#    # Check if files are correctly given
#    if (args.unmatedReads != False and args.mates != False) or \
#        (args.unmatedReads == False and args.mates == False):
#        print('Give either 1 single end file, or 2 paired end files!')
#        exit()
#              
#    # make a file name from the given path
#    if args.unmatedReads != False:
#        path = args.unmatedReads
#    else:
#        path = args.mates[0]
#    sample = path.split('/')[-1].split('.')[0].split('_')[0]
#    if 'trim' in path:
#        sample += 'trim'
#    
#    # check if bootstrapping is done, change the file name
#    if args.bootstrap != 0:
#        sample += '_bs'+args.bootstrap
#    return sample

def main():
    args = argument_details()
    with open(args.filename, 'r') as fin:
       with open(args.output, 'w') as fout:
          for line in fin:
              if line.strip():
                  fout.write(line)
    return None

if __name__ == "__main__":
    main()
