#!/usr/bin/env python3
"""
Does a blastp search on the inputted database (or nr as default)

Author: Stephen Coleman
Student ID: 940-309-160-050
"""

from subprocess import check_call, check_output
from os import getcwd
from os.path import exists
from argparse import ArgumentParser
from pathlib import Path
import tempfile
import csv

def argument_details():
    # argument parsing
    AP = ArgumentParser(description='This script is does a BLASTP search for the inputted sequences')
    AP.add_argument('sequences', help='fasta file to blast')
    AP.add_argument('output', help='output of blast search')
    AP.add_argument('-db', '--database', help='database to blast against.',
                    default='/local/data/blastdb/nr')
    AP.add_argument('-outfmt', '--output_format', help='Format of blastp output.',
                    default=7)
    AP.add_argument('-o', '--overwrite', help='Overwrite output file is already exists.',
                    default=7)
    args = AP.parse_args()
    return args


def blastp(database, query, output_to_file = False, output_file = None,
          overwrite = True, outfmt = 7):
    """Carries out a blastp search using blastp tool.
    Key inputs:
    database --- a blast database file outputted from makeblastdb.
    query --- query FASTA file.
    output_to_file --- boolean defining output as file or variable
    output_file --- only used if output_to_file is True; filename of
    desired output file
    overwrite --- boolean for overwriting pre-exisitng output file of the
    same name
    outfmt --- format of blastp output.
    """
    if output_to_file:
        if exists(output_file) and not overwrite:
            return output_file
        cmd = 'blastp -db {} -query {} -outfmt {} -out {} -num_alignments 1 -evalue 1e-5'.\
               format(database, query, outfmt, output_file)
    else:
        cmd = 'blastp -db {} -query {} -outfmt {} -num_alignments 1 -evalue 1e-5'.format(
            database, query, outfmt)
    if output_to_file:
        check_call(cmd, shell=True)
        return output_file
    printed_output = check_output(cmd, shell=True)
    return printed_output

def parse_blastp_output(blastp_output):
    """Yield list of tuples of values from BLASTP output.
    Key inputs:
    blastp_output --- output file from blastp program.
    """
    with open(blastp_output) as fo:
        for line in fo:
            if not line.strip():
                continue
            line.replace('\\n','')
            if 'Query:' in line:
                query = line.partition('Query: ')[2]
            elif 'Fields:' in line:
                fields = line.partition('Fields: ')[2].split(', ')
            elif '#' not in line:
                entries = line.split()
                entries[0] = query
                zipped = dict(zip(fields, entries))
                yield zipped


def main():
    args = argument_details()
    blast_res = blastp(args.database, args.sequences, output_to_file=True,
                       output_file=args.output, outfmt=args.output_format,
                       overwrite=args.overwrite)
    print(blast_res)


if __name__ == '__main__':
    main()
