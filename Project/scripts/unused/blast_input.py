#!/usr/bin/env python3
"""
Creates a fasta file of sequences for BLAST searching.

Uses the transcriptome and sleuth significant isoforms to find the sequences.

Author: Stephen Coleman
Student ID: 940-309-160-050
"""

from subprocess import check_call
from os import getcwd
from os.path import exists
from argparse import ArgumentParser
from pathlib import Path
import tempfile
import csv

def argument_details():
    # argument parsing
    AP = ArgumentParser(description='This script is implelements a pipeline of fastq_quality_trim, '
                                    + 'salmon and wasabi, using a csv file of samples to use and a '
                                    + 'given transcriptome for indexing')
    AP.add_argument('table', help='csv file of entries for blasting')
    AP.add_argument('transcriptome', help='used to find the relevant sequences the reads')
    AP.add_argument('output', help='fasta file of sequences from table')
    args = AP.parse_args()
    return args

def fasta_reader(filename):
    """Imports a FASTA file as a dict or list depending on user selection.

    Key inputs:
    filename --- the name of the FASTA file to be converted into a dict
    """
    result = {}
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                current = line[1:]
                result[current] = ''
            else:
                result[current] += line
    return result


def convert_csv_to_dict(csv):
    """Converts a csv table to a dict with table header as key.

    Kye inputs:
    csv --- csv filename
    """
    output = []
    with open(csv) as infile:
        for i, row in enumerate(infile):
            row = row.replace('"', '').strip()
            if row:
                entry = row.split(',')
                if i == 0:
                    entry[0] = 'Index'
                    key = entry
                    continue
                else:
                    my_dict = dict(zip(key, entry))
                output += [my_dict]
    return output

def write_fasta(record, output_name, mode = 'w', line_width = 80):
    """Write a FASTA file with option to append to exisitng file.

    Key inputs:
    record --- the dictionary the DNA sequence with sequence IDs as keys.  for entry in table_of_isoforms:
    output_name --- the name of the finished FASTA file sequence = transcript_fasta[entry['target id']
    mode --- the mode for writing to the FASTA file; 'w' writes a new file
    (overwriting if necessary) and 'a' allows appending to a pre-existing file.
    line_width --- the number of characters allowed on each line of the FASTA
    file.
    """
    assert mode == 'w' or mode == 'a', 'Mode must be \'w\' or \'a\'. '
    with open(output_name, mode) as f:
        for key in record.keys():
            title = '>{}\n'.format(key)
            f.write(title)
            range_sup = ((len(record[key]) // line_width) + 1) * line_width
            for i in range(0, range_sup + 1, line_width):
                f.write(record[key][i:min(i + line_width, len(record[key]))])
                f.write('\n')
    return


def find_sequence(table_of_isoforms, transcript_fasta):
    """Returns a dictionary of fasta sequences with keys as sequence titles.

    Key inputs:
    table_of_isoforms --- output from convert_csv_to_dict
    trancript_fasta --- dict with keys as FASTA entry names
    """
    output = {}
    for entry in table_of_isoforms:
        sequence = transcript_fasta[entry['target_id']]
        output[entry['target_id']] = sequence
    return output


def main():
    """Main program: read in arguments, create fasta file of sequences in table."""
    args = argument_details()
    transcriptome_fasta = fasta_reader(args.transcriptome)
    table_of_isoforms = convert_csv_to_dict(args.table)
    sequences_of_interest = find_sequence(table_of_isoforms, transcriptome_fasta)
    write_fasta(sequences_of_interest, args.output)

if __name__ == '__main__':
    main()
