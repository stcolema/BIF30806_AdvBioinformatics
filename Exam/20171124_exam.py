#!/usr/bin/env python3

"""
Author: Stephen Coleman
Student number: 940309-160-050
Python Exam BIF-30806
"""

import sys
import subprocess
import os.path

def write_fasta(record, output_name, mode = 'w', line_width = 80):
    """Write a FASTA file with option to append to exisitng file.

    Key inputs:
    record --- the dictionary the DNA sequence with sequence IDs as keys.
    output_name --- the name of the finished FASTA file
    mode --- the mode for writing to the FASTA file; 'w' writes a new file
    (overwriting if necessary) and 'a' allows appending to a pre-existing file.
    line_width --- the number of characters allowed on each line of the FASTA
    file.
    """
    assert mode == 'w' or mode == 'a', 'Mode must be \'w\' or \'a\'. '
    with open(output_name, mode) as f:
        for key in record.keys():
            title = key + '\n'
            f.write(title)
            range_sup = ((len(record[key]) // line_width) + 1) * line_width
            for i in range(0, range_sup + 1, line_width):
                f.write(record[key].upper()[i:min(i + line_width, len(record[key]))])
                f.write('\n')
            f.write('\n')
    return

def split_records_aug(filename):
    """Yields a generator of the different generated gene contained in filename.

    Key inputs:
    filename --- the file containing AUGUSTUS generated genes.
    """
    current = {}
    started = False
    count = 0
    with open(filename) as f:
        for line in f:
            if not line.strip():
                continue
            if 'start' in line:
                started = True
                continue
            if started:
                if not line.startswith('#'):
                    if count == 0:
                        line = line.replace('\t', ' ')
                        header = '>' + line.strip()
                        current[header] = ''
                    count += 1
                elif ']' in line:
                    current[header] += line.partition('# ')[2].partition(']')[0]
                    count = 0
                    started = False
                    yield current
                    current = {}
                else:
                    line.strip()
                    if '[' in line:
                        current[header] += line.strip().partition('[')[2]
                    elif line.startswith('# '):
                        current[header] += line.partition('# ')[2].strip()
    if current:
        yield current

def blast_database(target, dbtype, output_to_file = False, output_file = None,
          overwrite = False):
    """Creates a blast database using makeblastdb tool.

    Key inputs:
    target --- a FASTA file.
    dbtype --- molecule type; 'nucl' or 'prot'.
    output_to_file --- boolean defining output as file or variable
    output_file --- only used if output_to_file is True; filename of 
    desired output file
    overwrite --- boolean for overwriting pre-exisitng output file of the
    same name
    """
    if output_to_file:
        if os.path.exists(output_file) and not overwrite:
            return output_file
        cmd = 'makeblastdb -in {} -dbtype {} -out {}'.format(target, dbtype, output_file)
    else:
        cmd = 'makeblastdb -in {} -dbtype {}'.format(target, dbtype)
    printed_output = subprocess.check_output(cmd, shell=True)

    if output_to_file:
        return output_file

    return printed_output

def blastp(database, query, output_to_file = False, output_file = None,
          overwrite = False, outfmt = 7):
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
        if os.path.exists(output_file) and not overwrite:
            return output_file
        cmd = 'blastp -db {} -query {} -outfmt {} -out {} -num_alignments 1'.\
               format(database, query, outfmt, output_file)
    else:
        cmd = 'blastp -db {} -query {} -outfmt {} -num_alignments 1'.format(
            database, query, outfmt)

    printed_output = subprocess.check_output(cmd, shell=True)
    if output_to_file:
        return output_file
    return printed_output

def call_files():
    """Reads in filenames and command to write to file for main program."""
    try:
        predicted_proteins = sys.argv[1]
    except IndexError:
        predicted_proteins = input('Please input AUGUSTUS file for analysis: ')
    try:
        protein_db = sys.argv[2]
    except IndexError:
        protein_db = input('Please input a protein database file: ')

    try:
        output_file_aug_to_fasta = sys.argv[3]
        output_to_file = True
    except IndexError:
        output_to_file = input('Write output to file?'
                               + ' [Yes/No]: ')
        if output_to_file.upper() in 'YES':
            output_to_file = True
            output_file_aug_to_fasta = input('Please supply output file name '
                                             + 'for AUGUSTUS conversion to '
                                             + 'FASTA: ')
        else:
            output_to_file = False
            output_file_aug_to_fasta = None

    try:
        output_file_proteins_to_db = sys.argv[4]
    except IndexError:
        if output_to_file:
            output_file_proteins_to_db = input('Please supply output file name'
                                               + 'for blast database: ')
        else:
            output_file_proteins_to_db = None

    try:
        blastp_output = sys.argv[5]
    except IndexError:
        if output_to_file:
            blastp_output = input('Please supply output file name for blastp: ')
        else:
            blastp_output = None

    finally:
        if len(sys.argv) >= 7:
            overwrite = sys.argv[6]
        elif output_file and os.path.exists(output_file):
            overwrite = input('Output file already exists. Overwrite? '
                              + '[Yes/No]: ')
            if overwrite.upper() in 'YES':
                overwrite = True
            else:
                overwrite = False
        else: overwrite = False

    return (predicted_proteins, protein_db, output_file_aug_to_fasta, 
            output_file_proteins_to_db, blastp_output, 
            output_to_file, overwrite)

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

def query_length(query):
    """Returns the length of a query from blastp search.

    Key inputs:
    query --- a dict based on output of single alignment from blastp search
    """
    return int(query['q. end']) - int(query['q. start']) + 1

def query_coverage(query):
    """Returns coverage of alignment

    Key inputs:
    query --- dict from single alignment of blastp output.
    """
    length = query_length(query)
    coverage = (int(query['alignment length'])/length) * 100
    return coverage

def print_output(blast_results):
    """Print desired output of BLASTP in tab delimited format.

    Key inputs:
    blast_results --- list of dicts of blastp output
    """
    output = ''
    header = ['query ID', 'query length', 'subject ID', 'percent identity',
              'query coverage']
    line = '\t'.join(header) + '\n'
    output += line
    for query in blast_results:
        length = str(query_length(query))
        coverage = query_coverage(query)
        coverage = '{:.2}'.format(coverage)
        identity = '{:.2f}'.format(float(query['% identity']))
        line = '\t'.join([query['query acc.ver'], length, 
                         query['subject acc.ver'], identity,
                         coverage])
        line += '\n'
        output += line
    print(output)

def main():
    """BLASTP search using AUGUSTUS generated sequences on FASTA file.
    """
    count = 0

    # Read in the required files and filenames.
    predicted_proteins, protein_db, output_file_aug_to_fasta, \
        output_file_proteins_to_db, blastp_output, output_to_file, \
        overwrite = call_files()

    # Write all entries in the AUGUSTUS output to a FASTA file
    for record in split_records_aug(predicted_proteins):
        if count == 0:
            mode = 'w'
        else:
            mode = 'a'
        write_fasta(record, output_file_aug_to_fasta, mode)
        count += 1

    # Create a blast database and carry out a blastp search
    blast_db = blast_database(protein_db, 'prot', True,
                              output_file_proteins_to_db, overwrite)

    blastp_file = blastp(output_file_proteins_to_db, output_file_aug_to_fasta,
                         True, blastp_output, overwrite, 7)

    # Parse the blastp results for the desired information
    blast_results = parse_blastp_output(blastp_output)

    # Print the results
    print_output(blast_results)

if __name__ == '__main__':
    main()