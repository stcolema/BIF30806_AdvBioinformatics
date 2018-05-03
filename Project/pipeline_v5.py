#!/usr/bin/env python3
"""
Salmon indexing and quantifying of RNA-seq reads

Author: Thijs Maas
Student nr: 921116537120
"""

# Stephen copy
# Plan to read in csv file of files to use as Salmon input

from subprocess import check_call
from os.path import exists
from argparse import ArgumentParser
from pathlib import Path
import tempfile
import csv

def read_csv_list_to_dict(filename):
    """Given a csv with columns 'SRP', 'run accession', and 'library layout'
    (as a minimum) describing RNA-seq libraries, returns a dict of said descriptions.

    Key inputs:
    filename --- name of csv file
    """
    count = 0
    result = []
    with open(filename) as fin:
        reader = csv.reader(fin)
        for line in reader:
            if count == 0:
                line[0] = 'SRP'
                key = [item.strip() for item in line]
            else:
                entry = dict(zip(key, line))
                result += [entry]
            count += 1
    print(key)
    return result

def check_paried_end_reads(entry):
    """For dict describing RNA-seq library

    Key inputs:
    entry --- dict describing RNA-seq library
    """
    if 'PAIRED' in entry['library layout']:
        paired_end_reads = True
    else:
        paired_end_reads = False
    return paired_end_reads

def argument_details():
    # argument parsing
    AP = ArgumentParser(description='This script is used to first index, and \
    then quantify a csv file of RNA-seq reads library by pseudomapping using Salmon.')
    AP.add_argument('filename', help='file of RNA-seq libraries to apply salmon to')
    AP.add_argument('transcriptome', help='used for indexing the reads')
    AP.add_argument('-b', '--bootstrap', help='number of bootstraps, def=0', \
                    default=0) #optional
    AP.add_argument('-o', '--overwrite', help='overwrite pre-existing output, def=No [yes/no]', \
                    default='No') #optional
    AP.add_argument('-s', '--skip', help='skip files that fail, def=False', default=False)
    args = AP.parse_args()
    if args.overwrite in 'Yesyes' or args.overwrite in 'True':
        args.overwrite = True
    else:
        args.overwrite = False
    if args.skip in 'Yesyes' or args.skip in 'True':
        args.skip = True
    else:
        args.skip = False
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

def split_merged_pairs(filename):
    """Splits merged pair end reads"""
    cmd = 'python3 /local/data/course/project/groups/sewer/splittedfiles/' \
          + 'pairedsplitter.py {}'.format(filename)
    res = check_call(cmd, shell=True)
    return True

def unzip_fastq_file(unzipped_filename, file_location):
    """Unzips file and saves locally.

    Key inputs:
    unzipped_filename --- string of filenames to be unzipped without gz extension
    file_location --- directory name containing file to be unzipped
    """
    cmd = 'gunzip {1}/{0}.gz -c > {0}'.format(unzipped_filename, file_location)
    print("First gunzip this file...")
    check_call(cmd, shell=True)
    return True

def call_fastq_trimmer(dict_file, args):
    """Check if files are paired, unmerge if necessary and trim.

    Key inputs:
    dict_file --- dictionary describing RNA-seq file
    args --- command line input for program
    """
    split_merged = False
    paired_end = check_paried_end_reads(dict_file)
    directory = dict_file['SRP'].strip()
    filename = [dict_file['run accession'].strip(),]
    py_file = '/local/data/course/project/groups/sewer/FASTQCStuff/MakeFASTQCReport.py'
    min_read_length = 30 # float(dict_file['Read length']) * 0.9
    library_dir = '/local/data/course/project/RNAseq/{}'.format(directory)
    unzipped = False
    current_file = '{}/{}.fastq'.format(library_dir, filename[0])
    if Path('{}.gz'.format(current_file)).is_file():
        unzipped = unzip_fastq_file('{}.fastq'.format(filename[0]), library_dir)
        library_dir = '.'

    if paired_end:
        filename = ['{}_1'.format(filename[0]), '{}_2'.format(filename[0])]
        if not Path('/{}.fastq'.format(library_dir, filename[0])).is_file():
            if not Path('./{}.fastq'.format(filename[0])).is_file():
                if Path('/local/data/course/project/RNAseq/{}/{}.fastq'.format(directory, dict_file['run accession'].strip())).is_file():
                    split_merged = split_merged_pairs('/local/data/course/project/RNAseq/{}/{}.fastq'.format(directory, dict_file['run accession'].strip()))
                else:
                    for file in filename:
                        split_merged = split_merged_pairs('/local/data/course/project/RNAseq/{}/{}.fastq'.format(directory, file))

    if not args.overwrite:
        if Path('{}_Qtrimmed.fastq'.format(filename[0])).is_file():
            print('\n{}_Qtrimmed.fastq already exists. Moving to next file.\n'
                  .format(dict_file['run accession']))
            return None

    for file in filename:
        if split_merged:
            library_dir = '.'

        print(file)

        input_file = '{}/{}.fastq'.format(library_dir, file)
        output = '{}_Qtrimmed.fastq'.format(file)

        cmd = 'python3 {} {} {} {}'.format(py_file, input_file, min_read_length, 10)

        if args.skip:
            try:
                res = check_call(cmd, shell=True)
            except:
                return [directory, filename]
        else:
            res = check_call(cmd, shell=True)
    return None

def call_salmon(dict_file, args):
    """Call salmon.py on dict describing the RNA-seq library file.

    Key inputs:
    dict_file --- entry from output of read_csv_to_dict
    args --- arguments requested in argument_details function

    Returns None
    """
    if not args.overwrite and Path('{}trim_bs{}_quant'.format(dict_file['run accession'], 
                                                          args.bootstrap)
                                   ).is_dir():
        print('\n{}_bs{}_quant already exists. Moving to next file.\n'
              .format(dict_file['run accession'], args.bootstrap))
        return None
    paired_end = check_paried_end_reads(dict_file)
    paired = '-m '
    bootstrap = args.bootstrap
    transcriptome = args.transcriptome
    directory = dict_file['SRP']
    filename = dict_file['run accession']
    if not paired_end:
        paired = '-r '

#    library_dir = '/local/data/course/project/RNAseq/{}/{}' \
#    library_dir = '/local/data/course/project/groups/sewer/FASTQCtuff/{}' \
    library_dir = './{}' \
                  .format(filename)
    if paired_end:
        library_dir = '{}_1_Qtrimmed.fastq {}_2_Qtrimmed.fastq'.format(library_dir, library_dir)
    else:
        library_dir = '{}_Qtrimmed.fastq'.format(library_dir)
    py_file = '/local/data/course/project/groups/sewer/salmon/salmon.py'
    cmd = 'python3 {} {} '.format(py_file, transcriptome) \
          + '{}{} -b {}'.format(paired, library_dir, bootstrap)
    print(cmd)
    if args.skip:
        print(args.skip)
#        print('\n\n Stephen \n\n')
        try:
            res = check_call(cmd, shell=True)
        except:
            return [directory, filename]
    else:
        res = check_call(cmd, shell=True)
    return None

def call_wasabi():
    R_file = '/local/data/course/project/groups/sewer/wasabi/wasabi_no_warnings.R'
    cmd = 'Rscript {} .'.format(R_file)
    res = check_call(cmd, shell=True)
    return None

def main():
    args = argument_details()
    dict_of_libraries = read_csv_list_to_dict(args.filename)
    failures_fastq_trimming = []
    failures_salmon = []
    for entry in dict_of_libraries:
        failures_fastq_trimming += [call_fastq_trimmer(entry, args)]
        failures_salmon += [call_salmon(entry, args)]
    print('Failures for FASTQ trimming: ')
    for failure in failures_fastq_trimming:
        if failure != None:
            print(failure)
    print('Failures for salmon: ')
    for failure in failures_salmon:
        if failure != None:
            print(failure)
    call_wasabi()
    return None

if __name__ == "__main__":
    main()
