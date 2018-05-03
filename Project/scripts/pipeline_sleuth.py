#!/usr/bin/env python3
"""
Input csv file describing samples to use in pipeline.

Find files, unzip and unmerge as required. Apply fastq trimming,
Salmon indexing and quantifying of RNA-seq reads, wasabi and sleuth.

Author: Stephen Coleman
Student id: 940-309-160-050

Sample Command:
python3 /local/data/course/project/groups/FINAL/scripts/pipeline_sleuth.py <- this is the python file
/local/data/course/project/groups/sewer/transcriptome_data_tissue_mk2.csv <- this is the csv used as metadata for sleuth and instruction on which files to use in the pipeline
/local/data/course/project/transcriptomes/Duge_de_Bernonville_et_al/CDF97_v2.fa <- the transcriptome
 -b 50 <- the number of bootstraps for salmon
-s False <- skip failures (now defunct)
-o False <- overwrite files if found to already exist
"""

from subprocess import check_call
from os import getcwd
from os.path import exists
from argparse import ArgumentParser
from pathlib import Path
import tempfile
import csv

current_working_directory = getcwd()

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

def check_paired_end_reads(entry):
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
    AP = ArgumentParser(description='This script implements a pipeline of fastq_quality_trim, '
                                    + 'salmon and wasabi, using a csv file of samples to use and a '
                                    + 'given transcriptome for indexing')
    AP.add_argument('filename', help='csv file of RNA-seq libraries to apply pipeline to')
    AP.add_argument('transcriptome', help='used for indexing the reads')
    AP.add_argument('-b', '--bootstrap', help='number of bootstraps, def=0', \
                    default=0) #optional
    AP.add_argument('-o', '--overwrite', help='overwrite pre-existing output, def=No [yes/no]', \
                    default='No') #optional
    AP.add_argument('-s', '--skip', help='skip files that fail, def=False', default='no')
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

def remove_whitespace(filename):
   """Remove excess whitespace from file. Now surplus to requirements."""
   with open(filename) as fin, open('{}_temp', 'w') as fout:
       for line in fin:
           if line.strip():
               fout.write(line)
   return '{}_temp'.format(filename)

def split_merged_pairs(filename):
    """Splits merged pair end reads"""
    cmd = 'python3 /local/data/course/project/groups/sewer/FINAL/scripts/' \
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

def check_complete(filename_without_extension):
    """Checks for the presence of a 'complete' version of the file (see SRR1271857 for why)

    Key inputs:
    filename_without_extension --- string filename
    """
    if Path('{}_complete.fastq'.format(filename_without_extension)).is_file():
        return '{}_complete'.format(filename_without_extension), True
    else:
        return filename_without_extension, False

def check_all_paths(filename, directory, paired=False):
    """Returns the file name and location if the file is present in the local 
    directory or the given directory, with preference given to the local directory.

    Key inputs:
    filename --- string, name of file to be found lacking the file extension
    directory --- the directory to be searched after the local directory
    paired --- boolean indicating if the reads are paired or not. Default is False
    """
    file_used = []
    pair_found = False
    location = ''
    if paired:
        for i in range(1,3):
            if Path('./{}_complete_{}.fastq'.format(filename, i)).is_file():
                file_used += ['./{}_complete_{}'.format(filename, i)]
            elif Path('./{}_{}.fastq'.format(filename, i)).is_file():
                file_used += ['./{}_{}'.format(filename, i)]
            elif Path('{}/{}_complete_{}'.format(directory, filename, i)).is_file():
                file_used += ['{}/{}_complete_{}'.format(directory, filename, i)]
            elif Path('{}/{}_{}.fastq'.format(directory, filename, i)).is_file():
                file_used += ['{}/{}_{}'.format(directory, filename, i)]
            if len(file_used) == 2:
                pair_found = True
    if not file_used:
        if Path('./{}_complete.fastq'.format(filename)).is_file():
            file_used = ['./{}_complete'.format(filename)]
        elif Path('./{}.fastq'.format(filename)).is_file():
            file_used = ['./{}'.format(filename)]
        elif Path('{}/{}_complete'.format(directory, filename)).is_file():
            file_used = ['{}/{}_complete'.format(directory, filename)]
        elif Path('{}/{}.fastq'.format(directory, filename)).is_file():
            file_used = ['{}/{}'.format(directory, filename)]
    return file_used, pair_found

def check_all_paths_including_complete(filename, directory, paired=False):
    """Returns the file name and location if the file is present in the local
    directory or the given directory, with preference given to the local directory
    and files appended by '_complete'.

    Key inputs:
    filename --- string, name of file to be found lacking the file extension
    directory --- the directory to be searched after the local directory
    paired --- boolean indicating if the reads are paired or not. Default is False
    """
    file_used = []
    pair_found = False
    location = ''
    if paired:
        for i in range(1,3):
            if Path('./{}_complete_{}.fastq'.format(filename, i)).is_file():
                file_used += ['./{}_complete_{}'.format(filename, i)]
            elif Path('./{}_{}.fastq'.format(filename, i)).is_file():
                file_used += ['./{}_{}'.format(filename, i)]
            elif Path('{}/{}_complete_{}'.format(directory, filename, i)).is_file():
                file_used += ['{}/{}_complete_{}'.format(directory, filename, i)]
            elif Path('{}/{}_{}.fastq'.format(directory, filename, i)).is_file():
                file_used += ['{}/{}_{}'.format(directory, filename, i)]
            if len(file_used) == 2:
                pair_found = True
    if not file_used:
        if Path('./{}_complete.fastq'.format(filename)).is_file():
            file_used = ['./{}_complete'.format(filename)]
        elif Path('./{}.fastq'.format(filename)).is_file():
            file_used = ['./{}'.format(filename)]
        elif Path('{}/{}_complete'.format(directory, filename)).is_file():
            file_used = ['{}/{}_complete'.format(directory, filename)]
        elif Path('{}/{}.fastq'.format(directory, filename)).is_file():
            file_used = ['{}/{}'.format(directory, filename)]
    return file_used, pair_found

def call_fastq_trimmer(dict_file, args):
    """Check if files are paired, unmerge if necessary and trim.

    Key inputs:
    dict_file --- dictionary describing RNA-seq file
    args --- command line input for program
    """

    # Define some initial parameters and descriptions of the file and its location
    paired_end = check_paired_end_reads(dict_file)
    directory = dict_file['SRP'].strip()
    filename = [dict_file['run accession'].strip(),]
    py_file = '/local/data/course/project/groups/sewer/FINAL/scripts/' \
               + 'MakeFASTQCReport.py'
    min_read_length = 30 # float(dict_file['Read length']) * 0.9
    library_dir = '/local/data/course/project/RNAseq/{}'.format(directory)
    current_file = '{}/{}.fastq'.format(library_dir, filename[0])

    # Initialise output
    output = []

    # Find file in local directory or given directory (handles 'complete')
    file_path, pair_found = check_all_paths_including_complete(filename[0],
                                                               library_dir,
                                                               paired_end
                                                              )

    # Unzip file if required /requested
    if Path('{}.gz'.format(current_file)).is_file() and (not file_path
                                                         or args.overwrite
                                                        ):
            file_loc = current_file.rsplit('/', 1)[0]
            unzipped = unzip_fastq_file('{}.fastq'.format(filename[0]), file_loc)
            library_dir = current_working_directory
            file_path, pair_found = check_all_paths_including_complete(filename[0],
                                                                       library_dir,
                                                                       paired_end)

    # Split merged pair reads if required / requested
    if paired_end and (not pair_found or args.overwrite):
        file_to_unmerge, unused_flag = check_all_paths_including_complete(filename[0],
                                                                          library_dir,
                                                                          False)
        if file_to_unmerge:
            split_merged = split_merged_pairs('{}.fastq'.format(file_to_unmerge[0]))
            file_path, pair_found = check_all_paths(filename[0], library_dir, paired_end)

    for index, file in enumerate(file_path):
        file_name = file.rsplit('/', 1)[-1]
        if not args.overwrite:
            try:
                f = open('{}/{}_Qtrimmed.fastq'.format(current_working_directory, file_name))
            #if Path('{}/{}_Qtrimmed.fastq'.format(current_working_directory, file_name )).is_file:
                print('\n{}/{}_Qtrimmed.fastq already exists. Moving to next file.\n'
                      .format(current_working_directory, file_name))
                output += ['{}_Qtrimmed.fastq'.format(file_name)]

                # Check if both files checked
                if index == (len(file_path) - 1):
                    return output
                else:
                    continue
            except:
                pass
        else:
            pass

        output += ['{}_Qtrimmed.fastq'.format(file_name)]
        cmd = 'python3 {} {}.fastq {} {}'.format(py_file, file, min_read_length, 10)

        if args.skip:
            try:
                res = check_call(cmd, shell=True)
            except:
                return [directory, file_name]
        else:
            res = check_call(cmd, shell=True)
    return output

def call_salmon(dict_file, args, fastq_trim_output_name):
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

    paired_end = check_paired_end_reads(dict_file)
    paired = '-m '
    bootstrap = args.bootstrap
    transcriptome = args.transcriptome
    directory = dict_file['SRP']

    filename = dict_file['run accession']
    if not paired_end:
        paired = '-r '

    if paired_end:
        library_dir = '{} {}'.format(fastq_trim_output_name[0],
                                     fastq_trim_output_name[1])
    else:
        library_dir = '{}'.format(fastq_trim_output_name[0])
    py_file = '/local/data/course/project/groups/sewer/FINAL/scripts/salmon.py'
    cmd = 'python3 {} {} '.format(py_file, transcriptome) \
          + '{}{} -b {}'.format(paired, library_dir, bootstrap)

    if args.skip:
        try:
            res = check_call(cmd, shell=True)
        except:
            return [directory, filename]
    else:
        res = check_call(cmd, shell=True)
    return None

def call_wasabi():
    R_file = '/local/data/course/project/groups/sewer/FINAL/scripts/wasabi.R'
    cmd = 'Rscript {} .'.format(R_file)
    res = check_call(cmd, shell=True)
    return None

def call_sleuth(args):
    """Call sleuth.R as a subprocess on the salmon output"""
    directory_of_salmon = current_working_directory
    R_file = '/local/data/course/project/groups/sewer/FINAL/scripts/sleuth.R'
    cmd = 'Rscript {} {} {}.'.format(R_file, directory_of_salmon, args.filename)
    res = check_call(cmd, shell=True)
    return None

def main():
    """Main program; request inputs and applies fastq trimming, salmon,
    wasabi and sleuth to the samples listed in the input csv file
    """
    # Read in arguments
    args = argument_details()

    # Convert the csv file describing the samples to analyse to a dict
    dict_of_libraries = read_csv_list_to_dict(args.filename)

    # Initialise lists of failures
    failures_fastq_trimming = []
    failures_salmon = []
    outright_failures = []

    # Cycle through each sample contained in input csv
    # Find the fastq files, unzip and unmerge as required
    # Trim and apply salmon
    for entry in dict_of_libraries:
        if args.skip:
            try:
                fastq_output = call_fastq_trimmer(entry, args)
                failures_salmon += [call_salmon(entry, args, fastq_output)]
            except:
                outright_failures += [entry]
        else:
            fastq_output = call_fastq_trimmer(entry, args)
            failures_salmon += [call_salmon(entry, args, fastq_output)]

    print('Failures for FASTQ trimming: ')
    for failure in failures_fastq_trimming:
        if failure != None:
            print(failure)
    print('Failures for salmon: ')
    for failure in failures_salmon:
        if failure != None:
            print(failure)

    call_wasabi()
    call_sleuth(args)
    return None

if __name__ == "__main__":
    main()
