#!/user/bin/env python3
"""
Creates a metadata table from an input table and whatever directories are actually available
(based on an earlier stage of the sleuth.R program which used all the '_quant' directories present
in the given directory. This reduces the metadata table down to the if it has more entries that
are not present in the directory.

Could be useful as part of a pipeline which allows the user to decide if it should continue
regardless of some files failing steps.

Author: Stephen Coleman
Student ID: 940-309-160-050
"""


from argparse import ArgumentParser
from subprocess import check_call

def argument_details():
    # argument parsing
    AP = ArgumentParser(description='This script is used to create a metadata table for sleuth analysis.')
    AP.add_argument('directory', help='directory to search for quant directories within')
    AP.add_argument('metadata', help='csv file of metadata to use')
    AP.add_argument('quant_output', help='filename to save list of quant directories to')
    AP.add_argument('metadata_output', help='filename to save relevant metadata to')
    args = AP.parse_args()
    return args

def find_all_quant_dir(args):
    """Find all directories containing the word quant"""
    cmd = 'find {} -maxdepth 1 -name "*quant" > {}'.format(args.directory, args.quant_output)
    check_call(cmd, shell=True)

def create_metadata_table(args):
    """Writes a csv file of feasible samples for sleuth based on those requested in 
    args.metadata"""
    quant_dirs = open(args.quant_output).readlines()
    quant_dirs = '\t'.join(quant_dirs)
    with open(args.metadata) as metadata_f, \
         open(args.metadata_output, 'w') as metadata_fout:
        for index, line in enumerate(metadata_f):
            line_loc = line.strip().split(',')
            if line_loc[2] in quant_dirs or index == 0:
                metadata_fout.write(line)
            else:
                continue

def main():

    args = argument_details()
    find_all_quant_dir(args)
    create_metadata_table(args)

if __name__ == '__main__':
    main()
