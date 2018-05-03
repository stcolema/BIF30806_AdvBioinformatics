#!/usr/bin/env python3
"""
Script to split a merged paired-end file into two separate files

Author: Thijs Maas
Student nr: 921116537120
"""
from argparse import ArgumentParser
import sys
from datetime import datetime as dt

def splitter(line):
    """Returns the left and right readline seprerated from a merged line"""
    if line.startswith('@') and 'length=' in line:
        newlen = int(line.split('length=')[-1]) / 2
        head1 = line.split('length=')[0] + '1 length='+str(newlen)
        head2 = line.split('length=')[0] + '2 length='+str(newlen)
        return (head1, head2)
    elif line.startswith('+') and (len(line) == 1 or 'lenght' in line):
        return ('+', '+')
    else:
        if line.strip():
            line1 = line[0:len(line)//2]
            line2 = line[len(line)//2:]
            return (line1, line2)

def progress(c, start, total, lpr):
    """Return % done and ETA"""
    if (c+1) % (lpr*250) == 0:
        done = (c+1) / lpr
        now = dt.now()
        elaps = now - start
        percentage = done/total
        time_left = (elaps/done*total) - elaps
        sys.stdout.write("\r Split {} reads of {:0.0f} -> {:.2%} - Time left {}\
        ".format(done, total, percentage, str(time_left).split('.')[0]))
    else:
        return None

if __name__ == "__main__":

    # argument parsing
    AP = ArgumentParser(description='This script takes an fastq or fasta merged\
    paired end file and splits it into two files')
    AP.add_argument('pairedfile', help='merged paired-end FASTQ or FASTA file')
    args = AP.parse_args()

    print("Reading file...")

    with open(args.pairedfile, 'r') as f:
        # get file type
        if f.readlines()[2].startswith('+'):
            filetype = 'fastq'
            lpr = 4 # lines per read
        elif f.readlines()[2].startswith('@'):
            filetype = 'fasta'
            lpr = 2
        else:
            raise ReadError
        f.seek(0)

        # get file length
        total = len(f.readlines()) / lpr
        print("The file has {:,} reads".format(total))
        f.seek(0)

        # make names from file _1 _2
        name = args.pairedfile.split('/')[-1]
        name_1 = name.split('.fastq')[0].split('.fasta')[0]+'_1.'+filetype
        name_2 = name.split('.fastq')[0].split('.fasta')[0]+'_2.'+filetype
        print("File {} is going to be split into: \n{} \n{}".format(\
            name, name_1, name_2))

        # split file and write to new files
        with open(name_1, 'w') as w1, open(name_2, 'w') as w2:
            startt = dt.now()
            for c, line in enumerate(f.readlines()):
                print(splitter(line)[0].strip(), file=w1)
                print(splitter(line)[1].strip(), file=w2)

                # write progress
                progress(c, startt, total, lpr)
    # finish up
    sys.stdout.flush()
    sys.stdout.write("\r Split {:0.0f} reads of {:0.0f} -> {:.2%}\
                    ".format(total, total, 1))
    sys.stdout.flush()
    endt = dt.now()
    elaps = endt - startt
    print()
    print("100% done! {:0.0f} reads were splitted in {}.".format(total, \
    str(elaps).split('.')[0]))
