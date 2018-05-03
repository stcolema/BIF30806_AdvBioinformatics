#!/usr/bin/env python3

"""
Author: Stephen Coleman
Student number: 940309-160-050
Python Mock Exam BIF-30806

Receive an assembly of yeast chromosome 3 (made using Velvet)
Compare it to the sequence of chromosome 3 in the reference genome using LASTZ
Report the regions in the reference genome that are not covered by the Velvet
assembly. 
Compute and output some statistics on both the reference and Velvet assemblies.

Sample report on uncovered regions:
0:10 ACGTGGACAT
23:25 GT
45:46 C

Assignment
Write a python script that performs the following tasks:
1. report the assembly size, N50 size, and N50 index for each assembly
2. compare the two assemblies using lastz
3. report the regions from the reference genome that are not covered by the Velvet
assembly

N50 is contig length of smallest contig required in sum of contigs (in desc.
order) to cover 50% of the assembly.
Aaembly size is the length of the assembly?
N50 index?

"""


import sys
import subprocess
import os.path
import re

def fasta_reader(filename, output = dict):
    """Imports a FASTA file as a dict or list depending on user selection.
    
    Key inputs:
    filename --- the name of the FASTA file to be converted into a dict
    output --- either dict or list; defines the format of the FASTA file after
    parsing.
    """
    assert output in [dict, list], 'Output must be dict or list. '

    if output == dict:
        result = {}
    else:
        result = []
        current = -1
    # with open(filename) as f:
    # for line in f:
    for line in filename:    
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            
            if output == dict:
                current = line[1:]
                result[current] = {}
            else:
                current += 1
                result += [{}]
                result[current]['NAME'] = line[1:]
            result[current]['SEQUENCE'] = ''
        else:
            result[current]['SEQUENCE'] += line
    return result

def sequence_length(fasta):
    """Adds a entry {'SEQ LENGTH': length of sequence} to a FASTA dict / list.
    
    Key inputs:
       fasta --- a dictionary or list output of the function fasta_reader as 
    performed on a FASTA file
    """
    input_type = type(fasta)
    assert input_type in [dict, list], 'input_type must be list or dict. '

    if input_type == dict:
        iterate_list = fasta.keys()
    else:
        iterate_list = range(len(fasta))
    for seq in iterate_list:
        fasta[seq]['SEQ LENGTH'] = len(fasta[seq]['SEQUENCE'])
    return fasta

def sort_list_on_key(lst, sort_key, descending = True):
    """Returns list of dicts sorted based on sort_key.

    Key inputs:
       lst --- list of dicts
       sort_key --- key common to all entries in lst on which to sort lst
       descending --- boolean for descending or ascending order
    """
    assert type(lst) == list, 'lst input must be of a list of dicts. '
    try:
        return sorted(lst, key=lambda elem: elem[sort_key], reverse=descending)
    except KeyError as err:
        raise KeyError('Key not in elements in list. ') from err
    except Exception as err:
        raise Exception(err) from err
    
def assembly_length(assembly):
    """Returns the combined sequence length of a list of contigs.

    Key inputs:
       assembly --- list of dicts with key 'SEQ LENGTH'.
    """
    assert type(assembly) == list, 'assembly must be a list of dicts. '
    try:
        return sum(assembly[i]['SEQ LENGTH'] for i in range(len(assembly)))
    except KeyError as err:
        raise KeyError('Key not in dict; please ensure all entries of list'\
                        + ' include a \'SEQ LENGTH\' value. ')
    except Exception as err:
        raise Exception(err)

def n50(assembly):
    """Returns a list of contigs ordered on sequence length, n50 index and n50

    Key inputs:
       assembly --- a list of dicts containing at least a key 'SEQ LENGTH'.
    """
    new_list = sort_list_on_key(assembly, 'SEQ LENGTH', True)
    n50_len = assembly_length(assembly) / 2
    count = 0
    for i in range(len(new_list)):
        count += new_list[i]['SEQ LENGTH']
        if count > n50_len:
            return new_list, i, new_list[i]['SEQ LENGTH']

    raise Exception('N50 not found. Please investigate.')

def lastz(target, query, output_to_file = False, output_file = None,
          overwrite = False):
    """Run lastz program on two fasta files, returning lastz output.

    Key inputs:
       target --- a FASTA file.
       query --- a FASTA file of Velvet output for mapping to sequence.
       output_to_file --- boolean defining output as file or variable
       output_file --- only used if output_to_file is True; filename of 
    desired output file
       overwrite --- boolean for overwriting pre-exisitng output file of the
    same name
    """
    if output_to_file:
        if os.path.exists(output_file) and not overwrite:
            return output_file
        cmd = 'lastz {} {} --output={} --format=general'.format(target, query, output_file)
    else:
        cmd = 'lastz {} {} --format=rdotplot'.format(target, query)
    
    # Return output of lastz.
    printed_output = subprocess.check_output(cmd, shell=True)
    if output_to_file:
        return output_file
    return printed_output

def number_fasta_entries(filename):
    """Returns number of unique entries in FASTA file.

    Key inputs:
       filename --- name of input file
    """
    file_contents = open(filename).read() # file type always small
    contents = set(file_contents.split('>'))
    return len(contents)

def call_files():
    """Reads in files and inputs for LASTZ and main program."""

    try:
        target_file = sys.argv[1]
    except IndexError:
        target_file = input('Please input target file for analysis: ')
    try:
        query_file = sys.argv[2]
    except IndexError:
        query_file = input('Please input a query for related reads: ')

    try:
        output_file = sys.argv[3]
        output_to_file = True
    except IndexError:
        output_to_file = input('Write output to file?'
                               + ' [Yes/No]: ')
        if output_to_file.upper() in 'YES':
            output_to_file = True
            output_file = input('Please input output file name for LASTZ '
                                + 'analysis: ')
        else:
            output_to_file = False
            output_file = None
    finally:
        if len(sys.argv) >= 5:
            overwrite = sys.argv[4]
        elif output_file and os.path.exists(output_file):
            overwrite = input('Output file already exists. Overwrite? '
                              + '[Yes/No]: ')
            if overwrite.upper() in 'YES':
                overwrite = True
            else:
                overwrite = False
        else: overwrite = False

    return target_file, query_file, output_file, output_to_file, overwrite

def covered_region(lastz_output):
    """Returns list of tuples of covered regions from LASTZ output.

    Key inputs:
    lastz_output --- list of lines of LASTZ output.
    """
    covered_region = []
    for i, line in enumerate(lastz_output):
        if line.startswith('#'):
            continue
        else:
            loc_covered_region = line.strip().split()
            loc_covered_region = loc_covered_region[4:6]
            loc_covered_region = tuple(map(int, loc_covered_region))
            covered_region += [loc_covered_region]

    # Sort list based on start positions (in ascending order) 
    # and subsort based on end positions (in descending order).
    # Thus if positions overlap, those which cover more of the 
    # target query will come first - this is useful in the next
    # step of analysis.
    covered_region.sort(key = lambda elem: (elem[0], -elem[1]), reverse=False)
    return covered_region


def uncovered_region(sequence_length, covered_region):
    """Returns a list of tuples (start, end) of the uncovered regions.

    Key inputs:
    sequence_length --- the target sequence length
    covered_region --- list of tuples of regions covered by contigs.
    """
    uncovered_region = []
    last_end = 0
    for index, start, end in enumerate(covered_region):
        if index >= 1:
            if start > last_end:
                uncovered_region += [(last_end, start)]
                last_end = end
            elif start == last_end:
                last_end = end
    if covered_region[-1][1] != sequence_length:
        uncovered_region += [(covered_region[-1][1], sequence_length)]
    return uncovered_region

def length_sequences_fasta_file(filename):
    """Returns the length of sequences in a fasta file as a list on numbers.

    Key inputs:
    filename --- the FASTA file name.
    """
    sequence_length = []
    count = -1
    with open(filename) as f:
        for line in f:
            line.strip()
            if not line:
                continue
            if line.startswith('>'):
                sequence_length += [0]
                count += 1
            else:
                sequence_length[count] += len(line)
    return sequence_length

def uncovered_sequences(filename, uncovered_region):
    """Returns tuple of (start, end, sequence) of uncovered regions

    Key inputs:
    filename --- name of target fasta file
    uncovered_region --- list of tuples of uncovered regions in sequence
    """
    current_pos = 0
    region = 0
    next_line = False
    line_end = False
    list_uncovered_sequences = []
    with open(filename) as f:
        for line in f:
            while line_end:
                if uncovered_region[region][0] in range(current_pos, current_pos+len(line)):
                    if uncovered_region[region][1] in range(current_pos, current_pos+len(line)):
                        uncovered_sequence = line[uncovered_region[region][0]-current_pos : uncovered_region[region][1]-current_pos+1]
                        region += [uncovered_sequence]
                        line_end = True
                    else:
                        uncovered_sequence = line[uncovered_region[region][0]-current_pos :]
                        next_line = True
                        line_end = False
                elif uncovered_region[region][1] in range(current_pos, current_pos+len(line)):
                    uncovered_sequence += line[:uncovered_region[region][1]-current_pos+1]
                    next_line = False
                    region += 1
                    line_end = True
                    region += [uncovered_sequence]
                elif next_line:
                    uncovered_sequence += line[:]
                    line_end = False
                else:
                    line_end = False
    return tuple(zip(uncovered_region, region))


def uncovered_sequences_mk2(filename, uncovered_region):
    """Returns tuple of (start, end, sequence) of uncovered regions

    Key inputs:
    filename --- name of target fasta file
    uncovered_region --- list of tuples of uncovered regions in sequence
    """
    target_seq = fasta_reader(filename, list)
    result = []
    for start, end in uncovered_region:
        result += [target_seq[0]['SEQUENCE'][start:end+1]]
    return tuple(zip(uncovered_region, result))

def main():
    """Main program to carry out all desired steps. """
    target_file, query_file, output_file, output_to_file, overwrite \
     = call_files()
    with open(query_file) as f:
        list_of_lines = list(f.read().splitlines())
    fasta = fasta_reader(list_of_lines, list)
    fasta = sequence_length(fasta)
    n50_tuple = n50(fasta)
    lastz_output = lastz(target_file, query_file, 
                         output_to_file = output_to_file,
                         output_file = output_file, 
                         overwrite = overwrite
                        )
    if output_to_file:
        with open(output_file) as out_file:
            covered = covered_region(out_file)
    else:
        covered = covered_region(lastz_output)
    target_length = length_sequences_fasta_file(target_file)
    uncovered_region = (target_length, covered)
    with open(target_file) as f:
        uncovered_region = uncovered_sequences_mk2(f, uncovered_region)
    print(format_final_output(uncovered_region))

def format_final_output(uncovered_region):
    """Returns a string in the format <start>:<end>    <sequence>
    """
    final_output = ''
    for start, end, sequence in uncovered_region:
        final_output += ('{}:{}\t{}\n'.format(start, end, sequence))
    return final_output

if __name__ == '__main__':
    main()