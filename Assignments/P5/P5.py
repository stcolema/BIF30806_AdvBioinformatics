#!/usr/bin/env python3

"""
Author: Stephen Coleman
Student number: 940309-160-050
Python Mock Exam BIF-30806
"""

import sys
import subprocess
import os.path
import re

def fasta_reader(filename):
    """Converts a FASTA file into a dict with keys 
    corresponding to the sequence names and entries
    a single string of nucleotides or amino acids.
    
    Key inputs:
    filename --- the name of the FASTA file to be 
    converted into a dict.
    """
    result = {}
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                current = line[1:]
                result[current] = {}
                result[current]['SEQUENCE'] = ''
            else:
                result[current]['SEQUENCE'] += line
    return result
        
def sequence_length(fasta_dict):
    """Adds a key of 'SEQ LENGTH' to the dictionary encoding the FASTA 
    file and gives a value of the length of the sequence string to this 
    key.
    
    Key inputs:
    fasta_dict --- a dictionary output of the function fasta_reader as 
    performed on a FASTA file.
    """
    for seq in fasta_dict.keys():
        fasta_dict[seq]['SEQ LENGTH'] = len(fasta_dict[seq]['SEQUENCE'])
    return fasta_dict

def hamming_dist(sequence1, sequence2):
    """Calculates the hamming distance between two sequences of equal 
    length.
    
    Key inputs:
    sequence1 --- a string of DNA or Amino Acidss represented by single 
    characters.
    sequence2 --- a string similar to sequence1 and of the same length.
    """
    assert len(sequence1) == len(sequence2), 'Unequal sequence length. ' \
        + '{} compared to {}. '.format(len(sequence1), len(sequence2))
           
    dist = 0
    for sym1, sym2 in zip(sequence1, sequence2):
        if sym1 != sym2:
            dist += 1

    # for pos in range(len(sequence1)):
    #     if sequence1[pos] != sequence2[pos]:
    #         dist += 1
            
    return dist

def needle(sequence1, sequence2, gap_open_pen, gap_extend_pen, output_name):
    """Run needle program on two fasta files.

    Key inputs:
    sequence1 --- a FASTA file.
    sequence2 --- another FASTA file.
    gap_open_pen --- the penalty to be used by NEEDLE for penalising alignment
    score due to the presence of gaps.
    gap_extend_pen --- the penalty for each additional contiguous gap.
    output_name --- the filename for the NEEDLE output (must be .needle).
    """
    cmd = 'needle {} {} -gapopen {} -gapexten {} -outfile {}'\
          .format(sequence1, sequence2, gap_open_pen, gap_extend_pen, \
                  output_name)
    # Below line would include output of needle - instead we write it to a file
    # and use that.
    # printed_output = subprocess.check_output(cmd, shell=True) # Unnecessary

    res = subprocess.check_call(cmd, shell=True)
    return res

#def percentage_identity(sequence1, sequence2):

def needle_seperate_records(filename):
    """Seperates out each alignment in the needle output
    into seperate entries in a list.
    
    Key inputs:
    filename --- the needle output filename.
    """
    
    #Each entry in the needle output contains two of the below strings
    entry_break = '#======================================='
    
    #Document ends with two of the below strings
    end = '#---------------------------------------'
    
    start = False
    entries = []
    curr_entry = 0
    current = []
    
    with open(filename) as f:
        for line in f:
            
            #The output file begins with some information about the needle run.
            #The function only begins collecting data after this initial 
            #description.
            if entry_break in line or end in line:
                if start:
                    count += 1
                else:
                    start = True
                    count = 0
            
            if start:
                current += [line.strip()]        
            
                if count == 2:
                    entries += [current]
                    curr_entry += 1
                    count = 0
                    current = []                
    return entries

def extract_pairwise_alignments(needle_entry):
    """For each needle output entry returns a dictionary
    with keys corresponding to the sequence names and values of the 
    aligned sequences.
    
    Key inputs:
    needle_entry --- a single pariwise alignment from a needle output 
    as formatted by needle_seperate_records().
    """
    result = {}
    for line in needle_entry:
        if line.startswith('#'):
            if '1:' in line:
                sequence1 = line[4:].strip()
                result[sequence1] = ''
            elif '2:' in line:
                sequence2 = line[4:].strip()
                result[sequence2] = ''
            continue
            
        if sequence1 in line:
            entry = line.strip()
            entry = entry.replace(sequence1, '')
            symbols = re.compile(r'[^a-zA-Z\-]')
            entry = symbols.sub('', entry)
            result[sequence1] += entry
        
        if sequence2 in line:
            entry = line.strip()
            entry = entry.replace(sequence2, '')
            symbols = re.compile(r'[^a-zA-Z\-]')
            entry = symbols.sub('', entry)
            result[sequence2] += entry
        result['Sequence1'] = '{}'.format(sequence1)
        result['Sequence2'] = '{}'.format(sequence2)
    return result
    
def percent_identity(pairwise_alignment):
    """Returns the percent identity of a pairwise alignment.
    
    Key inputs:
    pairwise_alignment --- a dictionary of two aligned sequences. 
    Each entry is a sequence with gaps represented by '-' as required 
    to ensure alignment and equal length.
    """
    sequences = list(pairwise_alignment.keys())
    sequences.remove('Sequence1')
    sequences.remove('Sequence2')
    assert len(pairwise_alignment[sequences[0]]) ==  \
           len(pairwise_alignment[sequences[1]]), 'Length of sequences'\
                                                   + ' different. '

    hamm_dist = hamming_dist(pairwise_alignment[sequences[0]], \
                             pairwise_alignment[sequences[1]])
                             
    alignment_length = len(pairwise_alignment[sequences[0]])
    identity = ((alignment_length - hamm_dist) / alignment_length) * 100
    return hamm_dist, identity

def print_output(comparison):
    """Print tab delimited table from dicitonary.

    Key inputs:
    comparison --- dictionary of dictionaryies with keys of 
    each alignment comparison and sub-keys including: 
    Sequence1 : Name of sequence as recorded in FASTA file
    <sequence1_name> : sequence
    Sequence2: Name of aligned sequence
    <sequence2_name> : second sequence
    Hamm : hamming distance between aligned sequences
    Ident : percentage identity between aligned sequences
    """
    cols = ['Sequence1', 'Length', 'Sequence2', 'Length', 'Hamm',
     'Ident',]
    header = '\t'.join(cols)
    output = header + '\n'
     
    for line in comparison:
        
        sequence1_name = comparison[line]['Sequence1']
        sequence1_length = len(comparison[line][sequence1_name].replace('-', ''))
        
        sequence2_name = comparison[line]['Sequence2']
        sequence2_length = len(comparison[line][sequence2_name].replace('-', ''))
        
        hamm = comparison[line]['Hamm']
        ident = comparison[line]['Ident']
        
        new_line = [sequence1_name, sequence1_length, sequence2_name, \
         sequence2_length, hamm, '{:.2f}'.format(ident)]
        
        new_line = map(str, new_line)
        new_line = '\t'.join(new_line)
        
        output += new_line + '\n'
    
    print(output)

if __name__ == '__main__':
    
    # Try to read command line filenames into fasta_reader function
    try:
        reference_file = sys.argv[1]
        related_file = sys.argv[2]
        reference = fasta_reader(reference_file)
        related = fasta_reader(related_file)
    except IndexError:
        reference_file = input('Please input reference file for analysis: ')
        related_file = input('Please input a file for related sequences: ')
        reference = fasta_reader(reference_file)
        related = fasta_reader(related_file)
    except Exception as err:
        raise Exception('{}'.format(err))
    
    # Add the sequence length to each entry in the dictionaries
    reference = sequence_length(reference)
    related = sequence_length(related)
    
    # Define the inputs for the needle program
    try:
        gap_open_penalty = sys.argv[3]
        gap_extension_penalty = sys.argv[4]
    except:
        gap_open_penalty = 8
        gap_extension_penalty = 0.5
    
    output_file = 'out.needle' # Check if this exists?
    
    # Call needle
    needle_check = needle(reference_file, related_file, gap_open_penalty, 
           gap_extension_penalty, output_file)
    
    needle_records = needle_seperate_records(output_file)
    
    comparison = {}
    for record in needle_records:
        entry = extract_pairwise_alignments(record)
        sequence1 = entry['Sequence1']
        sequence2 = entry['Sequence2']
        key = '{} : {}'.format(sequence1, sequence2)
        comparison[key] = entry
        hamm_dist, identity = percent_identity(entry)
        comparison[key]['Hamm'] = hamm_dist
        comparison[key]['Ident'] = identity
    
    print_output(comparison)