#!/usr/bin/env python3
"""
Author: Lisa van der Goes
Student number: 961023268010
Script for P5
"""

from sys import argv

def get_infile():
    """ Returns input files ref and related in fasta format
    """
    if len(argv) < 3:
        InputFile = input('What are the input files: ref.fasta, \
        related.fasta?')
    else:
        InputFile_ref = argv[1]
        InputFile_related = argv[2]
    return InputFile_ref, InputFile_related
    
def parse_fasta(InputFile):
    """ Returns dictionary with name as key and sequence as label
    
    InputFileName -- file in fasta format
    """
    seqs = {}
    with open(InputFile) as InFile:
        for line in InFile:
            if not line.strip():
                continue
            if line.startswith('>'):
                label = line.strip()[1:].split()[0]
                seqs[label] = ""
            else:
                seqs[label] += line.strip()
        return seqs
    

def calc_length_seq(seqs):
    """ Returns a list of names and sequences
    
    seqs -- a dictionary with names ad keys and sequences as labels
    """
    tuple_sequences = seqs.items()
    names_lengths = []
    for i in tuple_sequences:
        length_seq = len(i[1])
        names_lengths.append([i[0], length_seq])
    return names_lengths
    
def calc_hamm_dist (sequence1, sequence2):
    """ Return a integer of the Hamming length
    
    sequence1 -- string of ref sequence
    sequence2 -- string of related sequence
    """
    if len(sequence1) == len(sequence2):
        for i in range(len(sequence1)):
            if sequence1[i] != sequence2[i]:
                count_hamm += 1
    return count_hamm
    
def percent_of_identity(sequence, hamming_distance):
    """ Return a float of the identity persentage
    
    sequence -- string of one of the input sequences
    hamming_distance -- integer of the calculated hamming distance
    """
    alignment_length = len(sequence)
    number_of_identical_positions = alignment_length - hamming_distance
    identity = number_of_identical_positions / alignment_length * 100
    return identity

def parse_needle(InputNeedle, ref_name, related_name):
    """ Returns the ref and related sequences as  a list
    
    InputNeedle -- input needle file
    ref_name -- string of ref name
    related_name -- string of related name
    """
    #not ready yet
    with open(InputNeedle) as InNeedle:
        seqs = ''
        for line in InNeedle:
            if not line.strip():
                continue
            if line.startswith(ref_name):
                seq_ref = line
                seq_ref = line.strip().split()[2:-1]
                line = Inneedle.readline().readline().strip.split()[2:-1]
                tot_line += line
                ref_seq = parsedfile.append(tot_line)
            if line.startswith(related_name):
        
        #~ return ref_seq, related_seq
                
                
def give_output(related_lengths, ref_length):
    """Prints the output table
    
    related_lengths -- lists of names and sequence lengths
    ref_length -- list of ref name and the sequence length
    """
    print ('Sequence1\tLength\tSequence2\tLength\tHamm\tIdent\n')
    for i in related_lengths:
        #~ ref_seq, related_seq = parse_needle('out.needle', \
        #~ ref_length[0], i[0])
        #~ hamm_dist = calc_hamm_dist(ref_seq, related_seq)
        #~ identity = percentent_of_identity(related_seq, hamm_dist)
        #~ print (ref_length[0],'\t',ref_length[1],'\t', i[0],'\t', i[1],\
        #~ '\t', hamm_dist, '\t', identity, '\n')
        
        #test print
        print (ref_length[0],'\t',ref_length[1],'\t',i[0],'\t',i[1],'\n')        

if __name__ == '__main__':
    
    
    # get the input files
    InputFile_ref, InputFile_related = get_infile()
    
    # parse the two files
    related_seqs = parse_fasta(InputFile_related)
    ref_seq = parse_fasta(InputFile_ref)
    
    # calculate the lengths of the files
    related_lengths = calc_length_seq(related_seqs)
    ref_length = calc_length_seq(ref_seq)
    for i in ref_length:
        reflength = i # otherwise list in list
    
    give_output(related_lengths, ref_length)
