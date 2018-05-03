#!/usr/bin/env python3

import re
import sys
import subprocess
import os.path

translation_table = {'UUU': 'F', 'CUU': 'L', 'AUU': 'I', 'GUU': 'V', 
					 'UUC': 'F', 'CUC': 'L', 'AUC': 'I', 'GUC': 'V',  
					 'UUA': 'L', 'CUA': 'L', 'AUA': 'I', 'GUA': 'V', 
					 'UUG': 'L', 'CUG': 'L', 'AUG': 'M', 'GUG': 'V', 
					 'UCU': 'S', 'CCU': 'P', 'ACU': 'T', 'GCU': 'A', 
					 'UCC': 'S', 'CCC': 'P', 'ACC': 'T', 'GCC': 'A', 
					 'UCA': 'S', 'CCA': 'P', 'ACA': 'T', 'GCA': 'A', 
					 'UCG': 'S', 'CCG': 'P', 'ACG': 'T', 'GCG': 'A', 
   					 'UAU': 'Y', 'CAU': 'H', 'AAU': 'N', 'GAU': 'D', 
   					 'UAC': 'Y', 'CAC': 'H', 'AAC': 'N', 'GAC': 'D', 
   					 'UAA': 'Stop', 'CAA': 'Q', 'AAA': 'K', 'GAA': 'E', 
   					 'UAG': 'Stop', 'CAG': 'Q', 'AAG': 'K', 'GAG': 'E', 
   					 'UGU': 'C', 'CGU': 'R', 'AGU': 'S', 'GGU': 'G', 
   					 'UGC': 'C', 'CGC': 'R', 'AGC': 'S', 'GGC': 'G', 
   					 'UGA': 'Stop', 'CGA': 'R', 'AGA': 'R', 'GGA': 'G', 
   					 'UGG': 'W', 'CGG': 'R', 'AGG': 'R', 'GGG': 'G'}


def dna_to_rna(dna_sequence):
	"""Convert DNA sequence to RNA.

	Key inputs:
	sequence --- DNA sequence as a string."""
	return dna_sequence.replace('T', 'U')

def find_open_reading_frames(AA_sequence, STOP_codon_representative = '*'):
	"""Finds the open reading frames in an AA sequence
	splitting the possible ORFs at Stop codons. Stop 
	codons are recognised by a user defined input 
	(default is *)

	Key arguments:
	AA_sequence --- a string of Amino Acids represented
	by a single character or 'Stop'.

	Optional arguments:
	STOP_codon_representative --- the string which STOP codons are represented
	by. This is the symbol the program recognises to split ORFs at."""
	assert isinstance(AA_sequence, str), 'AA_sequence must be a string.'
	#print('AA seq len: {}'.format(len(AA_sequence)))
	#print('AA sequence: {}'.format(AA_sequence))
	open_reading_frames = {}
	for position in range(len(AA_sequence)):
		if AA_sequence[position] == 'M':
			key = 'Frame_{}'.format(len(open_reading_frames) + 1)
			current_frame = AA_sequence[position]
			#print(AA_sequence[position:])
			for continued_pos in range(position + 1, len(AA_sequence)):
				if AA_sequence[continued_pos] not in STOP_codon_representative:
					current_frame += AA_sequence[continued_pos]
				else:
					key = 'Frame_{}'.format(len(open_reading_frames) + 1)
					open_reading_frames[key] = current_frame
					break
	return open_reading_frames

def reverse_complement(sequence, DNA = True):
	"""Finds the reveerse complement of a DNA or RNA sequence.

	Key inputs:
	sequence --- a nucleic acid sequence with each base represented
	by single characters.
	DNA --- a boolean denoting sequence being DNA or RNA (false is RNA).
	"""
	if DNA:
		base_pairs = {'A':'T', 'G':'C', 'T':'A', 'C':'G'}
	else:
		base_pairs = {'A':'U', 'G':'C', 'U':'A', 'C':'G'}
	seq_loq = sequence
	for base in base_pairs.keys():
		seq_loq = seq_loq.replace(base, base_pairs[base].lower())
	complement = seq_loq.upper()
	return complement[::-1]

def find_reading_frames(sequence, DNA = True):
	"""From a sequence of nucleotides represented by a single
	character, returns the 6 possible reading frames.

	Key inputs:
	dna_sequence --- a string of bases represented by single 
	characters. 
	DNA --- a boolean to denote DNA or RNA sequence."""
	assert isinstance(sequence, str), 'Input must be a string.'
	reading_frames = {}
	for i in range(1,4):
		key_positive = i
		key_reverse = -i
		reading_frames[key_positive] = sequence[i-1:]
		reading_frames[key_reverse] = reverse_complement(sequence, DNA)[i-1:]
	return reading_frames

def find_potential_peptides(dna_sequence, translation_table = translation_table,
	represent_STOP_codon_by = '*'):
	"""Finds the open reading frames of a DNA sequence.

	Key inputs:
	dna_sequence --- a string representing a dna sequence
	represent_STOP_codon_by --- the string used to reprsent STOP
	codons in the translated sequence. Default is '*'."""
	reading_frames = find_reading_frames(dna_sequence)
	rna_sequences = {}
	open_reading_frames = {}
	num_frames_on_strand = int(len(reading_frames)/2)
	for frame_index in range(-num_frames_on_strand, num_frames_on_strand + 1):
		if frame_index == 0:
			continue
		rna_sequences[frame_index] = dna_to_rna(reading_frames[frame_index])
		aa_sequence = translate_codon(rna_sequences[frame_index], 
			translation_table, represent_STOP_codon_by)
		open_reading_frames[frame_index] = find_open_reading_frames(aa_sequence, 
			represent_STOP_codon_by)
	return open_reading_frames

def translate_codon(RNA_sequence, translation_table, include_STOP_as = '*'):
	"""Translates an RNA sequence into amino acids, with an option to include 
	STOP codons as a specific single symbol (default is '').

	Key inputs:
	RNA_sequence --- a string of RNA bases (ACGU)
	translation_table --- a dictionary with keys of codons and values the
	associated amino acid.
	include_STOP_as --- the string to represent Stop codons (default is '*').
	"""
	protein = ''
	for i in range(0, len(RNA_sequence) - 2, 3):
		codon = RNA_sequence[i:i+3]
		amino_acid = translation_table[codon]
		if amino_acid == 'Stop':
			amino_acid = include_STOP_as
		protein += amino_acid
	return protein

def fasta_to_dict(fasta_file_name):
	"""Writes FASTA file to a dictionary with each key the name
	of the sequence in the FASTA file and the corresponding value
	the associated sequence converted to a single string.

	Key arguments:
	fasta_file_name --- filename of FASTA file that contains the 
	sequences for analysis.
	"""
	if type(fasta_file_name) != str:
		raise TypeError('Inputs must be strings.')

	try:
		fasta_dict = {}
		with open(fasta_file_name) as fasta:
				for line in fasta:
					if re.match(r'^>', line):
						key = line.strip()[1:]
						fasta_dict[key] = ''
						loc_line = ''
					else:
						loc_line += line.strip()
						fasta_dict[key] = loc_line
	except IOError:
		raise IOError('File inaccesible')
	except FileNotFoundError:
		raise FileNotFoundError('File does not exist.')
	except Exception as err:
		raise Exception("{}. \nExiting function.".format(err))
	finally:	
		return fasta_dict

def create_translation_table(filename):
	"""Given a txt file of translations, returns a dictionary with 
	keys corresponding to the codons that translate and values of the
	associated amino acids.

	Key arguments:
	filename --- a text file formatted 'codon AA codon AA' allowing
	for whitespace between entries."""
	translation_table = {}
	with open(filename) as source:
		for line in source:
			list_of_translations_local = line.split()
			for i in range(len(list_of_translations_local) - 1):
				elem = list_of_translations_local[i]
				next_elem = list_of_translations_local[i + 1]
				if len(elem) == 3:
					translation_table[elem] = next_elem
	return translation_table

def fasta_line(fasta_file_name, output_file):
	"""Writes FASTA file to text file with each sequence on one line (i.e. as 
	a single continuous string).

	Key arguments:
	fasta_file_name --- filename of FASTA file that contains the sequences 
	for analysis
	output_file --- the name of the file to write the manipulated strings to
	"""
	if type(fasta_file_name) != str or type(output_file) != str:
		raise TypeError('Inputs must be strings.')

	try:
		with open(fasta_file_name) as fasta:
			with open(output_file, 'w+') as output:
				for line in fasta:
					if re.match(r'^>', line):
						try:
							loc_line += '\n'
							output.write(loc_line)
						except:
							pass
						loc_line = ''
						output.write(line)
					else:
						loc_line += line.strip()
				output.write(loc_line)
	except IOError:
		raise IOError('File inaccesible')
	except FileNotFoundError:
		raise FileNotFoundError('File does not exist.')
	except Exception as err:
		print("{}. \nExiting function.".format(err))
		sys.exit()
	finally:	
		return

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
    except KeyError:
        raise KeyError('Key not in elements in list. ')
    except Exception as err:
        raise Exception(err)
    
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
        cmd = 'lastz {} {} --output={}'.format(target, query, output_file)
    else:
        cmd = 'lastz {} {}'.format(target, query)
    
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
