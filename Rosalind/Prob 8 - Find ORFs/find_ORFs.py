#!/usr/bin/env python3

import re, sys


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

if __name__ == '__main__':
	dna_seq = fasta_to_dict(sys.argv[1])
	represent_STOP_codon_by = '*'

	reading_frames = {}
	rna_sequences = {}
	ORFs = {}

	for key in dna_seq.keys():
		ORFs[key] = find_potential_peptides(dna_seq[key])


	unique = set([])

	for key in ORFs.keys():
		# print('Key is: {}'.format(key))
		for reading_frame in ORFs[key].keys():
			# print('Reading frame {}'.format(reading_frame))
			for open_reading_frame in ORFs[key][reading_frame]:
				unique.add(ORFs[key][reading_frame][open_reading_frame])
				# if ORFs[key][reading_frame][open_reading_frame] not in unique:
				# 	unique += [ORFs[key][reading_frame][open_reading_frame]]
			# 	print(ORFs[key][reading_frame][open_reading_frame])
			# print('\n')

	for peptide in unique:
		print(peptide)


	