#!/usr/bin/env python3

import sys

def create_translation_table(filename):
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


def translate_codon(RNA_sequence, translation_table):
	protein = ''
	for i in range(0, len(RNA_sequence) - 2, 3):
		codon = RNA_sequence[i:i+3]
		amino_acid = translation_table[codon]
		if amino_acid == 'Stop':
			amino_acid = ''
		protein += amino_acid
	return protein

if __name__ == '__main__':
	translation_table = create_translation_table(sys.argv[1])
	with open(sys.argv[2]) as f:
		for line in f:
			print(translate_codon(line, translation_table))