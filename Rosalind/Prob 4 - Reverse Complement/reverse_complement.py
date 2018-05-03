#!/usr/bin/env python3

import sys


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


if __name__ == '__main__':
	with open(sys.argv[1]) as f:
		for line in f:
			print(reverse_complement(line))