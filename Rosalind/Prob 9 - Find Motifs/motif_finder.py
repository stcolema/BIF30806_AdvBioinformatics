#!/usr/bin/env python3

import sys

def motif_position(sequence, motif, start = 1):
	"""Find position of motif in sequence. Returned position begins counting 
	from start.

	Key inputs:
	sequence --- a string 
	motif --- substring to search for
	start --- integer limited to 0 or 1; defines result (i.e. is the first
	position 0 or 1)."""

	assert isinstance(sequence, str), 'Sequence input must be a string. ' \
									  'Is type {}'.format(type(sequence))
	assert isinstance(motif, str), 'Motif input must be a string. ' \
								   'Is tpye {}'.format(type(motif))
	assert len(sequence) >= len(motif)
	if start not in [0,1]:
		raise ValueError('start should be 0 or 1.')

	all_found = False
	curr_pos = -1
	positions = []
	while not all_found:
		positions += [sequence.find(motif, curr_pos + 1) + start]
		curr_pos = positions[-1]
		if positions[-1] == -1 + start:
			positions = positions[:-1]
			all_found = True
	return positions

if __name__ == '__main__':
	lines =[]
	with open(sys.argv[1]) as f:
		for line in f:
			lines += [line.strip()]
	lines.sort(key = len)
	motif = lines[0]
	sequence = lines[1]
	positions = map(str, motif_position(sequence, motif))
	print(' '.join(positions))