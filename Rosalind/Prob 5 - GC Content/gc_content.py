#!/usr/bin/python3
# From a Fasta file with multiple sequences, finds that which has the the highest GC content

import re, sys

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
						loc_line += line.replace('\n', '')
				output.write(loc_line)
	except IOError:
		raise IOError('File inaccesible')
	except FileNotFoundError:
		raise FileNotFoundError('File does not exist.')
	finally:	
		return

def find_GC_content(fasta_file_name):
	"""Returns a dict with keys corresponding with the sequence titles from the
	 FASTA file	and values of the associated sequence's GC content. 

	Keyword arguments:
	fasta_file_name --- the name of the FASTA file to open for searching 
	for sequences."""
	with open(fasta_file_name) as fasta:
		GC_content = {}
		for line in fasta:

			# Each line (bar the last) ends with '\n'
			loc_line = line.replace('\n', '')

			# Finds '>' at opening of line (FASTA seq title)
			if re.match(r'^>', loc_line):
				GC_content[loc_line] = 0
				G_count = 0
				C_count = 0
				count = 0
				current = loc_line
			else:
				G_count += loc_line.count('G')
				C_count += loc_line.count('C')
				count += len(loc_line)
				GC_content[current] =  float((G_count + C_count)) / count
	return GC_content


if __name__ == '__main__':
	GC_content = find_GC_content(sys.argv[1])

	# Finds key with highest associated value
	highest_GC_content_seq = max(GC_content, key=GC_content.get)
	print('{}: {:.6%}'.format(highest_GC_content_seq,
		GC_content[highest_GC_content_seq]))
	fasta_line(sys.argv[1], sys.argv[2])