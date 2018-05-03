#!/usr/bin/env python3

"""
Author: Stephen Coleman
Student number: 940309-160-050
Contains functions to parse and provide statistics for a file formatted in the
style used by GenBank.
"""

import sys, re

def find_key_orig(record, key = 'ACCESSION'):
	"""Finds key and associated value from record.

	Key inputs:
	record --- a list of lines to search through
	key --- the section of information of interest. Only works
	for values of ACCESSION, LOCUS and OGANISM."""
	if key == 'ACCESSION' or key == 'LOCUS':
		strip_str = ' '
		num_strips = 1
	elif key == 'ORGANISM':
		strip_str = '  '
		num_strips = 2
	else:
		raise ValueError('key should be one of ACCESSION, LOCUS or ORGANISM')
	for line in record:
		if key in line:
			entry_lst = line.split(strip_str, num_strips)
			entry_lst = entry_lst[-2:]
			entry_lst[1] = entry_lst[1].strip().split(' ')
			entry_lst[1] =[item for item in entry_lst[1] if item != '']
			entry_lst[1] = ' '.join(entry_lst[1])
	entry_dict = {entry_lst[0]: entry_lst[1]}
	return entry_dict


def find_key(record, key, num_strips = 1, strip_str = ' '):
	"""Finds key and associated value from record.

	Key inputs:
	record --- a list of lines to search through
	key --- the section of information of interest."""
	part_of_section = False
	for line in record:
		if key in line or part_of_section:
			if re.match(r'^\s{0,2}[A-Z]+', line) and key not in line:
				part_of_section = False
				break
			if not part_of_section:
				entry_lst = re.sub('  +', strip_str, line)
				entry_lst = entry_lst.split(strip_str, num_strips)
				entry_lst = entry_lst[-2:]
				entry_lst[1] = entry_lst[1].strip().split(' ')
				entry_lst[1] = [item for item in entry_lst[1] if item != '']
				entry_lst[1] = ' '.join(entry_lst[1])
			else:
				loc_entry = re.sub('  +', strip_str, line)
				entry_lst[1] += '\n' + loc_entry.strip()
			
			part_of_section = True
	entry_dict = {entry_lst[0]: entry_lst[1]}
	return entry_dict

def find_sequence(record):
	"""Within a single entry in GenBank format, finds the DNA sequence

	Key inputs:
	record --- a single entry in GenBank format containing a DNA sequence
	preceded by ORIGIN and ending with // on a new line."""
	seq = ''
	in_seq = False
	for line in record:
		if 'ORIGIN' in line:
			in_seq = True
		elif line.startswith('//'):
			in_seq = False
		elif in_seq:
			current = [base for base in line.strip() if base.isalpha()]
			seq += ''.join(current)
	return seq

def split_records_gen(filename):
	"""Returns a generator of the different GenBank entries
	contained within filename.

	Key inputs:
	filename --- the file containing GenBank entries.
	"""
	current = []
	with open(filename) as f:
		for line in f:
			current += [line]
			if line.startswith('//'):
				yield current
				current = [] # Ask Sandra faoi
	if current:
		yield current

def split_records(filename):
	"""Returns a list of lists of the entries within filename.

	Key inputs:
	filename --- the name of the file containing the GenBank entries.
	"""
	record = []
	current = []
	with open(filename) as f:
		for line in f:
			current += [line]
			if line.startswith('//'):
				record += [current]
				current = []
	return record

def find_GC_content(sequence):
	"""Returns a dict with keys corresponding with the sequence titles from the
	 FASTA file	and values of the associated sequence's GC content. 

	Keyword arguments:
	fasta_file_name --- the name of the FASTA file to open for searching 
	for sequences."""
	G_count = 0
	C_count = 0
	count = 0
	loc_line = sequence.strip().upper()
	G_count += loc_line.count('G')
	C_count += loc_line.count('C')
	count += len(loc_line)
	GC_content =  float((G_count + C_count)) / count
	return GC_content

def write_fasta(record, output_name, mode = 'w', line_width = 80):
	"""Write a FASTA file containing the DNA sequence contained 
	within record.

	Key inputs:
	record --- the dictionary containing the VERSION and SOURCE (as 
	used in GenBank format) and the DNA sequence (under key SEQUENCE).
	output_name --- the name of the finished FASTA file
	mode --- the mode for writing to the FASTA file; w writes a new file
	(verwriting if necessary) and a allows appending to a pre-existing file.
	line_width --- the number of characters allowed on each line of the FASTA
	file.
	"""
	assert mode == 'w' or mode == 'a', 'Mode must be \'w\' or \'a\'. '
	with open(output_name, mode) as f:
		title = '>{} [{}]\n'.format(record['VERSION'], record['SOURCE'])
		f.write(title)
		range_sup = ((len(record['SEQUENCE']) // line_width) + 1) * line_width
		for i in range(0, range_sup + 1, line_width):
			f.write(record['SEQUENCE'].upper()[i:min(i + line_width, len(record['SEQUENCE']))])
			f.write('\n')
		f.write('\n')
	return

def write_fasta_single(record, output_name, mode = 'w'):
		"""Write a FASTA file containing the DNA sequence contained 
	within record.

	Key inputs:
	record --- the dictionary containing the VERSION and SOURCE (as 
	used in GenBank format) and the DNA sequence (under key SEQUENCE).
	output_name --- the name of the finished FASTA file
	mode --- the mode for writing to the FASTA file; w writes a new file
	(verwriting if necessary) and a allows appending to a pre-existing file.
	"""
	assert mode == 'w' or mode == 'a', 'Mode must be \'w\' or \'a\'. '
	with open(output_name, mode) as f:
		title = '>{} [{}]\n'.format(record['VERSION'], record['SOURCE'])
		f.write(title)
		f.write(record['SEQUENCE'].upper())
		f.write('\n\n')
	return

def genbank_stats(record, cols):
	"""Creates a rather ugly table of the entries in cols contained in record.

	Key inputs:
	record --- a dictionary with a (sub-)set of keys corresponding to the values
	contained in cols.
	cols --- a list of the desired statistics outputted as the column headers.
	"""
	line = '\t'.join(cols)
	print(line)
	line = []
	for col in cols:
		line += [record[col]]
	line = map(str, line)
	line = '\t'.join(line)
	print(line)

if __name__ == '__main__':
	record = split_records(sys.argv[1])

	# data = [None for i in range(len(record))] # Does not work if record is a generator.
	data = [] # Needs to be a list to order by GC content
	# data = {}
	count = 0
	for entry in record:
		data += [{}]
		# data[count] = {}
		data[count].update(find_key(entry, 'ACCESSION'))
		data[count].update(find_key(entry, 'ORGANISM', 2))
		data[count].update(find_key(entry, 'VERSION'))
		data[count].update(find_key(entry, 'SOURCE'))
		data[count]['SEQUENCE'] = find_sequence(entry)
		data[count]['GC Content'] = find_GC_content(data[count]['SEQUENCE'])
		data[count]['GC Content'] = '{:.2%}'.format(data[count]['GC Content'])
		data[count]['Seq Length'] = '{}bp'.format(len(data[count]['SEQUENCE']))
		count += 1

	# Sort the list in order of descending GC content
	new_list = sorted(data, key=lambda k: k['GC Content'], reverse=True)

	#Define boolean which declares if this entry is the first being written
	#to the FASTA file
	first = True
	for entry in new_list:
		#If not the first entry append entry to pre-existing file
		if first:
			mode = 'w'
		else:
			mode = 'a'
		write_fasta_single(entry, sys.argv[2], mode)
		first = False

	#Want to print table with headers as declared in cols_of_interest
	cols_of_interest = ['ACCESSION', 'SOURCE', 'GC Content', 'Seq Length']
	for entry in new_list:
		genbank_stats(entry, cols_of_interest)
		print('\n')