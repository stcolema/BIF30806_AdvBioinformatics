#!/usr/bin/python3
# Call from command line with argument of file name
# Outputs the number of A C G T bases in the file (in said order)

import sys

if __name__ == '__main__':
    nucleotides = 'AGCT'
    dna_count = dict.fromkeys(nucleotides, 0)
    with open(sys.argv[1]) as f:
         for line in f:
             for ch in line:
                 if ch in nucleotides:
                     dna_count[ch] += 1
    print('{A} {C} {G} {T}'.format(**dna_count))
