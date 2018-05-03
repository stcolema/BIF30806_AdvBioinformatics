#!/usr/bin/python3
# Trnascribes DNA string into RNA string and prints each sequence
# Does this by converting all occurences of T to U

import sys

if __name__ == '__main__':
    with open(sys.argv[1]) as f:
        for line in f:
            new_line = line.replace('T', 'U')
            print(new_line)
