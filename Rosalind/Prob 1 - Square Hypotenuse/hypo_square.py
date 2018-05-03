#!/usr/bin/python3
# From the included text file, calculates the square of the required 
# hypothesis of the two integers contained therein to form a right 
# angled triangle.

if __name__ == '__main__':
    with open('rosalind_ini2.txt') as f:
        for line in f:
             arguments = line.split()
             arguments = list(map(int, arguments))
             assert len(arguments) == 2, 'Requires 2 arguments. '

             # The lines 11:15 are equivalent
             squares = list(map(lambda x: x ** 2, arguments))
             squares = [i ** 2 for i in arguments]
             squares = []
             for arg in arguments:
                 squares += [arg * arg]
                 
             hypo = sum(squares)
    print(hypo)
