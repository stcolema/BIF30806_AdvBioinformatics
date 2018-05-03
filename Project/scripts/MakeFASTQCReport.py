#!/usr/bin/env python3

"""
Author: Luuk Peters (& Thijs Maas)
Student number: 950513-649-020

Description: 
"""

from sys import argv
from subprocess import check_call

def RunFASTQ_Trimmer(file_name, _quality_threshold, _min_length):
    """
    The first Line of a docstring should ttopell what the function returns    
    
    """    
    # check file type
    if pathToFile.endswith('.fastq.gz'):
        filetype = 'gzipped'
    if pathToFile.endswith('.fastq'):
        filetype = 'fastq'
    
    # Run command
    if filetype is 'fastq':
        print("Quality timming...")
        cmd = 'fastq_quality_trimmer -i {} -o {} -t {} -l {}'.\
        format(file_name, (file_name.split("/")[-1].strip(".fastq") + \
        "_Qtrimmed.fastq"), _quality_threshold, _min_length)
        check_call(cmd, shell=True) 
    
    if filetype is 'gzipped':
        tempfile = 'filefortrimming.temp'
        cmd = 'gunzip {} -c > {}'.format(file_name, tempfile)
        print("First gunzip this file...")
        check_call(cmd, shell=True)
        print("Quality timming...")
        cmd = 'fastq_quality_trimmer -i {0} -o {1} -t {2} -l {3} | rm {0}'.format(\
        tempfile, (file_name.split("/")[-1].strip(".gz.fastq") + \
        "_Qtrimmed.fastq"), _quality_threshold, _min_length)
        check_call(cmd, shell=True)  
    return None
        

if __name__ == "__main__":
    
    pathToFile = argv[1]
    print(pathToFile)
    quality_threshold = argv[2]
    min_length = argv[3]    
    
    RunFASTQ_Trimmer(pathToFile, quality_threshold, min_length)    

    
    
    
