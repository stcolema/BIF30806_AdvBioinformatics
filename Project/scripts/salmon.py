#!/usr/bin/env python3
"""
Salmon indexing and quantifying of RNA-seq reads

Author: Thijs Maas
Student nr: 921116537120
"""
from subprocess import check_call
from os.path import exists
from argparse import ArgumentParser
import tempfile

def sindex(transcriptome, outdir):
    """Runs salmon index and write to outdir.
    
    Keyword arguments:
    transcriptome -- string, fasta file of the transcriptome
    outdir -- string, directory name that will be made.
    Returns:
    None
    """
    if exists(outdir+'/sa.bin'):
        print('Index was already present, now continuing with quantification!')
        return None
    path = '/local/prog/Salmon-latest_linux_x86_64/bin/'
    cmd = "salmon index -t {} -i {}".format(transcriptome, outdir)
    check_call(path+cmd, shell=True)
    print("Index succesfully made, now continuing with quantification!")
    return None
    
def squant(index, outdir, bootstrap, unmated=False, pairs=False):
    """Runs salmon quant on unmated to paired end files.
    
    Keyword arguments:
    index -- string, folder name of index location of this RNA-seq library.
    outdir -- string, directory name that will be made.
    bootstrap --int, number of bootstraps done.
    unmated -- string, a path to a single-end RNA-seq fastq file
    pais -- string, two paths to paired-end RNA-seq fastq files
    Returns:
    None
    """
    #if exists(outfile):
    #    return None
    path = '/local/prog/Salmon-latest_linux_x86_64/bin/'
    if pairs == False:
        cmd = "salmon quant -i {} -r {} -o {}".format(index, unmated, outdir)
        
    elif pairs != False:
        cmd = "salmon quant -i {} -1 {} -2 {} -o {} ".format(index, pairs[0], pairs[1], outdir)

    cmd2 = " -p 12 --libType A --numBootstraps {}".format(bootstrap)
    check_call(path+cmd+cmd2, shell=True)
    return None
    
    
if __name__ == "__main__":
    
    # argument parsing
    AP = ArgumentParser(description='This script is used to first index, and \
    then quantify a RNA-seq reads library by pseudomapping using Salmon. The \
    user can use either single-end or paired-end reads, not both.')
    AP.add_argument('transcriptome', help='used for indexing the reads')
    AP.add_argument('-r', '--unmatedReads', help='file with unpaired reads', \
                    default=False) #optional
    AP.add_argument('-m', '--mates', help='space seperated list of paired \
                    reads', nargs=2, default=False) #optional
    AP.add_argument('-b', '--bootstrap', help='number of bootstraps, def=0', \
                    default=0) #optional
    args = AP.parse_args()
    
    # Check if files are correctly given
    if (args.unmatedReads != False and args.mates != False) or \
        (args.unmatedReads == False and args.mates == False):
        print('Give either 1 single end file, or 2 paired end files!')
        exit()
              
    # make a file name from the given path
    if args.unmatedReads != False:
        path = args.unmatedReads
    else:
        path = args.mates[0]
    sample = path.split('/')[-1].split('.')[0].split('_')[0]
    if 'trim' in path:
        sample += 'trim'
    
    
    # check if bootstrapping is done, change the file name
    if args.bootstrap != 0:
        sample += '_bs'+args.bootstrap
    
    # check which of 2 options are given
    if args.unmatedReads != False:
        print("Salmon for single end reads")
        sindex(args.transcriptome, 'trans_index')
        squant('trans_index', sample+'_quant', args.bootstrap, \
                   unmated=args.unmatedReads)
        
    elif args.mates != False:
        print("Salmon for paired end reads")
        sindex(args.transcriptome, 'trans_index')      
        squant('trans_index', sample+'_quant',  args.bootstrap, \
                   pairs=args.mates)
    
    
    
