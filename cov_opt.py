#! /usr/bin/env python

'''
This is a program that takes 2 inputs:
1. A .fa file of the following format:

>0 64 92
CCATAAACAAGTGCATCAACCATGCACCAACCTGTAATCAACCTGTACCATAAACAAGTGCATC

Where >0 is the ID, 64 is the sequence length and 92 is its k-mer coverage.
Columns following the coverage one may be included, but will not be considered.
2. A coverage cutoff, where sequences with coverage above that value are
counted and written to stdout. 
The coverage cutoff and the sum of the length of the selected sequences is also written to stdout.
'''

#Modules
import sys
import os
import pdb
import matplotlib.pyplot as plt
import numpy
from matplotlib.backends.backend_pdf import PdfPages

#Functions
def read_fasta(scaffolds):
    '''A function that opens a file and creates a dict of it
    input = file
    output = seq_dict (dict)
    '''
    seq_dict = {}        
    for line in scaffolds:
        line = line.rstrip()

        if len(line) == 0: #Makes sure list index is not out of range
            continue
        elif line[0] == '>':
            val='' #Resets val to an empty string
            key = line
            seq_dict[key]=val
        else:
            val = line
            seq_dict[key] +=val
                            
    return seq_dict


def above_cutoff(seq_dict, cutoff):
    '''A function that takes a dict containing sequences and ids
    and a cutoff value. It then selects the sequences based on
    the cutoff value and the 3d column of the id and prints them
    to stdout in fasta format. It also sums the lengths of the
    selected sequences and prints it to stdout.
    input = dict, float
    output = fasta sequences and lists of coverages before/after selection
    '''

    #To store the sum of the lengths of the selected sequences
    length_calculation = 0 
    #To count number of scaffolds after selection
    n_after = 0
    
    for key in seq_dict:
        split_id = key.split(' ')
        coverage = float(split_id[2])
        if coverage >= cutoff:
            n_after+=1
            length_calculation += float(split_id[1])
        else:
                continue
    
    print cutoff, n_after,length_calculation




#Main program

try:
    scaffolds = open(sys.argv[1], 'r') #opens file given in argument
    cutoff = float(sys.argv[2]) #Assigns second argument to cutoff
    #Checks if file is empty
    if os.stat(sys.argv[1]).st_size == 0:
       print "Empty file."
except IOError:
    print 'Cannot open', sys.argv[1],
    print 'or cutoff is not float'
else:
    try:
        seq_dict=read_fasta(scaffolds) #creates dict of file
    except IOError:
        print 'Could not read sequences in', sys.argv[1]
try:
    above_cutoff(seq_dict, cutoff)
except IOError:
    print 'Could not select sequences based on cutoff:', cutoff
    

