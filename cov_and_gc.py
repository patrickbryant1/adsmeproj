#! /usr/bin/env python

'''
This is a program that takes 2 inputs:
1. A .fa file of the following format:

>0 64 92
CCATAAACAAGTGCATCAACCATGCACCAACCTGTAATCAACCTGTACCATAAACAAGTGCATC

Where >0 is the ID, 64 is the sequence length and 92 is its k-mer coverage.
Columns following the coverage one may be included, but will not be considered.
2. A coverage cutoff, where sequences with coverage above that value are
selected and written to stdout.
The sum of the length of the selected sequences is also written to stdout.

The scaffolds are also selected based on <= 40 % GC-content. The GC distributions
before and after selection are plotted in a histogram.
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

def compute_gc(dna_seq,length):
    '''Function that computes GC-content for eacg DNA seq.
    Input=dna_seq (string)
    Output=gc_content (float)
    '''
    number_gc=dna_seq.count('C')+dna_seq.count('G') #counts instances of C and G
 
    gc_content=round(float(number_gc)/length, 3) #computes GC-content
 
    return gc_content
 

def select(seq_dict, cutoff):
    '''A function that takes a dict containing sequences and ids
    and a cutoff value. It then selects the sequences based on
    the cutoff value (the 3d column of the id) and <= 40% GC
    and prints them to stdout in fasta format. It also sums the
    lengths of the sequences before and after selection and returns them.
    input = dict, float
    output = fasta sequences, lists of gc before/after selection, length
    before/after selection.
    '''

    length_before = 0 #To store the sum of the lengths of the selected sequences
    length_after = 0
    gc_before = [] #Empty list to store vals for create_histo
    gc_after = [] 

    for key in seq_dict:
        split_id = key.split(' ')
        coverage = float(split_id[2])
        length = int(split_id[1])
        length_before += length
        gc_content = compute_gc(seq_dict[key], length)
        gc_before.append(gc_content)
        if coverage >= cutoff and gc_content <= 0.40:
            length_after += length
            gc_after.append(gc_content)
            print key, '\n', seq_dict[key]
        else:
            continue
    

    return (gc_before, gc_after, length_before, length_after)


def create_histo(gc_before, gc_after, length_before, length_after):
    '''A function that creates a histogram of a given list of values.
    Input = lists
    Output = histogram
    '''

    plt.hist(gc_before, histtype='stepfilled', color='b', alpha=0.5, label='Before') #plots the coverage list as a transparent (alpha=0.5) histogram
    plt.hist(gc_after, histtype='stepfilled', color='r', label='After') #plots the coverage list as a histogram 
    plt.title('GC histogram')
    plt.xlabel('%GC')
    plt.ylabel('Frequency')
    plt.text(0.45,700, length_before)
    plt.text(0.35,50, length_after)
    
    pp = PdfPages('gc_hist')
    plt.savefig(pp, format = 'pdf')
    pp.close()
    
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
    (gc_before, gc_after, length_before, length_after) = select(seq_dict, cutoff)
except IOError:
    print 'Could not select sequences based on cutoff:', cutoff
try:
    create_histo(gc_before, gc_after, length_before, length_after)
except IOError:
    print 'Could not create histogram'

