#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 28 11:41:14 2018

@author: shanedenecke
"""

from Bio import SeqIO
import sys
from Bio.SeqIO.FastaIO import FastaWriter
import os 

if len(sys.argv)>3:
    seqs=list(SeqIO.parse(sys.argv[1], 'fasta'))
    filt_low=int(sys.argv[2])
    filt_high=int(sys.argv[3])
elif len(sys.argv)==3:
     seqs=list(SeqIO.parse(sys.stdin, 'fasta'))
     filt_low=int(sys.argv[1])
     filt_high=int(sys.argv[2])
    
## filter fasta file
new=[]
for i in seqs:
    if (len(i.seq)>filt_low) and (len(i.seq)<filt_high):
        new.append(i)


## Write temporary file
handle = open('temp.fa', "w")
writer = FastaWriter(handle, wrap=0)
writer.write_file(new)
handle.close()


## Read in temporary file and print properly formatted fasta
x = open("temp.fa", "r")
y=x.readlines()
z=''.join(y)
if z[-1]=='\n': 
    z=z[:-1]
print (z)
os.remove("temp.fa") 