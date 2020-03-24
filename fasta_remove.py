#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  9 18:22:03 2019

@author: shanedenecke
"""

import os
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter
import sys 



os.path.dirname(os.path.abspath(__file__))
#os.chdir('/home/shanedenecke/Dropbox/wp5_midgut_uptake/ABC_phylo/ABC_id/Nv/trinity/alignments')

#sys.argv=['','./Dm_Nv_combined_ABC_NBD_filtered.fa','./Manual_domain_missing.txt']

##define argument variables
if len(sys.argv)>2:
    fa=list(SeqIO.parse(str(sys.argv[1]), 'fasta'))
    rem=[line.rstrip('\n') for line in open(sys.argv[2])]
#elif len(sys.argv)==2:
    #fa=list(SeqIO.parse(str(sys.stdin, 'fasta'))
    #rem=[line.rstrip('\n') for line in open(sys.argv[1])] 


clean=[]
for record in fa: 
   if record.id not in rem:
       clean.append(record)

## Write temporary file
handle = open('temp.fa', "w")
writer = FastaWriter(handle, wrap=0)
writer.write_file(clean)
handle.close()

## Read in temporary file and print properly formatted fasta
x = open("temp.fa", "r")
y=x.readlines()
z=''.join(y)
if z[-1]=='\n': 
    z=z[:-1]
print (z)
os.remove("temp.fa") 
