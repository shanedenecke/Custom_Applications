#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Sat May 18 17:20:39 2019

@author: shanedenecke
"""

import os
from Bio import SeqIO
import sys 
import argparse

os.path.dirname(os.path.abspath(__file__))
#os.chdir('/home/shanedenecke/Dropbox/Misc/random/unigene_py_test')

#sys.argv=['','HomSap_unigene.faa','HomSap_SLC_genes.txt']


CLI=argparse.ArgumentParser()
CLI.add_argument("-fasta",type=str,default='HomSap_unigene.faa',help='Add fasta file')
CLI.add_argument("-unigenes",type=str,default='HomSap_SLC_genes.txt',help='Add unigene list')
args = CLI.parse_args()



### standard input
#if not sys.stdin.isatty():
#    args = CLI.parse_args()
#    codes=[line.rstrip('\n') for line in sys.stdin] 
#elif sys.stdin.isatty():
#    args = CLI.parse_args()
#    CLI.add_argument("-unigenes",type=str,default='HomSap_SLC_genes.txt',help='Add unigene list')
#    codes=[line.rstrip('\n') for line in open(args.unigenes)]


### import data
recs=SeqIO.to_dict(SeqIO.parse(args.fasta, 'fasta'),key_function=lambda rec: rec.description)
codes=[line.rstrip('\n') for line in open(args.unigenes)]

#print(codes)

##define argument variables
#if len(sys.argv)>2:
#    recs=SeqIO.to_dict(SeqIO.parse(sys.argv[1], 'fasta'),key_function=lambda rec: rec.description)
#    codes=[line.rstrip('\n') for line in open(sys.argv[2])]

#if len(sys.argv)==2:
#    recs=SeqIO.to_dict(SeqIO.parse(sys.stdin, 'fasta'))
#    codes=[line.rstrip('\n') for line in open(sys.argv[1])]




final_dict={}
for code in codes:   
    uni={k:v for (k,v) in recs.items() if code in k}
    if len(uni)==0:
        pass
    elif len(uni)==1:
        final_dict.update(uni)
    elif len(uni)>1:
        mat={k:len(v.seq) for (k,v) in recs.items() if code in k}
        maximum=sorted(mat.values())[-1]
        f={k:v for k,v in mat.items() if maximum==v}
        indkey=list(f.keys())[0]
        final={indkey:recs[indkey]}
        final_dict.update(final)
        

a=list(final_dict.values())
with open('temp.fa','w') as f:
    for i in a:
        f.write('>'+i.description+'\n'+str(i.seq)+'\n')


x = open("temp.fa", "r")
y=x.readlines()
z=''.join(y)
if z[-1]=='\n': 
    z=z[:-1]
print (z)
os.remove("temp.fa") 




        

        
## Write temporary file


## Read in temporary file and print properly formatted fasta


