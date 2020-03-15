#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 15 15:29:00 2020

Parse PFAM domain search from HMM and extract particular domain from output.
e.g. hmmsearch --domtblout ABC

@author: shanedenecke
"""


### Packages
import os
import pandas as pd
import argparse
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter
from pathlib import Path


#### set home directory
#os.chdir('/home/shanedenecke/Dropbox/quick_temp')
home = str(Path.home())



################################################## Read in ARGS
CLI=argparse.ArgumentParser()
CLI.add_argument("-table",type=str,default='test3.txt',help='Path to file containing hmmsearch output in --domtable format')
CLI.add_argument("-fasta",type=str,default='SpoFru_China_unigene.faa',help='Fasta that you want to subset from. Can be entire proteome or subset list of genes')
args = CLI.parse_args()



#### Import basic data
hmmsearch=pd.read_csv(str(args.table),comment='#',header=None,sep='\s+',usecols=[0,19,20])
hmmsearch.columns=['geneid','start','stop']
base_fa=list(SeqIO.parse(args.fasta,'fasta'))
#recs=SeqIO.to_dict(SeqIO.parse(str(args.fasta), 'fasta'),key_function=lambda rec: rec.description)
unigenes=list(set(hmmsearch['geneid']))
sub_fa=[x for x in base_fa if x.id in unigenes]


final_fa=[]
###parse table to include only first instance of domain
for i in unigenes:
    sub_base=hmmsearch[hmmsearch['geneid']==i]
    minstart=min(sub_base['start'])
    sub_real=sub_base[sub_base['start']==minstart].iloc[:]
    
    
    target_fa=[x for x in sub_fa if x.id==sub_real.iloc[0,:]['geneid']][0]
    target_fa.seq=target_fa.seq[int(sub_real['start']):int(sub_real['stop'])]
    
    final_fa.append(target_fa)
    

## write temporary file
handle = open('temp.fa', "w")
writer = FastaWriter(handle, wrap=0)
writer.write_file(final_fa)
handle.close()


x = open("temp.fa", "r")
y=x.readlines()
z=''.join(y)

if z[-1]=='\n': 
    z=z[:-1]
print (z)
os.remove("temp.fa") 