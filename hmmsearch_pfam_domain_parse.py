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
#os.chdir('/home/sdenecke/Transporter_ID/ABC_id')
home = str(Path.home())



################################################## Read in ARGS
CLI=argparse.ArgumentParser()
CLI.add_argument("-table",type=str,default='./Filter/HMM_PF00005_output.tsv',help='Path to file containing hmmsearch output in --domtable format')
CLI.add_argument("-fasta",type=str,default='./Filter/ABC_preliminary_total.faa',help='Fasta that you want to subset from. Can be entire proteome or subset list of genes')
args = CLI.parse_args()



#### Import basic data
hmmsearch=pd.read_csv(str(args.table),comment='#',header=None,sep='\s+',usecols=[0,19,20])
hmmsearch.columns=['geneid','start','stop']
base_fa=list(SeqIO.parse(args.fasta,'fasta'))
recs=SeqIO.to_dict(SeqIO.parse(str(args.fasta), 'fasta'),key_function=lambda rec: rec.id)
unigenes=list(set(hmmsearch['geneid']))
sub_fa=[x for x in base_fa if x.id in unigenes]
sub_fa_names=[x.id for x in sub_fa]

final_fa=[]
###parse table to include only first instance of domain
for i in unigenes:
    sub_base=hmmsearch[hmmsearch['geneid']==i]
    #nam=str([x for x in sub_fa_names if x==i])
    
    try:
        target_fa=recs[i]
    except KeyError:
        #print(i+' Not in fasta file')
        continue
    
    if sub_base.shape[0]==0:
        #print(i+' Has no Relevant domain')
        continue
    elif sub_base.shape[0]>0:
        #target_fa=target_fa[0]
        minstart=min(sub_base['start'])
        sub_real=sub_base[sub_base['start']==minstart].iloc[:]
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
