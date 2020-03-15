#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 20 22:32:44 2018
FASTA renamer

Argument1: Fasta file
Arugment2: Dictionary

@author: shanedenecke
"""
import os
import pandas as pd
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter
import sys 


os.path.dirname(os.path.abspath(__file__))
#os.chdir('/home/shanedenecke/fa_rename_opt')

#sys.argv=['','/home/shanedenecke/fa_rename_opt/AcrEch.fasta','/home/shanedenecke/fa_rename_opt/AcrEch_OrthoDB_key_DICT2.txt']

##define argument variables
if len(sys.argv)>2:
    df=pd.read_csv(sys.argv[2],sep=',',header=0,dtype=str)
    fa=list(SeqIO.parse(str(sys.argv[1]), 'fasta'))
    raw=SeqIO.parse(str(sys.argv[1]), 'fasta')
elif len(sys.argv)==2:
    df=pd.read_csv(sys.argv[1],sep=',',header=0,dtype=str)
    fa=list(SeqIO.parse(sys.stdin, 'fasta'))
    raw=list(SeqIO.parse(sys.stdin, 'fasta'))

#names=list(df['name'])
codes=list(df['code'])
#fa_dict=dict(zip(df['name'],df['code']))
#j=[x for x in fa if x.id=='XP_021190993.1'][0]
renamed_fa=[]

df2=dict(zip(df.code, df.name))

for i in fa:
    if i.id in df2.keys():
        i.id=df2[i.id]
        i.description=""
 
    
## Write temporary file
handle = open('temp.fa', "w")
writer = FastaWriter(handle, wrap=0)
writer.write_file(fa)
handle.close()

## Read in temporary file and print properly formatted fasta
x = open("temp.fa", "r")
y=x.readlines()
z=''.join(y)
if z[-1]=='\n': 
    z=z[:-1]
print (z)
os.remove("temp.fa") 
