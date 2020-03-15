#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 16 14:19:12 2018

Argument 1/Stdin Raw IP scan output ideally from whole genome
Argument 2 Fasta of proteins to be subsetted
Arghument 3 Domain that you are looking for

@author: shanedenecke
"""
##import libraries
import os
import pandas as pd
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter
import sys 

##Set Directory
os.path.dirname(os.path.abspath(__file__))
#os.chdir('/data2/shane/Transporter_ID/ABC_id/')

#sys.argv=['','./Filter/Unsorted_clean/ABC_unsorted_total_length.faa','./Filter/IPSCAN.tsv','PF00005']


##define argument variables
if len(sys.argv)>3:
    coord=pd.read_csv(sys.argv[2],sep='\t',header=None,usecols=[0,2,3,4,6,7],names = ["gene", "len","database", "code","start", "end"])
    prot_list=list(SeqIO.parse(str(sys.argv[1]), 'fasta'))
    domain=str(sys.argv[3])
elif len(sys.argv)==3:
    coord=pd.read_csv(str(sys.argv[1]),sep='\t',header=None,usecols=[0,2,3,4,6,7],names = ["gene", "len","database", "code","start", "end"])
    prot_list=list(SeqIO.parse(sys.stdin, 'fasta'))
    domain=str(sys.argv[2])
  

## Filter IP Scan Output
abc=[x.id for x in prot_list]
ip_gene=coord[coord.gene.isin(abc)]
ip_fam=ip_gene[ip_gene['code']==domain]

#a=[x for x in domain_fa if x.id=='MicDem__ABC_Unsorted___LOC103573255']
## subset out domain of interest only taking first instance of domain
domain_fa=prot_list[:]
delist=[]
prot_names=[x.id for x in domain_fa]
for s in domain_fa:
    ind=prot_names.index(s.id)
    gene=s.id
    #sub1=ip_fam[ip_fam['gene']==i].sort_values(by='start').head(1) ## subset and sort ip data for each gene
    if gene in list(ip_fam['gene']):
        start=int(min(ip_fam.loc[ip_fam['gene']==gene]['start']))-1
        end=int(min(ip_fam.loc[ip_fam['gene']==gene]['end']))-1
        s.seq=s.seq[start:end]
    else:
        delist.append(s.id)
       

final=[x for x in domain_fa if x.id not in delist]

## write temporary file
handle = open('temp.fa', "w")
writer = FastaWriter(handle, wrap=0)
writer.write_file(final)
handle.close()


x = open("temp.fa", "r")
y=x.readlines()
z=''.join(y)

if z[-1]=='\n': 
    z=z[:-1]
print (z)
os.remove("temp.fa") 

