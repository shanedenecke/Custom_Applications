#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 20:13:49 2020

OrthoFind Output parse 

@author: shanedenecke
"""

import os
import pandas as pd
from Bio import SeqIO
import argparse
import shutil


########## Read in ARGS
CLI=argparse.ArgumentParser()
CLI.add_argument("-outdir",type=str,default='./outdir',help='Choose directory for 1:1 orthologue sequences')
CLI.add_argument("-inputdir",type=str,help='Output directory from Orhtofinder. Normal filetree')
#CLI.add_argument("-ortho_table",type=str,default='Orthogroups.tsv',help='Table output from orthofinder. Shouldbe titled Orthogroups.tsv and found in the Orthogrups folder')
#CLI.add_argument("-single_copy",type=str,default='Orthogroups_SingleCopyOrthologues.txt',help='Single copy orthologues file from Orthofinder. Should be labeled Orthogroups_SingleCopyOrthologues.txt')
CLI.add_argument("-total_fasta",type=str,default='total_fasta',help='All possible fasta sequences from the analysis in a single file')
CLI.add_argument("-maxseqs",type=int,default=1000000000,help='maximum number of sequences to retrieve')


args = CLI.parse_args()


### create ouput directory
#os.chdir('/home/shanedenecke/Dropbox/quick_temp/Ofinder_parse')
try:
    shutil.rmtree(args.outdir) ### remove directory if already exists 
except:
     pass

os.makedirs(args.outdir)


#### import single copy orthologues
with open(args.inputdir+'/Orthogroups/Orthogroups_SingleCopyOrthologues.txt') as f:
    single_ogs=f.read().splitlines()
    
### Import orthogroups
og_table=pd.read_csv(args.inputdir+'/Orthogroups/Orthogroups.tsv',sep='\t')
og_table2=og_table[og_table.Orthogroup.isin(single_ogs)]
all_ids=[]

#### get subset of all IDs to speed up process
for i in og_table2.columns:
    all_ids.append(list(og_table2.loc[:,i]))
flat_ids=[item for sublist in all_ids for item in sublist]

###import fasta parser 
recs=SeqIO.parse('total_fasta.faa','fasta')
recs_reduce=[x for x in recs if x.id in flat_ids]

single_ogs=single_ogs[0:min(args.maxseqs,len(single_ogs))]
for group in single_ogs:
    row=og_table2[og_table2.Orthogroup==group]
    og=row.iloc[0]['Orthogroup']
    minus_og=row.drop('Orthogroup',axis=1)
    idlist=list(minus_og.iloc[0])
    seqlist=[x for x in recs_reduce if x.id in idlist]
    
    SeqIO.write(seqlist,args.outdir+'/'+og+'.faa','fasta')