#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 20:13:49 2020

OrthoFind Output parse 

################### Error at line 57 on server. Also need to add in some lines at the beginning to arrange output in more general format (avoid date folder)


@author: shanedenecke
"""

import os
import pandas as pd
from Bio import SeqIO
import argparse
import shutil


########## Read in ARGS
CLI=argparse.ArgumentParser()
CLI.add_argument("-outdir",type=str,default='./Hemispec/one_to_one',help='Choose directory for 1:1 orthologue sequences')
CLI.add_argument("-inputdir",type=str,default='./Hemispec/',help='Output directory from Orhtofinder. Normal filetree')
CLI.add_argument("-total_fasta",type=str,default='./Hemispec/tempseqs/total_proteome.faa',help='All possible fasta sequences from the analysis in a single file')
CLI.add_argument("-maxseqs",type=int,default=1000000000,help='maximum number of sequences to retrieve')

args = CLI.parse_args()


### create ouput directory
#os.chdir('')
try:
    shutil.rmtree(args.outdir) ### remove directory if already exists 
except:
     pass

os.makedirs(args.outdir)
os.system('cp '+args.inputdir+'/*/ '+args.inputdir+'/') ### extract contents



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

###Concatonate fasta files if needed. Otherwise just import file
if os.path.isdir(args.total_fasta):
    recs=[]
    filenames = [x for x in os.listdir(args.total_fasta) if '.fa' in x]
    for i in filenames:
        seqIO.parse(i,'fasta')
    
    
    with open(args.total_fasta+'/TOTAL_FASTA.faa', 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                for line in infile:
                   outfile.write(line)
else:
    recs=SeqIO.parse(args.total_fasta,'fasta')
    recs_reduce=[x for x in recs if x.id in flat_ids]


single_ogs=single_ogs[0:min(args.maxseqs,len(single_ogs))]
for group in single_ogs:
    row=og_table2[og_table2.Orthogroup==group]
    og=row.iloc[0]['Orthogroup']
    minus_og=row.drop('Orthogroup',axis=1)
    idlist=list(minus_og.iloc[0])
    seqlist=[x for x in recs_reduce if x.id in idlist]
    
    SeqIO.write(seqlist,args.outdir+'/'+og+'.faa','fasta')