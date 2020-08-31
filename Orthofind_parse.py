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
CLI.add_argument("-outdir",type=str,help='Choose directory for 1:1 orthologue sequences')
CLI.add_argument("-indir",type=str,help='Output directory from Orhtofinder. Top level of Orthofinder outpu')
CLI.add_argument("-total_fasta",type=str,help='All possible fasta sequences from the analysis in a single file')
CLI.add_argument("-maxseqs",type=int,default=1000000000,help='maximum number of sequences to retrieve')
CLI.add_argument("-mode",type=str,default='id',help='Either seq or id. Do you want seqeunces to be returned or a table with IDs?')

args = CLI.parse_args()




### create ouput directory
#os.chdir('/mnt/disk/shane/Transporter_ID/Arthropod_ABC_pipeline/')
#args.outdir='./CAFE/Hemimetabola_taxid_codesc/one_to_one'
#args.indir='./CAFE/Hemimetabola_taxid_codes/orthofinder_temp'
#args.total_fasta='./CAFE/Hemimetabola_taxid_codes/tempseqs/'
#args.maxseqs=100000
#args.mode='seq'

try:
    shutil.rmtree(args.outdir) ### remove directory if already exists 
except:
     pass

os.makedirs(args.outdir)
os.system('cp -r '+args.indir+'/*/Orthogroups '+args.indir+'/') ### extract contents



#### import single copy orthologues
with open(args.indir+'/Orthogroups/Orthogroups_SingleCopyOrthologues.txt') as f:
    single_ogs=f.read().splitlines()
    
### Import orthogroups
og_table=pd.read_csv(args.indir+'/Orthogroups/Orthogroups.tsv',sep='\t')
og_table2=og_table[og_table.Orthogroup.isin(single_ogs)]
if 'NezVir_unigene' in og_table2.columns: og_table2.NezVir_unigene=[x.split(' ')[0] for x in og_table2.NezVir_unigene]
all_ids=[]


if args.mode=='seq': ############# BROKEN
    #### get subset of all IDs to speed up process
    for i in og_table2.columns:
        all_ids.append(list(og_table2.loc[:,i]))
    flat_ids=[item for sublist in all_ids for item in sublist]
    #flat_ids=[x.split(' ')[0] for x in flat_ids]
    ###Concatonate fasta files if needed. Otherwise just import file
    if os.path.isdir(args.total_fasta):
        recs_reduce=[]
        filenames = [args.total_fasta+'/'+x for x in os.listdir(args.total_fasta) if '.fa' in x]
        for x in filenames:
            species=os.path.basename(x).replace('_unigene.faa','')
            proteome=list(SeqIO.parse(x,'fasta'))
            for seq in proteome:
                if seq.id in flat_ids: 
                    seq.description=species
                    recs_reduce.append(seq)
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
        for i in range(0,len(seqlist)):
            seqlist[i].id=seqlist[i].description
        
        SeqIO.write(seqlist,args.outdir+'/'+og+'.faa','fasta')

elif args.mode=='id':
    og_table2.to_csv(args.outdir+'/one_to_one_orthology_table.tsv',sep='\t',index=False)

else:
    print('Please select a mode!')
