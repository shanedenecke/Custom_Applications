#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 10:02:53 2020

Script to subset fasta into a unigene format based on GFF file

@author: shanedenecke
"""

########## Import libraries ##########
import os
from Bio import SeqIO
import argparse
import re
import pandas as pd
from tqdm import tqdm


########## Set working direcoty ##########
os.path.dirname(os.path.abspath(__file__))
#os.chdir('/home/shanedenecke/Dropbox/quick_temp/spofru_test')


########## Add in arguments ##########
CLI=argparse.ArgumentParser()
CLI.add_argument("-fasta",type=str,default='./GCF_011064685.1_ZJU_Sfru_1.0_protein.faa',help='Fasta file downloaded from NCBI')
CLI.add_argument("-gff",type=str,default='./GCF_011064685.1_ZJU_Sfru_1.0_genomic.gff',help='GFF file downloaded from NCBI genome page')
#CLI.add_argument("-genName",type=str,default='gene',help='Field identifier for gene in GFF file')
#CLI.add_argument("-outName",type=str,default='unigene.faa',help='name of output file')


args = CLI.parse_args()


########## Import raw files ##########
gff_raw=pd.read_csv(args.gff,sep='\t',comment='#',header=None)
recs=list(SeqIO.parse(args.fasta, 'fasta'))
recs=SeqIO.to_dict(SeqIO.parse(args.fasta, 'fasta'),key_function=lambda rec: rec.id)




########## Process GFF into key with gene protein length
gff_sub=gff_raw[gff_raw[2]=='CDS'][[2,3,4,8]]
gff_sub[8]=[re.sub('ID=cds-([X|Y]P_[0-9]+..).+gene=(LOC[0-9]+).+$','\\1__\\2',x) for x in gff_sub[8]]
gff_sub[['protein','gene']]=gff_sub[8].str.split('__',expand=True)
gff_sub=gff_sub[~gff_sub['gene'].isna()]
gff_sub['length']=abs(gff_sub[3]-gff_sub[4])

idkey=gff_sub[['protein','gene']]
lenkey=gff_sub[['protein','length']]

########## Find longest isoform ##########
lenkey2=lenkey.groupby('protein').sum()
lenkey2['protein']=lenkey2.index
lenkey2=lenkey2.reset_index(drop=True)

########## merge lenth and id keys to get unique protein list ##########
fullKey=pd.merge(lenkey2,idkey,on='protein')
longest=fullKey.loc[fullKey.groupby('gene')['length'].idxmax()]
longestDict=pd.Series(longest.gene.values,index=longest.protein).to_dict()
unicodes=list(longest['protein'])


########## Subset Fasta file ##########

for i in tqdm(unicodes):
    temp=recs[i]
    temp.id=i+'__'+longestDict[i]
    with open('temp.fa','a') as f:
                f.write('>'+temp.id+'\n'+str(temp.seq)+'\n')

x = open("temp.fa", "r")
open("temp.fa", "r")
y=x.readlines()
z=''.join(y)
if z[-1]=='\n': 
    z=z[:-1]
print (z)
os.remove("temp.fa") 
