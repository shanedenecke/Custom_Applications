#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 19:51:40 2020

CRISPR_Target_Finder2.py



@author: shanedenecke
"""

from Bio import SeqIO
import argparse
import re
import os
import pandas as pd
import statistics as stat
import sys
import easygui as eg

### Parse arguments
CLI=argparse.ArgumentParser()
CLI.add_argument("-target",type=str,default='user_input',help='Path to the target sequence in fasta format')
CLI.add_argument("-genome",type=str,default='user_input',help='Genome to xref for off targets')
CLI.add_argument("--init",type=str,default='NN',help='Must be two characters. e.g. GG. Dots represent anything')
CLI.add_argument("--mismatch",type=str,default='4',help='Number of maximum mismatches')
CLI.add_argument("--seed",type=int,default=1,help='Number of maximum seed mismatches')
CLI.add_argument("--distal",type=int,default=2,help='Number of maximum mismatches')
CLI.add_argument("--PAM",type=str,default='NGG',help='Number of maximum mismatches')

args = CLI.parse_args()
### Import arguments
#
#

if args.target=='user_input':
    target_base=list(SeqIO.parse(eg.fileopenbox(title='CHOOSE AN INPUT SEQUENCE FOR TARGETTING',default='/home/shanedenecke/Dropbox/*'),'fasta'))[0]
else: 
    target_base=list(SeqIO.parse(args.target,'fasta'))[0]

if args.genome=='user_input':
    genome=eg.fileopenbox(title='CHOOSE A GENOME TO CHECK FOR OFF TARGETS',default='/home/shanedenecke/Dropbox/*')
else: 
    genome=args.genome


init=args.init.replace('N','.')
mismatch=args.mismatch
pam=20*'N'+args.PAM

## create target search seqeunces
target_seq=str(target_base.seq)
revcomp=str(target_base.seq.reverse_complement())

## create target id for future naming use
target_id=target_base.id




forw=re.findall(init+'[A-Z]{19}GG', target_seq, re.I)
reve=re.findall(init+'[A-Z]{19}GG', revcomp, re.I)
allsites=forw+reve

for i in range(0,len(allsites)):
    allsites[i]=allsites[i]+' '+mismatch

input_file=[genome,pam]+allsites


with open('cas_offtarget_input.txt','w') as hand:
     hand.writelines("%s\n" % x for x in input_file)
     
     
os.system("/home/shanedenecke/Applications/Cas-offfinder/cas-offinder \
          cas_offtarget_input.txt \
          C \
          cas_offtarget_rawoutput.txt")



raw=pd.read_csv('cas_offtarget_rawoutput.txt',sep='\t',header=None)
raw.columns=['target_seq','Chr','loc','off_target','strand','mismatch']
raw['uni_id']=range(0,raw.shape[0])

offtar=list(raw['off_target'])
uniseq=list(set(raw['target_seq']))

tardir={}
for i in uniseq:
    tardir.update({i:'Sequence_'+str(uniseq.index(i))})

pure=raw[raw['mismatch']==0]
tar_chr=stat.mode(pure['Chr'])
tar_range=[stat.mean(pure['loc'])-len(target_seq),stat.mean(pure['loc'])+len(target_seq)]



clean=pd.DataFrame()
for i in range(raw.shape[0]):
    sub=raw.iloc[i,:]
    offtar=sub['off_target']
    it=re.finditer(r"[a-z]", offtar) ## create iterator
    inds=[m.start(0) for m in it] ### get all mismatch indicies
    d=len([x for x in inds if x<10])
    sed=len([x for x in inds if x>=10])
    seqn=tardir[sub.target_seq]
    #dist.append(len([x for x in inds if x<10]))
    #seed.append(len([x for x in inds if x>=10]))
    #seqnum.append(tardir[sub.target_seq])
    if ((sub['Chr']==tar_chr) & ((sub['loc']>tar_range[0]) & (sub['loc']<tar_range[1]))):
        tar='Good_Target'
    else:
        tar='Off_Target'
    clean=clean.append(pd.DataFrame({'Distal':[d],'Seed':[sed],'Seq_number':[seqn],'Category':[tar],"uni_id":[sub['uni_id']]}),ignore_index=True)

total=pd.merge(clean,raw,on='uni_id').sort_values(['Seq_number','Category'])

filt=total[(total['Distal']<=args.distal) & (total['Seed']<=args.seed)]
#only_perfect=total[(total['Distal']==0) & (total['Seed']==0)]
os.remove('cas_offtarget_input.txt')
os.remove('cas_offtarget_rawoutput.txt')

filt.to_csv(sys.stdout,header=True,index=False,sep='\t')

