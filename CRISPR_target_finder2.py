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
from Bio.Seq import Seq
from Bio import SeqFeature as sf 
from Bio.Alphabet import generic_dna

### Parse arguments
CLI=argparse.ArgumentParser()
CLI.add_argument("-target",type=str,default='user_input',help='Path to the target sequence in fasta format')
CLI.add_argument("-genome",type=str,default='user_input',help='Genome to xref for off targets')
CLI.add_argument("--init",type=str,default='NN',help='Must be two characters. e.g. GG. Dots represent anything')
CLI.add_argument("--mismatch",type=str,default='4',help='Number of maximum mismatches')
CLI.add_argument("--seed",type=int,default=1,help='Number of maximum seed mismatches')
CLI.add_argument("--distal",type=int,default=2,help='Number of maximum mismatches')
CLI.add_argument("--PAM",type=str,default='NGG',help='Number of maximum mismatches')
CLI.add_argument("--annot",type=str,help='genbank annotation file')
CLI.add_argument("--gff",type=str,help='GFF file for genome in order to determine off target significance ')
CLI.add_argument("--output",type=str,default='./',help='Output directory')


args = CLI.parse_args()
### Import arguments
#
#args.target='./target_region.fna'
#args.genome='./SpoFru_China_genome.fna'
#

if args.target=='user_input':
    target_base=list(SeqIO.parse(eg.fileopenbox(title='CHOOSE AN INPUT SEQUENCE FOR TARGETTING',default='/home/shanedenecke/Dropbox/quick_temp/CRISPR_targets/*'),'fasta'))[0]
else: 
    target_base=list(SeqIO.parse(args.target,'fasta'))[0]

if args.genome=='user_input':
    genome=eg.fileopenbox(title='CHOOSE A GENOME TO CHECK FOR OFF TARGETS',default='/home/shanedenecke/Dropbox/quick_temp/CRISPR_targets/*')
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

#### process gff file
if args.gff is not None:
    gff=pd.read_csv(args.gff,sep='\t',header=None).iloc[:,[0,2,3,4]]
    gff.columns=['chr','type','start','stop']
    
clean=pd.DataFrame()
for i in range(raw.shape[0]):
    sub=raw.iloc[i,:]
    offtar=sub['off_target']
    it=re.finditer(r"[a-z]", offtar) ## create iterator
    inds=[m.start(0) for m in it] ### get all mismatch indicies
    d=len([x for x in inds if x<10])
    sed=len([x for x in inds if x>=10])
    seqn=tardir[sub.target_seq]
    if ((sub['Chr']==tar_chr) & ((sub['loc']>tar_range[0]) & (sub['loc']<tar_range[1]))):
        tar='Good_Target'
    else:
        tar='Off_Target'
    
    
    ### attempt to use GFF file (if this exists) to add "target_feature"
    if args.gff is not None:
        loc=sub['loc']
        gffsub=gff[(gff.chr==sub['Chr']) & ((gff.start>loc-20000) & (gff.stop<loc+20000))]
        for j in range(gffsub.shape[0]):
            if gffsub.iloc[j,:].start < loc & gffsub.iloc[j,:].stop< loc:
                feature=gffsub.iloc[j,:].type
            else:
                feature='intergenic'
        ## append to overall dataframe for this target
        clean=clean.append(pd.DataFrame({'Distal':[d],'Seed':[sed],'Seq_number':[seqn],'Category':[tar],
                                     "uni_id":[sub['uni_id']],"Target_feature":[feature]}),ignore_index=True)
    else:
        clean=clean.append(pd.DataFrame({'Distal':[d],'Seed':[sed],'Seq_number':[seqn],'Category':[tar],
                                     "uni_id":[sub['uni_id']]}),ignore_index=True)

    
    
    
total=pd.merge(clean,raw,on='uni_id').sort_values(['Seq_number','Category'])

filt=total[(total['Distal']<=args.distal) & (total['Seed']<=args.seed)]


if args.annot is not None:
    with open(args.annot, "r") as handle:
        annot=list(SeqIO.parse(handle,'genbank'))[0]
    
    
    annot_frame=filt[filt.Category=='Good_Target']
    for i in range(annot_frame.shape[0]):
        sub=annot_frame.iloc[i,:]
        
        baseseq=Seq(sub.off_target)
        searchseq=str(baseseq.reverse_complement())+'|'+str(baseseq)    
        
        inds=[[m.start(0),m.end(0)] for m in re.finditer(searchseq, str(annot.seq))][0] ## find indexs of matches
        
        feat_loc = sf.FeatureLocation(inds[0],inds[1])
        my_feature = sf.SeqFeature(feat_loc,type='sgRNA',id=sub.Seq_number,strand=1)
        
        annot.features.append(my_feature)
    
    
    
    annot.seq.alphabet = generic_dna
    with open(args.output+'sgRNA_genbank_annotation.gb', 'w') as handle:
        SeqIO.write(annot,handle,'genbank')


#### add in ugene featuers???
#with open(args.output+'sgRNA_genbank_annotation.gb', 'r') as handle:
#    content = [x.strip() for x in handle.readlines()] 
    
    
#only_perfect=total[(total['Distal']==0) & (total['Seed']==0)]
#os.remove('cas_offtarget_input.txt')
#os.remove('cas_offtarget_rawoutput.txt')

filt.to_csv(sys.stdout,header=True,index=False,sep='\t')

