#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 11:08:02 2020

@author: shanedenecke
"""

from Bio import SeqIO
from Bio import pairwise2
import os
import subprocess 
import pandas as pd
#from Bio.SubsMat.MatrixInfo import blosum62
#from Bio.SubsMat.MatrixInfo import mclach
from Bio.SubsMat.MatrixInfo import levin
import argparse


#os.chdir('/home/shanedenecke/Dropbox/quick_temp')


########## Read in ARGS
CLI=argparse.ArgumentParser()
CLI.add_argument("-start_organism",type=str,default='/home/shanedenecke/Dropbox/omics_projects/omics_ref/Nezara_viridula/unigene/NezVir_unigene.faa',help='')
CLI.add_argument("-target_organism",type=str,default='/home/shanedenecke/Dropbox/omics_projects/omics_ref/Homo_sapiens/unigene/HomSap_unigene.faa',help='')
CLI.add_argument("--output",type=str,default='./pairwise_output.tsv',help='')
CLI.add_argument("--threads",type=str,default='1',help='number of threads')

args = CLI.parse_args()

#args.start_organism='/home/shanedenecke/Dropbox/quick_temp/test_start.faa'

### Import proteomes
start_proteome=list(SeqIO.parse(args.start_organism,format='fasta'))
target_proteome=list(SeqIO.parse(args.target_organism,format='fasta'))
#start_name=os.path.basename(args.start_organism)
#target_name=os.path.basename(args.target_organism)

### Run blast
with subprocess.Popen(['blastp','-query',args.start_organism,'-db',args.target_organism,
                           '-outfmt','6 qseqid sseqid evalue qcovs pident','-num_threads',args.threads],stdout=subprocess.PIPE) as proc:
        blast_out=pd.read_csv(proc.stdout,sep='\t',header=None)
       

blast_out.to_csv('Blast_output.csv',index=False)
#blast_out=pd.read_csv('Blastout.csv',header=None)
### Filter blast and target _proteome
uniblast=blast_out.drop_duplicates(subset=0,keep='first')
start_ids=[x.id for x in start_proteome]
target_proteome_lean=[x for x in target_proteome if x.id in list(uniblast.iloc[:,1])]


d={}
for i in start_ids[0:500]:
    if i not in list(uniblast.iloc[:,0]):
        d.update({i:'No Match'})
        pass
    else:
        row=uniblast[uniblast[0]==i]
        target_seq=[str(x.seq) for x in target_proteome_lean if x.id==row.iloc[:,1].values[0]][0]
        start_seq=[str(x.seq) for x in start_proteome if x.id==i][0]
        #alscore = pairwise2.align.globaldx(target_seq, start_seq, blosum62,score_only=True)
        align = pairwise2.align.globalxx(target_seq, start_seq)
        score=align[-1][-3]
        l=align[-1][-1]
        
        #alscore=pairwise2.align.globalxx(target_seq, start_seq,score_only=True)/len(target_seq)
        #corrected=alscore/(len(start_seq)*2)
        corrected=score/l
        d.update({i:corrected})


final_out=pd.DataFrame(list(d.items()))
final_out.columns=['Gene_id','Identity_score'] 
final_out.to_csv(args.output,sep='\t',index=False)
#os.system('cat '+args.output)
