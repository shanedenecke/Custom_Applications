#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 16:22:36 2020

@author: shanedenecke
"""

import subprocess as sp
import pandas as pd
from Bio import SeqIO
import argparse
import os
import re

#os.chdir('/home/shanedenecke/Dropbox/omics_projects/aw1_transcriptome/Spodoptera_marker_sequences_Kathrin')

CLI=argparse.ArgumentParser()
CLI.add_argument("-input",type=str,help='List of sequences that you start with')
CLI.add_argument("-target",type=str,help='Proteome of species you want to target against')
CLI.add_argument("-seqs",type=str,help='List of unique identifiers for genes that you are interested in')
CLI.add_argument("--qcov",type=int,default=30,help='Threshold for query cover in blast')
CLI.add_argument("--evalue",type=int,default=1e-3,help='Threshold for evalue in blast')
CLI.add_argument("--threads",type=int,default=2,help='Number of threads')
CLI.add_argument("--outdir",type=str,default='Recip_blast_output',help='name of output directory')
CLI.add_argument("--blast_type",type=str,default='blastp',help='Type of Blast that is performed. Can be blastx, blastp, blastx, or tblastn')
CLI.add_argument("--target_spname",type=str,default='Species1',help='What is the name of the target species (column outputs)')
CLI.add_argument("--input_spname",type=str,default='Species2',help='What is the name of the starting species (column outputs)')


args = CLI.parse_args()

#os.chdir('/mnt/disk/shane/temp/raw_proteomes')
#args.target='./DroMel_unigene.faa'
#args.input='./SpoFru_unigene.faa'
#args.seqs='./SpoFru_ids.txt'
#args.blast_type='blastp'
#args.outdir='SpoFru_DroMel_recip'

### parse blast type
if args.blast_type=='blastp':
    blast1='blastp'
    blast2='blastp'
    db1='prot'
    db2='prot'
elif args.blast_type=='blastx':
    blast1='blastx'
    blast2='tblastn'
    db1='nucl'
    db2='prot'
elif args.blast_type=='blastn':
    blast1='blastn'
    blast2='blastn'
    db1='nucl'
    db2='nucl'
elif args.blast_type=='tblastn':
    blast1='tblastn'
    blast2='blastx'
    db1='prot'
    db2='nucl'

#print(blast1)
#print(blast2)
#print(db1)
#print(db2)

try:
    os.makedirs(args.outdir)
except:
    print('Directory Exists')

### make blast databases
input_dir=os.path.dirname(args.input)
input_base=os.path.basename(args.input)#.split('.')[0:-1] #base=''.join(base) ### extra code for collapsing list
if len([x for x in os.listdir(input_dir) if re.search(input_base+'.psq',x)])==0:
    sp.run(['makeblastdb','-in',args.input,'-parse_seqids','-dbtype',db1])

target_dir=os.path.dirname(args.target)
target_base=os.path.basename(args.target)
if len([x for x in os.listdir(target_dir) if re.search(target_base+'.psq',x)])==0:
    print('test')
    sp.run(['makeblastdb','-in',args.target,'-parse_seqids','-dbtype',db2])


### Subset out initial sequences
     
### read in sequence IDS
with open(args.seqs) as f:
    input_ids=[line.rstrip() for line in f]


### import data
recs=SeqIO.to_dict(SeqIO.parse(args.input, 'fasta'),key_function=lambda rec: rec.description)


### initialize final dictionary which will contain all outputs
final_dict={}


### iterate over unicodes and subset out longest isoform
for code in input_ids:
    
    ### remove starting > from id and find all instances of seqs in fasta
    code=re.sub('>','',code)
    uni={k:v for (k,v) in recs.items() if code in k}
    
    #### skip code if not found in fasta
    if len(uni)==0:
        pass
    
    #### If found  then take longest isoform
    elif len(uni)>0: 
        mat={k:len(v.seq) for (k,v) in recs.items() if code in k} ### make dictionary with length of each seqeunce
        maximum=sorted(mat.values())[-1] ### take longest one. If multiple will take random sequence
        f={k:v for k,v in mat.items() if maximum==v} ###
        indkey=list(f.keys())[0]
        final=recs[indkey]
        final_dict.update({final.id:final.seq})        


### Write to temporary file
with open(args.outdir+'/input_sequences.faa','w') as f:
    for k,v in final_dict.items():
        f.write('>'+k+'\n'+str(v)+'\n')


  
### Run initial blast
cmd=[blast1,'-query',args.outdir+'/input_sequences.faa','-db',args.target,'-qcov_hsp_perc',str(args.qcov),'-evalue',str(args.evalue),
                           '-outfmt','6 qseqid sseqid evalue qcovs pident','-num_threads',str(args.threads),'-max_target_seqs','4']

#with open("blast_error_output.txt", "w") as f:
#    try:
#        sp.check_call(cmd, stderr=f)
#        with sp.Popen(cmd,stdout=sp.PIPE) as proc:
#            initial_blast=pd.read_csv(proc.stdout,sep='\t',header=None)
#    except sp.CalledProcessError as e:
#        print(e)
#        exit(1)

### Run reciprocal blast
with sp.Popen(cmd,stdout=sp.PIPE) as proc:
        initial_blast=pd.read_csv(proc.stdout,sep='\t',header=None)



initial_blast.columns=['qseqid', 'sseqid', 'evalue', 'qcovs', 'pident']
initial_blast.to_csv(args.outdir+'/Initial_blast_raw.tsv')
initial_blast.loc[:,'sseqid']=[x.replace('emb|','').replace('|','') for x in initial_blast.loc[:,'sseqid']] ### Remove stupid emb| crap from sseqid

#### Subset out best hits from target proteome
initial_blast_final=initial_blast.sort_values(by='evalue').drop_duplicates(subset='qseqid',keep='first')
target_ids=initial_blast_final.loc[:,'sseqid'].values.tolist() #### chains a bunch of commands. Sorts by evalue drops dupicates and extracts sseqids
target_seqs=[x for x in SeqIO.parse(args.target,'fasta') if x.id in target_ids] 
SeqIO.write(target_seqs,args.outdir+'/Target_best_hits.faa','fasta') 


cmd2=[blast2,'-query',args.outdir+'/Target_best_hits.faa','-db',args.input,'-qcov_hsp_perc',str(args.qcov),'-evalue',str(args.evalue),
                           '-outfmt','6 qseqid sseqid evalue qcovs pident','-num_threads',str(args.threads),'-max_target_seqs','4']
### Run reciprocal blast
with sp.Popen(cmd2,stdout=sp.PIPE) as proc:
        recip_blast=pd.read_csv(proc.stdout,sep='\t',header=None)

recip_blast.columns=['qseqid', 'sseqid', 'evalue', 'qcovs', 'pident']
recip_blast.to_csv(args.outdir+'/Reciprocal_blast_raw.tsv')
recip_blast.loc[:,'sseqid']=[x.replace('emb|','').replace('|','') for x in recip_blast.loc[:,'sseqid']] ### Remove stupid emb| crap from sseqid

#### Subset out best hits from initial proteome
recip_blast_final=recip_blast.sort_values(by='evalue').drop_duplicates(subset='qseqid',keep='first')


### compare two blast_outputs
ind=dict(zip(initial_blast_final.qseqid,initial_blast_final.sseqid))
red=dict(zip(recip_blast_final.qseqid,recip_blast_final.sseqid))

final_d={}
for k,v in ind.items():
    #initial_id=k
    #initial_hit=v
    #recip_hit=red[v]
    #recip_hit=[k for k,v in ind.items() if v==recip_value][0]
    
   
    if k not in red.values():
        print('no reciprocal match for '+k)
    elif v not in red.keys():
    	print('reciprocal mismatch for '+k)
    elif k==red[v]:
        final_d.update({k:v})
    else:
         print('reciprocal mismatch for '+k)

      
final_frame=pd.DataFrame(list(final_d.items()))    
final_frame.columns=[args.input_spname,args.target_spname]
final_frame.to_csv(args.outdir+'/Final_reciprocal_best_hits.tsv',index=False)

try:
    os.remove(args.outdir+'/Final_reciprocal_hits.fasta')
except:
    pass

for i in range(final_frame.shape[0]):
    row=final_frame.iloc[i,:]
    start=row[args.input_spname]
    target=row[args.target_spname]
    fa_raw=[x for x in SeqIO.parse(args.outdir+'/Target_best_hits.faa','fasta') if x.id==target][0]
    fa_raw.id=args.target_spname+'__'+target+'__'+args.input_spname+'__'+start
    with open(args.outdir+'/Final_reciprocal_hits.fasta','a') as f:
        f.write('>'+fa_raw.id+'\n')
        f.write(str(fa_raw.seq)+'\n')





