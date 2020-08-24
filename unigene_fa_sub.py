#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Sat May 18 17:20:39 2019

@author: shanedenecke
"""

import os
from Bio import SeqIO
import argparse
import re
import sys
from tqdm import tqdm

os.path.dirname(os.path.abspath(__file__))
#os.chdir('/home/shanedenecke/Dropbox/quick_temp/abc_compare')



CLI=argparse.ArgumentParser()
CLI.add_argument("-fasta",type=str,default=(None if sys.stdin.isatty() else sys.stdin),help='Add fasta file')
CLI.add_argument("-mode",type=str,default='regex',help='Either "regex" or "file". Select regex if you want to filter the fasta all in one command. file takes a list of unique identifiers and uses that for filtering')
CLI.add_argument("-codes",type=str,help='Add unigene list. Can either be regex string or file of ids depending on mode')
CLI.add_argument("-outfmt",type=str,default='Full',help='"Full" or "Short". Do you want the final output to be the full description or only the ID')


args = CLI.parse_args()


#args.fasta='pxug_v01.faa'
#args.mode='regex'
#args.codes='(PXUG_V[0-9]_[0-9]+).+$'
#args.outfmt='Short'

### import data
recs=SeqIO.to_dict(SeqIO.parse(args.fasta, 'fasta'),key_function=lambda rec: rec.description)


### initialize final dictionary which will contain all outputs
final_dict={}


### obtain unicodes variable either from file or from regular expression depending on mode
if args.mode=='regex':
    codes=[re.sub(args.codes,'\\1',x) for x in recs.keys()]
    unicodes=set(codes)
else:
    unicodes=[line.rstrip('\n') for line in open(args.codes)]


### iterate over unicodes and subset out longest isoform
for code in tqdm(unicodes):
    
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
        if args.outfmt=='Full':
            with open('temp.fa','a') as f:
                f.write('>'+final.description+'\n'+str(final.seq)+'\n')
                #final_dict.update({final.description:final.seq})
        elif args.outfmt=='Short':
            with open('temp.fa','a') as f:
                f.write('>'+final.id+'\n'+str(final.seq)+'\n')
            #final_dict.update({final.id:final.seq})        


### Write to temporary file
#with open('temp.fa','w') as f:
#    for k,v in final_dict.items():
#        f.write('>'+k+'\n'+str(v)+'\n')

### Print file for output
x = open("temp.fa", "r")
y=x.readlines()
z=''.join(y)
if z[-1]=='\n': 
    z=z[:-1]
print (z)
os.remove("temp.fa") 



