#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 15:47:34 2020

@author: shanedenecke
"""

#import os
import pandas as pd
import re
import click


@click.command()
@click.option('-primer3', help='Table output from primer 3 command line program. should be list of names=values')
@click.option('-seqname', default='test_sequence', help='Sequence name to be used. Helps during output')
@click.option('-outfile', default='./Primer3_parse_output.tsv', help='Name of the tsv file to output')
@click.option('--number', default='2', help='How many primer pairs to take output for each sequence')
#os.chdir('/home/shanedenecke/Dropbox/omics_projects/aw1_transcriptome/Spodoptera_marker_sequences_Kathrin')

def primer3_parse(primer3, seqname, outfile,number):
    prinum=str(int(number)-1)
    raw_primer3 = [line.rstrip() for line in open(primer3) if re.search('_[0-'+prinum+']',line)]
    
    final=pd.DataFrame()
    for i in range(0,5):
        #i=0
        subset=[x for x in raw_primer3 if re.search('_'+str(i),x)]
        subset2=[x.replace('_'+str(i),'') for x in subset]
        #d={x.split('=')[0]:[x.split('=')[1],'primer_'+str(i)] for x in subset2}
        d={x.split('=')[0]:[x.split('=')[1]] for x in subset2}
        frame=pd.DataFrame(d)
        frame['Primer_number']='Primer_'+str(i)
        final=final.append(frame)
    
    final['Sequence_name']=seqname
    simple=final[['Sequence_name','Primer_number','PRIMER_LEFT_SEQUENCE','PRIMER_RIGHT_SEQUENCE','PRIMER_LEFT_TM','PRIMER_RIGHT_TM','PRIMER_PAIR_PRODUCT_SIZE']]
    simple.to_csv(outfile,index=None,sep='\t')

if __name__ == '__main__':
    primer3_parse()      
    



#d={x.split('=')[0]:x.split('=')[1] for x in subset2}
#frame=pd.DataFrame(d.items())
#frame.columns=['feature','value']
#frame['primer_number']='primer_'+str(i)
#l.append(frame)
#final=final.append(frame)


