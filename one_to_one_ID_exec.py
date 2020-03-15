#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 18 17:20:39 2019

@author: shanedenecke
"""

### Packages
import os
import pandas as pd
import argparse
import warnings
#import sys
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter
import gc 
import re
import shutil 
from pathlib import Path

home = str(Path.home())


warnings.filterwarnings("ignore")
def unicount(series):
    return len(list(set(series)))

def database_find(taxid,data):
    x=set(list(data[data['taxid']==taxid].loc[:,'database']))
    if 'ENSEMBL' in x:
        return('ENSEMBL')
    elif'NCBIgid' in x:
        return('NCBIgid')
    elif'UniProt' in x:
        return('UniProt')
    elif'NCBIproteinGI' in x:
        return('NCBIproteinGI')
    elif'InterPro' in x:
        return('InterPro')
    elif'Trinity' in x:
        return('Trinity')
    else:
        print('NO MATCH')


## Set working directory
#os.path.dirname(os.path.abspath(__file__))
#os.chdir('/home/shanedenecke/test')




################## ARGUMENT DEBUG

#sys.argv=['','Nezara',['85310_0','7070_0'],'id']
#if(sys.argv[1]=='Arthropod' or sys.argv[1]=='arthropod'):
#    odb_node=pd.read_csv('/home/shanedenecke/Applications/custom/OrthoDB/odb10v0_OG2genes.6656.tab',sep='\t',header=None)  
#if(sys.argv[1]=='Metazoa' or sys.argv[1]=='metazoa'):
#    odb_node=pd.read_csv('/home/shanedenecke/Applications/custom/OrthoDB/odb10v0_OG2genes.33208.tab',sep='\t',header=None)
#if(sys.argv[1]=='Nezara' or sys.argv[1]=='nezara'):
    #odb_node=pd.read_csv('/home/shanedenecke/Applications/custom/OrthoDB/Nezara_node_Arth_new.tsv',sep='\t',header=None)
#if len(sys.argv[2])==1:
#    taxid_list=[line.rstrip('\n') for line in open(sys.argv[2][0])]
#else:
#    taxid_list=sys.argv[2]
#arg3=sys.argv[3]


################################################## Read in ARGS
CLI=argparse.ArgumentParser()
CLI.add_argument("-node",type=str,default='Metazoa',help='Choose node to cluster species at. Can be either "Metazoa or Arthropod"')
CLI.add_argument("-taxid",nargs="*",type=str,default=['7227_0','9606_0'],help='put in space separate list of taxid. Should be in the format NUMBER_0. e.g. "7227_0" for Drosophila')
CLI.add_argument("-output",type=str,default='id',help='Choose the output type. For a table of IDs put "id". For sequences put "seq"')
args = CLI.parse_args()




####################################3




if(args.node=='Arthropod' or args.node=='arthropod'):
    odb_node=pd.read_csv(home+'/Applications/Custom_Applications/OrthoDB_source/odb10v0_OG2genes.6656.tab',sep='\t',header=None)    
if(args.node=='Metazoa' or args.node=='metazoa'):
    odb_node=pd.read_csv(home+'/Applications/Custom_Applications/OrthoDB_source/odb10v0_OG2genes.33208.tab',sep='\t',header=None)
if len(args.taxid)==1:
    taxid_list=[line.rstrip('\n') for line in open(args.taxid[0])]
else:
    taxid_list=args.taxid
arg3=args.output


###Key file
odb_key=pd.read_csv(home+'/Applications/Custom_Applications/OrthoDB_source/odb10v1_gene_xrefs_extended_arth_sub.tab',sep='\t',header=None)


if len(taxid_list)>3:
    taxid_names=''
else: 
    taxid_names="_".join(taxid_list)
######################################################### Proccess

######## Key data
### clean
odb_key.columns=['odb','xref','database']
odb_key[['taxid','gene']]=odb_key['odb'].str.split(':',expand=True)
### get databases for each species


############## Node data
### split odb column into taxid and ortho group
odb_node[['taxid','gene']]=odb_node[1].str.split(':',expand=True)
odb_node.columns=['OG','junk','taxid','gene']
### subset for only reelvant taxids



################################### Choose species for tree
odb_key_sp=odb_key[odb_key['taxid'].isin(taxid_list)]

###############
dbs=[database_find(x,odb_key_sp) for x in taxid_list]
dbs_dict=dict(zip(taxid_list,dbs))

print(dbs_dict)

b=[]
for k, v in dbs_dict.items():
   a=odb_key_sp[((odb_key['database']==v) & (odb_key['taxid']==k))]
   b.append(a)
odb_key_sub=pd.concat(b)
odb_key_sub=odb_key_sub.drop(['database','gene','taxid'],axis=1)
del odb_key_sp

tax_subset=odb_node[odb_node['taxid'].isin(taxid_list)]
### count number of genes in each orthogroup
ortho_counts=tax_subset.groupby('OG')['taxid'].agg({'total_genes':'count','unique_taxids':unicount}) ## aggregate function



#del odb_node 
#del odb_key
gc.collect()


######################### ONE TO ONE ORTHOLOGUES

### Filter for one to one orthologues
all_species=ortho_counts[ortho_counts['total_genes']==len(taxid_list)] ## filter for all genes that are present in all species listed in taxid file
one_to_one=all_species[all_species['unique_taxids']==len(taxid_list)] ## filter for OGs which have exactly same number of total genes as species they are present in
og_groups=[str(x) for x in list(one_to_one.index)]
og_groups=list(one_to_one.index)
og_genes=tax_subset[tax_subset['OG'].isin(og_groups)]
og_genes.rename(columns={'junk':'odb'},inplace=True)
#del tax_subset

if arg3=='id':
    og_final=og_genes.drop(['gene'],axis=1)
    m=pd.merge(og_final,odb_key_sub,on='odb',how='left')
    m['xref']=m.apply(lambda row: row['xref'] if type(row['xref'])==str else row['odb'],axis=1)
    m=m.drop(['odb'],axis=1)
    m=m.set_index('OG')
    
    l=[]
    for ind in list(set(m.index)):
        
        sub=m[m.index==ind]
        sub=sub.drop_duplicates(subset='taxid', keep="first") ### remove duplicates. Choose random of 2 annotated IDS
        
        if sub.shape[0]!=len(taxid_list):
            continue 
        
        d=dict(zip(list(sub['taxid']),list(sub['xref'])))
        
        #sp1=str(sub[sub['taxid']==taxid_list[0]]['xref'].item())
        #sp2=str(sub[sub['taxid']==taxid_list[1]]['xref'].item())
        #d={taxid_list[0]:sp1, taxid_list[1]: sp2}
        df=pd.DataFrame(d,index=[str(ind)])
        l.append(df)
        
    final=pd.concat(l)
    final['OG']=final.index
    
    final.to_csv('temp.tsv',index=False,sep='\t')
    os.system('cat temp.tsv')

print(final.shape[0])
print(final.shape[0])


if arg3=='seq':
    recs=SeqIO.to_dict(SeqIO.parse(home+'/Applications/Custom_Applications/OrthoDB_source/odb_arthropoda_augment.faa', 'fasta'),key_function=lambda rec: rec.description)
    recs2={k:v for (k,v) in recs.items() if re.sub(':.+$','',k) in taxid_list}
    og_list=list(og_genes['odb'])
    fa_sub={k:v for (k,v) in recs2.items() if k in og_list}
    del recs 
    del recs2
    gc.collect()
    
    #fa_sub={k:v for (k,v) in SeqIO.to_dict(SeqIO.parse('/home/shanedenecke/Applications/custom/OrthoDB/odb_all_metazoa.faa', 'fasta'),key_function=lambda rec: rec.description) if k in og_list}

### write ortho group fastas
    try:
        shutil.rmtree('./og_sequences'+taxid_names)
    except:
        pass
    
    os.mkdir('./og_sequences'+taxid_names)
    for i in og_groups:
        
        og_sub=list(og_genes[og_genes['OG']==i]['odb'])
        og_fasta={k:v for (k,v) in fa_sub.items() if k in og_sub}    
        final_recs=SeqIO.parse(og_fasta.values(),format='fasta')
        
        with open('temp.faa', 'w') as handle:
            SeqIO.write(og_fasta.values(), handle, 'fasta')   
            
        newlist=list(SeqIO.parse('temp.faa', 'fasta'))
        
        handle = open('./og_sequences'+taxid_names+'/'+str(i)+'.faa', "w")
        writer = FastaWriter(handle, wrap=0)
        writer.write_file(newlist)
        handle.close()
        
        os.remove("temp.faa") 
