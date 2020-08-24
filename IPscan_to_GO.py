#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 16:09:40 2020

@author: shanedenecke
"""

import pandas as pd
pd.set_option('display.max_columns', None)
import argparse

import re

CLI=argparse.ArgumentParser()
CLI.add_argument("-ip",type=str,help='Output of interpro scan with GO terms in 13th column')
args = CLI.parse_args()


args.ip=pd.read_csv('/home/shanedenecke/Dropbox/omics_projects/August20/GO_common/HelArm_ipscan2.tsv',sep='\t',header=None)


ip=args.ip

go_samp=[list(ip.iloc[0]).index(x) for x in list(ip.iloc[0]) if re.search('GO:',str(x))][0]

ip=ip[[0,go_samp]]
ip.columns=['Geneid','GO']

ip.drop_duplicates(inplace=True)

ip=ip[ip.GO.notna()]

single=ip[ip.GO.isin([x for x in ip.GO.values if not re.search('\|',x)])]
multi=ip[ip.GO.isin([x for x in ip.GO.values if re.search('\|',x)])]

l=pd.DataFrame()
for i in range(multi.shape[0]):
    golist=multi.iloc[i]['GO'].split('|')
    genlist=[multi.iloc[i]['Geneid']]*len(golist)
    sub=pd.DataFrame(list(zip(genlist,golist)),columns=['Geneid','GO'])
    l=l.append(sub,ignore_index=True)

final=pd.concat([single,l])

print(final.to_csv(sep='\t',index=False))

