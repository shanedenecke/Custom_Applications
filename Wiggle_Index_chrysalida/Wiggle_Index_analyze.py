#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 16:13:36 2021

@author: shanedenecke
"""

###################### Import Sequences ######################
import os
import argparse
import pandas as pd
import re
import shutil
import math
import matplotlib.pyplot as plt
import statistics as st


source_dir=os.path.dirname(os.path.abspath(__file__))

################### Parse Arguments ############################################
CLI=argparse.ArgumentParser()
CLI.add_argument("-wiggle_folder",type=str,default='user_input',help='Path to the folder containing Wiggle Videos')

args = CLI.parse_args()


###################### Import sequence folders ############################################
basefolder=args.wiggle_folder
os.chdir(basefolder)

###################### Move wiggle outputs into specific folder ######################
try:
    os.makedirs('Wiggle_output')
except:
    print('Wiggle ouput folder already made')
os.system('cp '+source_dir+'/Wiggle_Analyze/zCropped*/Output/Wiggle* '+basefolder+'/Wiggle_output/')



###################### Combine Wiggle index files ######################
l=[pd.read_csv('./Wiggle_output/'+x,sep='\t',header=0) for x in os.listdir('./Wiggle_output/')]
rawData=pd.concat(l)
rawData['Image Name']=[re.sub('(rep[0-9]+)','\\1__',x) for x in rawData['Image Name']]
rawData2=rawData[['Image Name','Wiggle Index']]
rawData2[['Date','Time','Genotype','Pesticide','Dose','Rep','Location']]=rawData2['Image Name'].str.split('__',expand=True)
rawData3=rawData2.drop('Image Name',axis=1)
rawData3['Location']=[x.replace('.tif','') for x in rawData3['Location']]

rawData3.to_csv('./Wiggle_output/Combined_Raw.csv',index=None)


rawdata4=rawData3.pivot_table(values='Wiggle Index',columns='Time',index=['Date','Genotype','Pesticide','Dose','Rep','Location']).reset_index(level=[0,1,2,3,4,5], drop=False)
#rawdata4.dropna()
rawdata4.to_csv('./Wiggle_output/rawdata4.csv',index=None)


times=[x for x in rawdata4.columns if re.search('min',x)]
for i in times:
    rawdata4[i+'_RMR']=rawdata4[i]/rawdata4['0min']

rmrs=[x for x in rawdata4.columns if re.search('RMR',x)]

rmr_table=rawdata4.melt(id_vars=['Date','Genotype','Pesticide','Dose','Rep','Location'],value_vars=rmrs,var_name='Time',value_name='RMR')
rmr_table['Time']=[x.replace('_RMR','') for x in rmr_table['Time']]

raw_table=rawdata4.melt(id_vars=['Date','Genotype','Pesticide','Dose','Rep','Location'],value_vars=times,var_name='Time',value_name='Raw')

final_wiggle=pd.merge(rmr_table,raw_table)


def stderr(x):
    return st.stdev(x)/math.sqrt(len(x))
    
aggregations={
    'RMR_mean': lambda x: st.mean(x),
    'RMR_stderr':lambda x: stderr(x)
    }

sum_wiggle=final_wiggle.groupby(['Date','Genotype','Pesticide','Time','Dose'],as_index=False)['RMR'].agg(aggregations)
sum_wiggle.to_csv('./Wiggle_output/Wiggle_Summary.csv',index=None)

sum_wiggle['Time']=[int(x.replace('min','')) for x in sum_wiggle["Time"]]
sum_wiggle=sum_wiggle.sort_values('Time').reset_index(drop=True)

################## Make plots ##################

##############doses=#################
try:
    shutil.rmtree('Wiggle_Plots')
    os.makedirs('Wiggle_Plots')
except:
    os.makedirs('Wiggle_Plots')

for j in set(sum_wiggle.Pesticide):
    plot_wiggle=sum_wiggle[sum_wiggle.Pesticide==j]
    for i in set(plot_wiggle.Dose):
        plot_wiggle=plot_wiggle[plot_wiggle.Dose==i]
        plt.figure(figsize=(24,12))
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.xlabel('\nTime (minutes)',fontsize=30)
        plt.ylabel('Relative Movement Ratio\n',fontsize=30)
        plt.title('Wiggle Index Pilot',fontsize=40)
        targets = set(list(plot_wiggle.Genotype))
        colors=['r','b','g']
         
         
         
         
        for target, color in zip(targets,colors):
            indicesToKeep = plot_wiggle.Genotype == target
            plt.scatter(plot_wiggle.loc[indicesToKeep, 'Time']
                        , plot_wiggle.loc[indicesToKeep, 'RMR_mean'], c = color, s = 150,label=target)
            plt.plot(plot_wiggle.loc[indicesToKeep, 'Time']
                        , plot_wiggle.loc[indicesToKeep, 'RMR_mean'], c = color)
            plt.errorbar(plot_wiggle.loc[indicesToKeep, 'Time'], 
                         plot_wiggle.loc[indicesToKeep, 'RMR_mean'], 
                         yerr = plot_wiggle.loc[indicesToKeep, 'RMR_stderr'], 
                         fmt ='.',c = color,capsize=10) 
        plt.legend(prop={'size': 25})
        plt.savefig('./Wiggle_Plots/'+j+'Wiggle_Plot.pdf')
        plt.show()

