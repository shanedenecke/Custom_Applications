#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 16:12:55 2021

@author: shanedenecke
"""

###################### Import Sequences ######################
import os
import argparse
import pandas as pd
import re
import subprocess as sp
import shutil


###################### Get Source Directory ######################
source_dir=os.path.dirname(os.path.abspath(__file__))
#source_dir='/home/shanedenecke/Applications/Custom_Applications/Wiggle_Index_local'


###################### Parse Arguments ############################################
CLI=argparse.ArgumentParser()
CLI.add_argument("-wiggle_folder",type=str,default='user_input',help='Path to the folder containing Wiggle Videos')
CLI.add_argument("-key",type=str,default='user_input',help='FilePath to CSV conatining Wiggle renaming file')

args = CLI.parse_args()

#args.wiggle_folder='/home/shanedenecke/Dropbox/Wiggle_Index/10Feb21/'
#args.key='/home/shanedenecke/Dropbox/Wiggle_Index/10Feb21_Wiggle_rename.csv'


###################### Import sequence folders ############################################
basefolder=args.wiggle_folder
os.chdir(basefolder)

###################### Sort wigglevids  ######################
wiggle_vids=os.listdir()
wiggle_vids.sort()

###################### Import Naming key ############################################
name_key=pd.read_csv(args.key,dtype=str)


###################### Clean names ######################
full = list(name_key.apply('__'.join, axis=1))
try:
	shutil.rmtree('renamed_videos')
except:
	print('Creating reanmed videos directory')

os.makedirs('renamed_videos')

for i in range(0,len(full)):
    oldname=wiggle_vids[i]
    newname=full[i]+'.mp4'
    
    if len([x for x in os.listdir(basefolder) if re.search(os.path.basename(oldname),x)])>0:
        shutil.copyfile(oldname,'./renamed_videos/'+newname)
        #os.rename(oldname,newname)
    

###################### Convert renamed videos to image stacks ######################
try:
    shutil.rmtree(source_dir+'/Wiggle_Analyze/')
except:
    print('Wiggle_Analyze Directory doesnt yet exist')
    os.makedirs(source_dir+'/Wiggle_Analyze/')
    

for i in os.listdir(basefolder+'/renamed_videos'):
    newdir=source_dir+'/Wiggle_Analyze/'+i.replace('.mp4','')
    try:
         os.mkdir(newdir) 
         sp.run(['ffmpeg','-i','./renamed_videos/'+i,'-vf','fps=25',newdir+'/'+i.replace('.mp4','')+'_%d.jpg']) 
    except: 
        print('directory exists')


