#!/usr/bin/Rscript

shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(readr))

## collate all SLC transporter lists from human and drosophila searches 

#sp.list=c()
#for(file in list.files('./proteomes')){
#  if(grepl(".fa",file)){
#    sp.list=c(sp.list,unlist(strsplit(file,'_'))[1])
#  }
#}
#sp.list=c('Mp','Tc','Ld','Lm','Hh','Bm','Ha','Bg','Ac','Am','Ag')


human.search=fread('./HUMAN_search/final_output/SLC_final_output.csv')
human.search$name=gsub("(SLC_.+_)[0-9]+$",'\\1',human.search$name)
dros.search=fread('./DROSOPHILA_search/final_output/SLC_final_output.csv')
dros.search$name=gsub("(SLC_.+_)[0-9]+$",'\\1',dros.search$name)




unique_hits=rbindlist(list(human.search,dros.search)) %>% unique()

l3=list()
for(k in unique(unique_hits$name)){
	sub=subset(unique_hits,name==k)
	mem=seq(1,nrow(sub),by=1)
	sub$name=paste(sub$nam,mem,sep='')
	l3[[k]]=sub
}
named.hits=rbindlist(l3)


named.hits=named.hits[!duplicated(named.hits$code),]
fwrite(named.hits,'./Preliminary_SLC_table.csv')



