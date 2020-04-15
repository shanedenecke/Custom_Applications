#!/usr/bin/env Rscript 
shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(readr))
shhh(library(stringr))

## already got lengths from expression analysis

args = commandArgs(trailingOnly=TRUE)
#args[1]='/home/sdenecke/Transporter_ID/ABC_id/CAFE/Ultrametric_tree/Hemipteran_og_sequences'
#print(args[1])
setwd(args[1])

seq.list=list()
names.list=list()

total_records <- function(x){
  as.integer(system2("wc",args = c("-l",x," | awk '{print $1}'"),stdout = TRUE))
}

files=list.files(full.names = T)[grepl('.phy$',list.files())] 
files=files[total_records(files)>1] %>% na.omit() %>% as.character()

#### catch error
if(length(unique(sapply(files,total_records)))!=1){
  rec.counts=sapply(files,total_records)
  numspec=max(rec.counts)
  files=names(rec.counts)[rec.counts==numspec]
}

for(i in files){
 phy=fread(i,sep=' ',header=F,skip=1) %>% arrange(V1)
 seq.list[[i]]=phy$V2
 names.list[[i]]=phy$V1
}

### #### catch error ###### are all the names of the list in the correct order?
if(length(unique(names.list))!=1){
  print('Sequences not in the right order')
  stop()
}

#### Concatanate all sequences 
concat=c()
num.seq=min(300,length(seq.list))
for(i in 1:num.seq){
  concat=paste0(concat,seq.list[[i]])  
}

### Add names to sequence file
tax_names=names.list[[1]]
combined.full=cbind(tax_names,concat) %>% data.table()

##Annotate sequence file
len=nchar(concat[1])
headstring=paste0(length(tax_names),len)
colnames(combined.full)=as.character(c(dim(combined.full)[1],len))


cat(format_tsv(combined.full))
#fwrite(combined.full,'Full_species.phy',sep=' ')
