#!/usr/bin/env Rscript 
shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(dplyr))
shhh(library(readr))

args=commandArgs(trailingOnly = T)
#args[1]='./counts_output/counts.txt'


counts.data=fread(as.character(args[1])) #%>% select(Geneid,Length,matches('out.bam'))

cols=colnames(counts.data)
iter=grep("out.bam",cols)
normlen=counts.data$Length/1000

tpm.data=select(counts.data,Geneid,Length)

for(i in iter){
  relevant.column=counts.data[[i]]
  rpk=relevant.column/normlen
  sc=sum(rpk)/1000000
  tpm=rpk/sc
  colnam1=unlist(strsplit(cols[i],'/'))
  colnam2=colnam1[length(colnam1)-1]
  tpm.data[[colnam2]]=tpm
}

cat(format_tsv(tpm.data))