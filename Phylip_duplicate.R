shhh <- suppressPackageStartupMessages
shhh(library(dplyr))
shhh(library(data.table))
shhh(library(stringr))
shhh(library(ape))
shhh(library(ggtree))
shhh(library(tidyr))
shhh(library(readr))

args = commandArgs(trailingOnly=TRUE)
#args[1]='~/Transporter_ID/ABC_id/phylo/ABC_total_NBDABCA_.faa.aln.trimm.phy'

input.phylip=fread(args[1])
names=input.phylip$V1
uni.phy=unique(names)


l=list()
for(i in uni.phy){
  sub=input.phylip[names==i]
  if(length(sub$V1)>1){
    sub$V1=paste0(sub$V1,'_',1:length(sub$V1))
  }
  l[[i]]=sub
}

final=rbindlist(l)
  
  cat(format_tsv(final,col_names = F))

#fwrite(b,a,sep=' ',col.names=F)
