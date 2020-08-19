#!/usr/bin/env Rscript  
shhh <- suppressPackageStartupMessages
shhh(library(dplyr))
shhh(library(data.table))
shhh(library(edgeR))
shhh(library(argparser))
shhh(library(readr))



p=arg_parser('EdgeR_command_line')
p <- add_argument(p, "--counts", help="path to counts file")
p <- add_argument(p, "--samples",help="List of samples. Must be partial match for your column names. E.g. if your columns are Sf_MG1,Sf_MG2,Sf_C1,Sf_C2 you would put Sf_MG,Sf_C")
p <- add_argument(p, "--geneids", help="column name of gene/transcript identifiers",default='Geneid')

argv=parse_args(p)

#setwd('/home/shanedenecke/Dropbox/omics_projects/August20/Lepidoptera_L2_plant/Tables/')
#argv$counts='./Ha_clean_counts.tsv'
#argv$samples="Ha_L2_P_C Ha_L2_P_MG"
#argv$geneids='Geneid'


samps=strsplit(argv$samples,split=' |,') %>% unlist()
#print(samps[1])
#print(samps[2])



counts.raw=fread(argv$counts)

######################################### Edge R
counts.de=counts.raw %>% select(argv$geneids,contains(samps))
rownames(counts.de)=counts.raw[[argv$geneids]]
counts.de[[argv$geneids]]=NULL

group=c()
for(i in samps){
  group=c(group,rep(i,length(grep(i,colnames(counts.de)))))
}

delist=DGEList(counts=counts.de,group=group)

### filter
keep=filterByExpr(delist)
delist=delist[keep, , keep.lib.sizes=FALSE]

## Add normalization
delist=calcNormFactors(delist)
delist = estimateCommonDisp(delist)
delist = estimateTagwiseDisp(delist)

##perform test
et = exactTest(delist, pair=samps)
tTags = topTags(et,n=NULL)

output=tTags$table 
row.index=rownames(output) %>% as.numeric()
genes=counts.raw[row.index][[argv$geneids]]
output[[argv$geneids]]=genes
allgenes=select(counts.raw,argv$geneids)

all.de=data.table(sampleA=samps[1], sampleB=samps[2], output) %>% 
  arrange(PValue) %>%  
  data.table()

cat(format_tsv(all.de))




