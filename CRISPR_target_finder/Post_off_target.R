## Make thiswork with bash script
#.libPaths('/home/shanedenecke/R/x86_64-redhat-linux-gnu-library/3.4')
args=commandArgs(trailingOnly = T)
library(dplyr)
## Set wd
#args[1]='/home/shanedenecke/Dropbox/target_discovery/SLC_targets/Ha_Ph_Genes/SLC26_CR'
print(getwd())
#setwd(as.character(args[1]))
#setwd('/home/shanedenecke/Dropbox/wp3_genetic_modification/Bactrocera_gene_drive/white')
## import libraries 
#library(stringr)
#library(dplyr)
#library(data.table)

##Get list of fasta files in directory
fil=list.files()[grepl(".fa$",list.files())]

## Create empty data frame 
targs=data.frame(name=fil,sequence=rep(0,length(fil)))

## Put sequence into corresponding file name 
for(a in fil){
  targs[which(targs$name==a),'sequence']=readLines(a)[2]
}

## Add necessary all info to targs data frame. Filter based on given critera

fin=list()
for (i in 2:length(list.dirs())){
  setwd(list.dirs()[i])
  data=read.csv("sgRNA.csv",stringsAsFactors = F)[,c('X..Location','Site','Mm.Type','Mm.All')] #%>% select( X..Location,Site,Mm.Type,Mm.All)
  #data=data[complete.cases(data),]
  temp.list=strsplit(data$Mm.Type,split="",fixed=T)
  temp.list2=lapply(temp.list, function(x) x[2:3])
  split=as.data.frame(t(rbind.data.frame(temp.list2)),row.names = F)
  #as.data.frame(str_split_fixed(data$Mm.Type,"",3))[,c(2,3)]
  names(split)=c('seed','distal')
    q=cbind(data,split,stringsAsFactors=F) 
    q$seed=as.numeric(as.character(q$seed))
    q$distal=as.numeric(as.character(q$distal))
    data=subset(q,Mm.All<5)
    data=subset(data,seed<=1)
    if(nrow(data)>1){
      index=fil[i-1]
      inrow=targs[which(targs$name==index),]
      data$name=inrow$name
      data$sequence=inrow$sequence
      fin[[i]]=data
    }
  setwd("..")
}

tot=do.call("rbind",fin)
tot$name=gsub("seq","",tot$name) %>% gsub(".fa","",.) %>% as.numeric()
tot=arrange(tot,name)
nam=unlist(strsplit(args[1],split="/"))[length(unlist(strsplit(args[1],split="/")))]

#dir.create(args[1])
#setwd(args[1])
write.csv(tot,paste0(fil[1],nam,".csv"),row.names = F)


uni=cbind(unique(tot$sequence),unique(as.character(tot$name)))
colnames(uni)=c('seq','name')
write.csv(uni,paste0(fil[1],nam,"unique_seq.csv"),row.names = F)

    
