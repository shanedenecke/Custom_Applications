### General CRISPR script

#.libPaths('/home/shanedenecke/R/x86_64-redhat-linux-gnu-library/3.5')

#args=c('/home/shanedenecke/Dropbox/target_discovery/SLC_targets/Ha_Ph_Genes/SLC26_CR','/home/shanedenecke/Dropbox/target_discovery/SLC_targets/Ha_Ph_Genes/SLC26_early_exons.fna')
args=commandArgs(trailingOnly = T)
setwd(as.character(args[1]))




#setwd('/home/shanedenecke/Dropbox/wp3_genetic_modification/Bactrocera_gene_drive/white')
#setwd('..')
#tar.seq="./bo_w_1_2.fa"

tar.seq=as.character(args[2])
fa=readLines(tar.seq)
fa=fa[-1]
clean.seq=paste(fa,collapse="")


len=23
leftlims <- 1:(nchar(clean.seq) - (len - 1))
rightlims <- len:nchar(clean.seq)



fasta_output=function(sequences,name){
  line1=paste(">sgRNA target sequence",name)
  for(i in 1:length(sequences)){
    first=paste(line1,i)
    second=sequences[i]
    writeLines(c(first,second),con=paste(name,i,".fa",sep=""),sep="\n")
  }
}

#sequence='ACTG'
rev_comp=function(sequence){
  split=unlist(strsplit(sequence,split=""))
  r=rev(split)
  key=list(A='T',C='G',T='A',G='C')
  rc=c()
  for (i in r){
    rc=c(rc,key[[i]])
  }
  #rc=comp(r,forceToLower = F)
  final=paste(rc,collapse="")
  return(final)
}



possible=mapply(substr,clean.seq, leftlims, rightlims,USE.NAMES=FALSE)
rc.pos=sapply(possible,rev_comp,USE.NAMES = F)


all.targets=c(possible,rc.pos)

cr.targets=c()
for(i in all.targets){
  if((grepl("^G",i) & grepl("GG$",i))){
    #if(grepl("GG$",i)){
      cr.targets=c(cr.targets,i)
    }
  #}
}

#nam=unlist(strsplit(args[1],split="/"))[length(unlist(strsplit(args[1],split="/")))]
#setwd('/home/shanedenecke/Dropbox/wp3_genetic_modification/Bactrocera_gene_drive/white')
fasta_output(sequences=cr.targets,name='seq')#=args[3])
