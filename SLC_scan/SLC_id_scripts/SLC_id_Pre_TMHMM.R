#!/usr/bin/Rscript


######################### SETUP

### Import packages
shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(dplyr))
shhh(library(seqinr))
shhh(library(argparser))

### functions
dash.remove=function(x){
  split=unlist(strsplit(x,'_'))
  toge=paste(split[1],split[2],sep='_')
  final=toge
  return(final)
}


### Set Source directory
getScriptPath <- function(){
    cmd.args <- commandArgs()
    m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
    script.dir <- dirname(regmatches(cmd.args, m))
    if(length(script.dir) == 0) stop("can't determine script dir: please call the script with Rscript")
    if(length(script.dir) > 1) stop("can't determine script dir: more than one '--file' argument detected")
    return(script.dir)
}
scriptPath=getScriptPath()
sourcePath=dirname(scriptPath)

### Get arguments. Only argument is the metadata file with columns "Species_name" and  "abbreviation"
p=arg_parser('SLC_pre_TMHMM')
p <- add_argument(p, "--prot", help="Input proteome")
p <- add_argument(p, "--ab", help="abbreviation")
argv=parse_args(p)

proteome=argv$prot
abbreviation=argv$ab


### Debugging
#setwd('/mnt/disk/shane/temp/SLC_test/bommor')
#scriptPath='/home/sdenecke/Applications/Custom_Applications/SLC_scan/SLC_id_scripts'
#sourcePath=dirname(scriptPath)
#proteome='/mnt/disk/shane/temp/SLC_test/bommor_unigene.faa'
#abbreviation="bommor"

print(sourcePath)
print(scriptPath)

#### Create new directories for outputs
dir.create('./TMHMM_filter',showWarnings=FALSE)

## read in metadata
meta.data=fread(paste0(sourcePath,'/SLC_id_reference/Arthropod_species_metadata.tsv'),header=T,sep='\t')


########## GET PRELIMINARY LIST AND SUBSET FASTA

dict=fread('./Preliminary_SLC_table.csv')
fam=sapply(dict$name,dash.remove)
dict$abbreviation=abbreviation
dict$family=fam
#dict2=merge(dict,select(meta.data,Species_name,abbreviation),by='abbreviation',use.names=T) %>% mutate(name)
dict2=dict
dict2$name=paste(dict2$abbreviation,dict2$code,dict2$name,sep='___')
dict2$name=gsub(' ','_',dict2$name)
dict2$code=as.character(dict2$code)
ind.rename.dict= dict2 %>% select(name,code)
writeLines(ind.rename.dict$code,'./TMHMM_filter/slc_unfiltered_codes.txt')
fwrite(ind.rename.dict,'./TMHMM_filter/temp_rename_dict.csv')
file.copy(proteome,'./TMHMM_filter/temp_proteome.faa',overwrite = T)

system(paste0(scriptPath,'/unigene_fa_sub.sh ./TMHMM_filter/temp_proteome.faa ./TMHMM_filter/slc_unfiltered_codes.txt > ./TMHMM_filter/SLC_unfiltered_all_raw.faa'))
system(paste0(scriptPath,'/fasta_rename.py ./TMHMM_filter/SLC_unfiltered_all_raw.faa ./TMHMM_filter/temp_rename_dict.csv >> ./TMHMM_filter/Renamed_unfiltered_SLC.faa 2> TM_errors.txt'))

#unfilterd_dict=rbindlist(l,use.names = T)
fwrite(dict2,'./TMHMM_filter/Renamed_unfiltered_SLC_dict.csv')


