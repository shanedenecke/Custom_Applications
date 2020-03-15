#! /bin/bash

#1st argument new directory. name of folder that doesn't yet exist
#2nd argument is location of a fasta file to be used as target sequence
#3rd argument is genome file
mkdir -m 777 $1
Rscript ~/Applications/custom/CRISPR_target_finder/General_CRISPR.R $1 $2 ## set path to custom R script 1
cd $1

#     cd /home/shanedenecke/Dropbox/wp3_genetic_modification/Bactrocera_gene_drive/white
for i in *.fa
do	
	perl ~/Applications/CasOT-1.0/casot.pl -s 1 -t=$i -g=$3 ## set path to reference genome
done

Rscript ~/Applications/custom/CRISPR_target_finder/Post_off_target.R $1 ## Set path to custom R script 2

  