#!/usr/bin/env bash


### add help 
if [ "$1" == "-h" ]; then
  echo "
  Welcome! This shell script is designed to map RNA-seq reads to an existing genome
  
  Arguments:
  
  -bam The BAM files that you will use for featurecounts. E.g. './Folder/*.bam' Need to put in double quotes for some reason
  -gff The genome annotation in GFF or GTF format
  -threads: Self explanatory
  -outdir: Output diretory where all of your outputs will be located. Note: They will be put into relevant subdiretories automatically
 "
  exit 0
fi



### add arguments
while [ "$#" -gt 0 ]; do
  case "$1" in
    -bam) BAM="$2"; shift 2;;
    -gff) GFF="$2"; shift 2;;
    -threads) THREADS="$2"; shift 2;;
    -outdir) OUTDIR="$2"; shift 2;;
    -gene) GENE="$2"; shift 2;;
  esac
done

#BAM=./HelArm/mapping/*.bam
#GFF=../genome_reference/Helicoverpa/harmi.gtf
#THREADS=4
#GENE
#OUTDIR=./HelArm/FeatureCounts

mkdir $OUTDIR
featureCounts -F 'GTF' -M -s 2 -T $THREADS -p -t CDS -g $GENE -a $GFF -o $OUTDIR/Counts_output.tsv $BAM > $OUTDIR/Counts_output_summary 2>&1
counts_to_tpm.R $OUTDIR/Counts_output.tsv > $OUTDIR/TPM_output.tsv

