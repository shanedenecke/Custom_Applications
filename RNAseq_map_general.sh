#!/usr/bin/env bash


### add help 
if [ "$1" == "-h" ]; then
  echo "
  Welcome! This shell script is designed to map RNA-seq reads to an existing genome
  
  Arguments:
  
  -genome The genome that you'll be mapping your reads to
  -fastq The zipped or unzipped fastq files that you want to analyze. This should take the form of a directory with a regex. E.g. './Folder/*.gz' Need to put in double quotes for some reason
  -threads: Self explanatory
  -outdir: Output diretory where all of your outputs will be located. Note: They will be put into relevant subdiretories automatically
 "
  exit 0
fi



### add arguments
while [ "$#" -gt 0 ]; do
  case "$1" in
    -genome) GENOME="$2"; shift 2;;
    #-fastq_dir) FASTQ_DIR="$2"; shift 2;;
    #-fastq_sub) FASTQ_SUB="$2"; shift 2;;
    -fastq) FASTQ="$2"; shift 2;;
    -threads) THREADS="$2"; shift 2;;
    -outdir) OUTDIR="$2"; shift 2;;
    -gc) GC="$2"; shift 2;;
  esac
done

#cd /mnt/disk/shane/Transcriptomics/L2_plant/
#GENOME=../genome_reference/Spodoptera/mar20/SpoFru_genome_mar20.fna
#FASTQ="./raw_zipped_fastq/*.gz"
#THREADS=8
#OUTDIR=/mnt/disk/shane/Transcriptomics/Sf_midgut_cell
#./scripts/RNAseq_map_general.sh -genome /mnt/disk/shane/Transcriptomics/genome_reference/Helicoverpa/GCF_002156985.1_Harm_1.0_genomic.fna -fastq "/mnt/disk/shane/Transcriptomics/L2_plant/HelArm/*.gz" -threads 6 -outdir HelArm_L2_plant


GENOME=$(readlink -f $GENOME)
FASTQ=$(readlink -f $FASTQ)
OUTDIR=$(readlink -f $OUTDIR)

mkdir $OUTDIR
cd $OUTDIR
mkdir mapping
mkdir FastQC
mkdir mapping_rate


### Make sure genome is unzipped
if [[ $GENOME =~ .gz$ ]];then
	echo 'zipped genome'
	gunzip -k $GENOME
	GENOME2=$(echo $GENOME | sed 's/.gz//g') 
	GENOME=$GENOME2
else
	echo 'unzipped genome'
fi
 
### Make database for genome
genomedir=$(dirname $GENOME)
if ls $GENOME'_index'* 1> /dev/null 2>&1; then
	echo 'Hisat database already made'
else
	hisat2-build -q -p $THREADS $GENOME $GENOME'_index'
fi



###CREATE SAMPLE LIST
samples=$(for i in $FASTQ; do echo $i | sed -E 's/_[0-9].fas.+$//g' ;done | sort -u)

# PERFORM SAMPLE MAPPING LOOP

for sample in $samples;do 
 b=$(basename $sample)
 ### check whether a previous output exists
 
 if ls ./mapping/$b*fastq; then
  echo 'Skipping Fastq processing'
 else
  ### PERFORM TRIMMING AND UNZIP FILES
  if [[ $FASTQ =~ .gz$ ]];then
   echo 'Zipped FastQ files'
   fastp -g --detect_adapter_for_pe --thread $THREADS --in1 $sample'_1.fastq.gz' --in2 $sample'_2.fastq.gz' --out1 ./mapping/$b'_1_trimmed.fastq.gz' --out2 ./mapping/$b'_2_trimmed.fastq.gz'
   unpigz -p $THREADS -c ./mapping/$b'_1_trimmed.fastq.gz' > ./mapping/$b'_1.fastq'
   unpigz -p $THREADS -c ./mapping/$b'_2_trimmed.fastq.gz' > ./mapping/$b'_2.fastq'
  else
   fastp -g --detect_adapter_for_pe --thread $THREADS --in1 $sample'_1.fastq' --in2 $sample'_2.fastq' --out1 ./mapping/$b'_1.fastq' --out2 ./mapping/$b'_2_trimmed.fastq'
  fi
  
 ###PERFORM FASTQC
 fastqc ./mapping/$b'_1.fastq' -q -o ./FastQC
 fastqc ./mapping/$b'_2.fastq' -q -o ./FastQC 
 	
 fi ### close loop to do Fasta pre processing
 
 echo 'Finished FASTQ Processing for sample '$b
 
 ###MAPPING
 hisat2 --dta-cufflinks -p $THREADS -x $GENOME'_index' -1 ./mapping/$b'_1.fastq' -2 ./mapping/$b'_2.fastq' -S ./mapping/$b.sam > ./mapping_rate/mapping_rate_$b.txt 2>&1 
 samtools view -b --threads $THREADS ./mapping/$b.sam > ./mapping/$b.bam
 samtools sort -m 1G --threads $THREADS ./mapping/$b.bam -o ./mapping/$b.sorted.bam

 echo 'Finished Mapping for sample '$b
 
 if [ "$GC" = "T" ] || [ "$GC" = '' ]; then
  ###GARBAGE COLLECT
  rm ./mapping/$b'_1.fastq'
  rm ./mapping/$b'_2.fastq'
 fi
 rm ./mapping/$b.sam	
 rm ./mapping/$b.bam
 rm ./mapping/*.gz
done

