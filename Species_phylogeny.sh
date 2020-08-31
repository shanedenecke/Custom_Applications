#!/usr/bin/env bash

## add help 
if [ "$1" == "-h" ]; then
  echo "
  Welcome ! This shell script is designed to create phyogenetic trees from species. It takes as inputs a list of species (either taxid codes or species names) and outputs a tree.
  
  Arguments:
  
  -taxid: File path for taxid codes. For OrthoDB must be taxid numbers (e.g. 12345_0). For  must be 6 letter species abbreviations
  -algo: Algorithm to find 1:1 orthologues. Can be either 'OrthoDB' or 'OrthoFinder'. These have predefined paths on chrysalida to precomputed files/proteoms
  -outgroup: OUTGROUP used for species tree if any
  -threads: Self explanatory
  -proteomes: directory containing all proteomes to copy from
  -outdir: directory for output
  -maxseqs: Maximum number of sequences to use for tree"
  exit 0
fi



while [ "$#" -gt 0 ]; do
  case "$1" in
    -taxid) TAXID="$2"; shift 2;;
    -algo) ALGO="$2"; shift 2;;
    -outgroup) OUTGROUP="$2"; shift 2;;
    -proteomes) PROTEOMES="$2"; shift 2;;
    -threads) THREADS="$2"; shift 2;;
    -outdir) OUTDIR="$2"; shift 2;;
    -maxseqs) MAXSEQS="$2"; shift 2;;
  esac
done

### Debugging
#TAXID=./CAFE/taxid_lists/Hemimetabola_taxid_codes.txt
#ALGO='Orthofinder'
#OUTGROUP='None'
#THREADS=10
#PROTEOMES=./proteomes
#OUTDIR=./CAFE

### set defaults
if [ -z $ALGO ];then ALGO="Orthofinder"; fi
if [ -z $OUTGROUP ];then OUTGROUP="None"; fi
if [ -z $PROTEOMES ];then PROTEOMES=~/Applications/Custom_Applications/OrthoDB_source/species_proteomes; fi
if [ -z $THREADS ];then THREADS=4; fi
if [ -z $OUTDIR ];then OUTDIR=./; fi



#Set up directory tree
base=$(echo $(basename $TAXID | cut -f 1 -d '.'))
basedir=$OUTDIR/$base
mkdir $basedir
mkdir $basedir/rax_output
mkdir $basedir/one_to_one

### output one to one orthologues either from orthoDB or OrthoFinder
if [ $ALGO == 'OrthoDB' ]; then
  echo 'Starting OrthoDB search'
  ~/Applications/Custom_Applications/odb_parse.py -taxid $TAXID -output seq --outdir $basedir/one_to_one --maxseqs 200
elif [ $ALGO == 'Orthofinder' ]; then
    echo 'Starting Orthofinder search'
    mkdir $basedir/tempseqs
    grep -f $TAXID ~/Applications/Custom_Applications/OrthoDB_source/taxid_sp_convert.tsv | cut -f 2 | while read i;do 
     cp $PROTEOMES/$i'_unigene.faa' ./$basedir/tempseqs/$i'_unigene.faa'
     sed -i 's/\.//g' ./$basedir/tempseqs/$i'_unigene.faa'
    done ### copy proteome files 
    
    ~/Applications/OrthoFinder/orthofinder -t $THREADS -og -f $basedir/tempseqs -o $basedir/orthofinder_temp
    #cat $basedir/tempseqs/*.faa > $basedir/tempseqs/total_proteome.faa ### create master proteome 
    python3 ~/Applications/Custom_Applications/Orthofind_parse.py -outdir $basedir/one_to_one -indir $basedir/orthofinder_temp -total_fasta $basedir/tempseqs/ -maxseqs $MAXSEQS -mode seq
fi

### perform alignments for all one to ones 
for x in  $basedir/one_to_one/*.faa
do
  cat $x | sed 's/\./_/g' | mafft --quiet --thread $THREADS - > $x'.aln'
  ~/Applications/trimal/source/trimal -automated1 -phylip_paml -in $x.aln -out $x.phy
done
  
#### merge all phylip files   
Rscript ~/Applications/Custom_Applications/Phylip_merge.R $basedir/one_to_one/ > $basedir/one_to_one/Full_species.phy
sed -i 's/J/A/g' $basedir/one_to_one'/Full_species.phy'
sed -i 's/\./A/g' $basedir/one_to_one'/Full_species.phy'

### Perform RaxML tree
  if echo $OUTGROUP | grep -Eq "None|none" ; then
  echo "Running without OUTGROUP"
  #~/Applications/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 100 -T $THREADS -m PROTGAMMAAUTO -s $basedir/one_to_one'/Full_species.phy' -n $base.nwk -w $(realpath $basedir/rax_output/)
  ~/Applications/raxml-ng --all --msa $basedir/one_to_one'/Full_species.phy' --prefix $(realpath $basedir/rax_output/$base.nwk) --threads $THREADS  --bs-trees autoMRE{200} --model LG+G8+F --redo
else
  echo "Running with OUTGROUP"
    ~/Applications/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 100 -T $THREADS -m PROTGAMMAAUTO -s $basedir/one_to_one'/Full_species.phy' -n $base.nwk -w $(realpath $basedir/rax_output/) -o $OUTGROUP 
fi
