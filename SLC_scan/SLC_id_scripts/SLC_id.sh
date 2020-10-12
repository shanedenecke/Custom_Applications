#!/usr/bin/env bash


### add help 
if [ "$1" == "-h" ]; then
  echo "
  Welcome ! This shell script is designed to search SLC transporters in non-model arthropods
  
  Arguments:
  
  -proteome: fasta file containing amino acid sequences to search 
  -threads: Self explanatory
  -outdir: Output diretory where all of your outputs will be located. Note: They will be put into relevant subdiretories automatically
  example
  -abb: Abbreviation for species
  ./SLC_ID_SCRIPTS/SLC_id_standalone/SLC_id_scripts/SLC_id.sh -proteomes $H/proteomes -busco_thresh 75 -threads $THREADS -outdir $H -metadata /mnt/disk/shane/Transporter_ID/SLC_id_pipeline/SLC_ID_SCRIPTS/GENERAL_REFERENCE/keys/Arthropod_species_metadata.tsv
  "
  exit 0
fi

### add arguments
while [ "$#" -gt 0 ]; do
  case "$1" in
    -proteome) PROTEOME="$2"; shift 2;;
    -threads) THREADS="$2"; shift 2;;
    -outdir) OUTDIR="$2"; shift 2;;
    -abb) ABB="$2"; shift 2;;
  esac
done

### Establish directory for scripts and reference files common to this shell script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"


### For debugging
#PROTEOME=./proteomes/HelZea_unigene.faa
#THREADS=14
#OUTDIR=/mnt/disk/shane/Transporter_ID/SLC_test/HelZea_test
#SCRIPT_DIR=/home/sdenecke/Applications/Custom_Applications/SLC_scan/SLC_id_scripts
#ABB='HelZea'


SOURCE_DIR="$(dirname "$SCRIPT_DIR")"


PROTEOME=$(readlink -f $PROTEOME)
OUTDIR=$(readlink -f $OUTDIR)
echo 'The proteome file is '$PROTEOME
echo 'The Output directory is '$OUTDIR
echo 'The source directory is '$SOURCE_DIR

#### Start analysis
mkdir -p $OUTDIR
cd $OUTDIR
mkdir -p $OUTDIR/SLC_search_intermediates
mkdir $OUTDIR/Human_search
mkdir $OUTDIR/Drosophila_search

target_species=$(echo $(basename $PROTEOME) | cut -d '_' -f 1) 
echo 'Target Species is '$target_species
$SCRIPT_DIR/SLC_HMM_Search.sh -database $SOURCE_DIR/SLC_id_reference/HomSap_Database -target $PROTEOME -out $OUTDIR/'HUMAN_search' -threads $THREADS
$SCRIPT_DIR/SLC_HMM_Search.sh -database $SOURCE_DIR/SLC_id_reference/DroMel_Database -target $PROTEOME -out $OUTDIR/'DROSOPHILA_search' -threads $THREADS
$SCRIPT_DIR/SLC_crossref_human_dros_searches.R ### Crossreference Human and Drosophila searches and output to "preliminary SLC dicts" folder


####################### 6) Post Process SLC tables
$SCRIPT_DIR/SLC_id_Pre_TMHMM.R --prot $PROTEOME --ab $ABB
tmhmm  ./TMHMM_filter/Renamed_unfiltered_SLC.faa > $OUTDIR/TMHMM_filter/SLC_TMHMM_full.txt
cat ./TMHMM_filter/SLC_TMHMM_full.txt | grep 'Number of predicted' | perl -pe 's/^..([A-z].+) Number of predicted TMHs:\s+(\S)/$1\t$2/g' > $OUTDIR/TMHMM_filter/SLC_TMHMM_scores.txt
$SCRIPT_DIR/SLC_id_post_TMHMM.R 




##### Clean everything up
mv $OUTDIR/TMHMM* $OUTDIR/SLC_search_intermediates
mv $OUTDIR/*_search $OUTDIR/SLC_search_intermediates
mv $OUTDIR/Final_outputs $OUTDIR/SLC_search_intermediates
rm TM_errors.txt
mv Preliminary_SLC_table.csv $OUTDIR/SLC_search_intermediates

cp ./SLC_search_intermediates/Final_outputs/Total_dictionary.csv ./$ABB'_SLC_dictionary.csv'
cp ./SLC_search_intermediates/Final_outputs/Total_raw_fasta.faa ./$ABB'_SLC_fasta.faa'


cd -


