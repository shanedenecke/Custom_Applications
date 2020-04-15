#!/usr/bin/env bash
### argument 1 list of taxid codes. for OrthoDB must be file of taxid_0 format. For Orthofinder must be file with 6 letter abbreviations 
### argument 2 list of source of orthology. Must be either OrthoDB or Orthofinder. 
### argument 3 outpgroups
#Species_phylogeny.sh -taxid_codes Hemispec.txt -ortho_algo Orthofinder -outgroups 7227_0 -threads 10

## add help 
if [ "$1" == "-h" ]; then
  echo "
  Welcome ! This shell script is designed to create phyogenetic trees from species. It takes as inputs a list of species (either taxid codes or species names) and outputs a tree.
  
  Arguments:
  
  -taxid_codes: File path for taxid codes. For OrthoDB must be taxid numbers (e.g. 12345_0). For Orthofinder must be 6 letter species abbreviations
  -ortho_algo: Algorithm to find 1:1 orthologues. Can be either 'OrthoDB' or 'OrthoFinder'. These have predefined paths on chrysalida to precomputed files/proteoms
  -outgroups: outgroups used for species tree if any
  -threads: Self explanatory"
  exit 0
fi



while [ "$#" -gt 0 ]; do
  case "$1" in
    -taxid_codes) taxids="$2"; shift 2;;
    -ortho_algo) algo="$2"; shift 2;;
    -outgroups) outgroups="$2"; shift 2;;
    #-output_directory) output_directory="$2"; shift 2;;
    -threads) THREADS="$2"; shift 2;;
  esac
done


#taxids='Hemispec.txt'
#ortho_algo='Orthofinder'
#outgrpuos='None'
#THREADS=10

#Set up directory tree
base=$(echo $(basename $taxids | cut -f 1 -d '.'))
mkdir $base
mkdir $base/rax_output
mkdir $base/one_to_one

### output one to one orthologues either from orthoDB or OrthoFinder
if [ $algo == 'OrthoDB' ]; then
  echo 'Starting OrthoDB search'
  ~/Applications/Custom_Applications/odb_parse.py -taxid $taxids -output seq --outdir $base/one_to_one --maxseqs 200
elif [ $algo == 'Orthofinder' ]; then
  echo 'Starting Orthofinder search'
  mkdir $base/tempseqs
  grep -f $taxids ~/Applications/Custom_Applications/OrthoDB_source/taxid_sp_convert.tsv | cut -f 2 | while read i;do cp ~/Applications/Custom_Applications/OrthoDB_source/species_proteomes/$i'_unigene.faa' ./$base/tempseqs/$i'_unigene.faa'; done ### copy proteome files 
  ~/Applications/OrthoFinder/orthofinder -og -f $base/tempseqs -o $base/orthofinder_temp
  #cat $base/tempseqs/*.faa > $base/tempseqs/total_proteome.faa ### create master proteome 
  python3 ~/Applications/Custom_Applications/Orthofind_parse.py -outdir $base/one_to_one -inputdir $base/orthofinder_temp -total_fasta $base/tempseqs/ -maxseqs 200
fi

### perform alignments for all one to ones 
for x in  $base/one_to_one/*.faa
do
  cat $x | sed 's/\./_/g' | mafft --quiet --thread $THREADS - > $x'.aln'
  ~/Applications/trimal/source/trimal -automated1 -phylip_paml -in $x.aln -out $x.phy
done
  
#### merge all phylip files   
Rscript ~/Applications/Custom_Applications/Phylip_merge.R $base/one_to_one/ > $base/one_to_one/Full_species.phy
sed -i 's/J/A/g' $base/one_to_one'/Full_species.phy'
sed -i 's/\./A/g' $base/one_to_one'/Full_species.phy'

### Perform RaxML tree
if echo $outgroups | grep -q "None|none" ; then
  echo "Running without outgroups"
  ~/Applications/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 100 -T $THREADS -m PROTGAMMAAUTO -s $base/one_to_one'/Full_species.phy' -n $base.nwk -w $(realpath $base/rax_output/)
else
  echo "Running with outgroups"
  ~/Applications/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 100 -T $THREADS -m PROTGAMMAAUTO -s $base/one_to_one'/Full_species.phy' -n $base.nwk -w $(realpath $base/rax_output/) -o $outgroups 
fi
