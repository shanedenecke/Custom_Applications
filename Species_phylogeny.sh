#!/usr/bin/env bash
### argument 1 list of taxid codes. for OrthoDB must be file of taxid_0 format. For Orthofinder must be file with 6 letter abbreviations 
### argument 2 list of source of orthology. Must be either OrthoDB or Orthofinder. 
### argument 3 outpgroups
while [ "$#" -gt 0 ]; do
  case "$1" in
    -taxid_codes) taxids="$2"; shift 2;;
    -ortho_algo) algo="$2"; shift 2;;
    -outgroups) outgroups="$2"; shift 2;;
    #-output_directory) output_directory="$2"; shift 2;;
    -threads) THREADS="$2"; shift 2;;
  esac
done

## add help 
if [ "$1" == "-h" ]; then
  echo "
  -taxid_codes: File path for taxid codes. For OrthoDB must be taxid numbers (e.g. 12345_0). For Orthofinder must be 6 letter species abbreviations
  -ortho_algo: Algorithm to find 1:1 orthologues. Can be either 'OrthoDB' or 'OrthoFinder'. These have predefined paths on chrysalida to precomputed files/proteoms
  -outgroups: outgroups used for species tree if any
  -threads: Self explanatory"
  exit 0
fi

#Set up directory tree
base=$(echo $(basename $taxids))
mkdir $base
mkdir $base/rax_output
mkdir $base/one_to_one

### output one to one orthologues either from orthoDB or OrthoFinder
if[ $algo = 'OrthoDB' ]; then
  ~/Applications/Custom_Applications/one_to_one_ID_exec.py -taxid $taxids -output seq -outdir $base/one_to_one
elif[ $algo = 'Orthofinder' ] 
  mkdir $base/tempseqs
  mkdir $base/orthofinder_temp
  cat $taxids | while read i; do cp ~/Applications/Custom_Applications/OrthoDB_source/species_proteomes/$i'_unigene.faa' ./tempseqs/$i'_unigene.faa'; done 
  ~/Applications/Orthofinder/orthofinder -f tempseqs -o orthofinder_temp
  cat tempseqs/*.faa > tempseqs/total_proteome.faa
  python3 ~/Applications/Custom_Applications/Orthofind_parse.py -outdir $base/one_to_one -inputdir orthofinder_temp -total_fasta tempseqs/total_proteome.faa
fi

### perform alignments for all one to ones 
for x in  $output_directory/*.faa
do
  cat $x | sed 's/\./_/g' | mafft --quiet --thread $THREADS - > $x'.aln'
  ~/Applications/trimal/source/trimal -automated1 -phylip_paml -in $x.aln -out $x.phy
done
  
#### merge all phylip files
Rscript ~/Applications/Custom_Applications/Phylip_merge.R $base/one_to_one > $base/one_to_one'/Full_species.phy'
sed -i 's/J/A/g' $base/one_to_one'/Full_species.phy'
sed -i 's/\./A/g' $base/one_to_one'/Full_species.phy'

### Perform RaxML tree
~/Applications/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 100 -T $THREADS -m PROTGAMMAAUTO -s $x.phy -n $b.nwk -w mkdir $base/rax_output -o $outgroups

