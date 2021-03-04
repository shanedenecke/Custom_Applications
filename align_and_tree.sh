### ###################### Create HMM profiles for each SLC family in species ###################
##


### add help 
if [ "$1" == "-h" ]; then
  echo "
  Welcome ! This shell script is designed to search SLC transporters in non-model arthropods
  
  Arguments:
  
  -proteins: Path to folder containing one or more Arthropod proteomes that you wish to search. Protoemes should be labled witht a 6 letter abbreviation followed by '_unigene.faa' e.g. DroMel_unigene.faa for Drosophila melanogaster 
  -threads: Self explanatory
  -outgroup: Output diretory where all of your outputs will be located. Note: They will be put into relevant subdiretories automatically
  "
  exit 0
fi

### add arguments
while [ "$#" -gt 0 ]; do
  case "$1" in
    -proteins) PROTEINS="$2"; shift 2;;
    -threads) THREADS="$2"; shift 2;;
    -outgroup) OUTGROUP="$2"; shift 2;;
  esac
done

## use these for debugging

################################################################################
base=$(echo $(basename $PROTEINS) | cut -d '.' -f 1)
mkdir $base


sed 's/-/_/g' $PROTEINS >  $base'/clean_'$base'.faa'

mafft --thread $THREADS $base'/clean_'$base'.faa' > $base'/'$base'.aln'
/home/sdenecke/Applications/trimal/source/trimal -automated1 -phylip_paml  -in $base'/'$base'.aln' -out $base'/'$base'.aln.trimm.phy'
#Rscript /home/sdenecke/Applications/Custom_Applications $base


raxfile=$base'/'$base'.aln.trimm.phy'
raxdir=$(pwd)'/'$base'/'



if [ ! -z $OUTGROUP ]  ]
   then
    echo 'without outgroup'
   #~/Applications/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 500 -T $THREADS -m PROTGAMMAAUTO -s $raxfile -n $base'.tre' -w $raxdir
 /home/sdenecke/Applications/raxml-ng --all --msa $raxfile --prefix $(realpath $raxdir/$base.nwk) --threads $THREADS  --bs-trees autoMRE{500} --model LG+G8+F --redo
 else
    echo 'include outgroup'
    /home/sdenecke/Applications/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 500 -T $THREADS -m PROTGAMMAAUTO -s $raxfile -n $base'.tre' -w $raxdir -o $OUTGROUP
 fi
  
