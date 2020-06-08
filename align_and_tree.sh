### ###################### Create HMM profiles for each SLC family in species ###################
##

## argument 1= name of fasta file
## argument 2= Number of threads 
## argument 3= outgroup 
#cd /data2/shane/Documents/transporter_candidate_phylogeny/fasta
#A=ABCB_family.faa
#B=24
### Align and produce raxML tree



## use these for debugging

################################################################################

base=$(echo $(basename $1) | cut -d '.' -f 1)
mkdir $base


sed 's/-/_/g' $1 >  $base'/clean_'$base'.faa'

mafft --thread $2 $base'/clean_'$base'.faa' > $base'/'$base'.aln'
~/Applications/trimal/source/trimal -automated1 -phylip_paml  -in $base'/'$base'.aln' -out $base'/'$base'.aln.trimm.phy'
Rscript ~/Applications/Custom_Applications $base


raxfile=$base'/'$base'.aln.trimm.phy'
raxdir=$(pwd)'/'$base'/'



if [ ! -z $3 ]  ]
   then
    echo 'without outgroup'
   ~/Applications/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 500 -T $2 -m PROTGAMMAAUTO -s $raxfile -n $base'.tre' -w $raxdir

 else
    echo 'include outgroup'
    ~/Applications/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 500 -T $2 -m PROTGAMMAAUTO -s $raxfile -n $base'.tre' -w $raxdir -o $3
 fi
  
