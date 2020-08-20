### Reciprocal blast function

## 1=starting database
## 2=target proteome
## 3=evalue cutoff
## 4=qcov cutoff
## 5=num_threads
#ARG3=${3:-1e-7}
#ARG4=${4:-60}



##initial
makeblastdb -in $2 -parse_seqids -dbtype prot 
mkdir blast
blastp -query $1 -db $2 -outfmt "6 qseqid sseqid evalue qcovs" -evalue $3 -qcov_hsp_perc $4 -num_threads $5 > ./blast/initial_blast.tsv
cat ./blast/initial_blast.tsv | cut -f 2 | sort -u > ./blast/initial_blast_names.txt
echo "Numer of genes in initial blast is: _______________" $(cat ./blast/initial_blast.tsv | cut -f 2 | sort -u | wc -l)

perl -i -pe 's/^emb.(.+).+$/$1/g' ./blast/initial_blast_names.txt
## extract initial genes
rm ./blast/initial_blast.faa
while read i
do
grep -A 1 $i $2 >> ./blast/initial_blast.faa
done < ./blast/initial_blast_names.txt

##recip blast
makeblastdb -in $1 -parse_seqids -dbtype prot 
blastp -query ./blast/initial_blast.faa -db $1 -outfmt "6 qseqid sseqid evalue qcovs" -evalue $3 -qcov_hsp_perc $4 -num_threads $5 > ./blast/recip_blast.csv
cat ./blast/recip_blast.csv | cut -f 1 | sort -u > ./blast/recip_blast_names.txt
echo "Numer of genes in Reciprocal blast" $(cat ./blast/recip_blast.csv | cut -f 1 | sort -u | wc -l)


##extract recip blast
rm ./blast/recip_blast.faa
while read i
do
grep -A 1 $i $2 >> ./blast/recip_blast.faa
done < ./blast/recip_blast_names.txt
