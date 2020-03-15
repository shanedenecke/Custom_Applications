##working
touch temp.fa
IFS=$'\n'; 
for next in $(cat $1)
do 
  esearch -db protein -query ${next} | efetch -db protein -format fasta >> $2
done
