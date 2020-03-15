#!/usr/bin/env bash
sed 's/\.//g' $1 | sed 's/\*//g' | sed 's/ /_/g' | sed 's/\\//g' | tr '[:lower:]' '[:upper:]' > temp.fasta
awk '/^>/{id=$0;getline;arr[id]=$0}END{for(id in arr)printf("%s\n%s\n",id,arr[id])}' temp.fasta 
rm temp.fasta
