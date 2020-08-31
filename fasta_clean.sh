#!/usr/bin/env bash



### add help 
if [ "$1" == "-h" ]; then
  echo "
  This command removes all crap (e.g. "-" "J" ".") from fasta files
  
  Arguments:
  
  -proteome: A proteome that you want to clean
  "
  exit 0
fi

### add arguments
while [ "$#" -gt 0 ]; do
  case "$1" in
    -proteome) PROTEOME="$2"; shift 2;;
  esac
done


codes="[^A|^C|^D|^E|^F|^G|^H|^I|^K|^L|^M|^N|^P|^Q|^R|^S|^T|^V|^W|^Y|^Z]"



sed -i -E "/^[^>]/s/$codes//g" $PROTEOME 
