#### tmhmm_filter_function

/home/pioannidis/Programs/tmhmm-2.0c/bin/tmhmm  $1 | grep "Number of predicted" |\
perl -pe 's/..(.+) Number of predicted TMHs:\s+(\S+)/$1\t$2/g' | awk '{ if ($2 >='$2') { print } }'
