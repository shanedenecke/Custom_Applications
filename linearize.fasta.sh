#!/usr/bin/env bash

### Create linearized Fastaa file from non linear one
### $1= input fasta; $2= output fasta

##### Linearize fasta file
sed -e 's/\(^>.*$\)/#\1#/' $1 | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' 