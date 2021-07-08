#!/bin/bash

for i in {1..1000};do
for filename in ./Data/alignments/*; do
  /Users/ashikajayanthy/Documents/bin/hmmemit "$filename" >> ./Data/hmm_sequences/$i.fasta
done
done
