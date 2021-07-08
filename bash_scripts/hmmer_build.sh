#!/bin/bash

for filename in ./Data/alignments/*; do
  $IFS = "."
  read -a strarr <<< "$filename"
  /Users/ashikajayanthy/Documents/bin/hmmbuild "${strarr[0]}" "$filename" || continue
done
