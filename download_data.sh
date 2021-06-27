#!/bin/bash

mkdir ./Data
cd Data
pwd=$(pwd)

declare -a taxonomic_group=( "viral" "archaea" "bacteria" "fungi" "invertebrate/Caenorhabditis_elegans" "invertebrate/Drosophila_melanogaster/" "vertebrate_other/Xenopus_tropicalis/" "vertebrate_mammalian/Mus_musculus" "vertebrate_mammalian/Homo_sapiens" )

for i in "${taxonomic_group[@]}"
do

curl 'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/${taxonomic_group}/assembly_summary.txt' | \
awk '{FS="\t"} !/^#/ {print $20} ' | \
sed 's|(ftp://ftp.ncbi.nlm.nih.gov/genomes/all/.+/)(GCF_.+)|\1\2/\2_genomic.fna.gz|' > genomic_file

mkdir ${taxonomic_group[@]}
cd ${taxonomic_group[@]}
wget --input genomic_file
cd pwd

done
