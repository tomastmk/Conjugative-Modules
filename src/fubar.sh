#!/bin/bash

alinhamentos_path=$1
trees_path=$2
output_path=$3

for file in "$alinhamentos_path"/*.fna; do
  echo "Processing $file"
  tree="$trees_path"/"$(basename "$file" .fna)".tree
  hyphy fubar \
    --code Mold-Protozoan-mtDNA \
    --alignment "$file" \
    --tree "$tree" \
    --output "$output_path"/"$(basename $file .fna).json"
  
done
