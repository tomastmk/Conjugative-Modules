#!/bin/bash

output="/home/ttamaki/Conjugative-Modules/results/enrichment/alignments/back_tracked_filtered"

mkdir -p $output

for f in /home/ttamaki/Conjugative-Modules/results/enrichment/alignments/back_tracked/*.fna; do 
  codon_filtered="$output/$(basename "$f")"
  maxalign-rs "$f" "$codon_filtered"
  trimal -maxidentity 0.99 -in "$codon_filtered" -out "$codon_filtered"

done



