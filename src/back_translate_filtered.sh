#!/bin/bash

nucleotide_fasta=$1
amino_acid_fasta=$2
output=$3

for i in "$nucleotide_fasta"/*.fna; do
    basename=$(basename "$i" .fna)

    seqkit grep \
        -f <(seqkit seq -n "$i") \
        "$amino_acid_fasta/$basename.faa" \
        > "$output/$basename.faa"
done
