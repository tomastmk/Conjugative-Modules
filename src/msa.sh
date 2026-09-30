#!/bin/bash

for i in /home/ttamaki/Conjugative-Modules/results/enrichment/clusters/proteins/*; do
  mafft $i >> /home/ttamaki/Conjugative-Modules/results/enrichment/alignments/$(basename "$i")
done
