#!/bin/bash

output_dir="$1"
alinhamentos_dir="$2"
mkdir -p $output_dir

for f in "$alinhamentos_dir"/*.faa; do
    echo "$f"
    FastTree -lg $f > "$output_dir"/"$(basename $f .faa).tree"
done
