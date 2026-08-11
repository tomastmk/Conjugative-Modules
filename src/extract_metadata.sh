while read accession; do
    pixi run efetch \
        -db nuccore \
        -id "$accession" \
        -format gb
done < /home/ttamaki/New/data/plasmid_accessions.txt > /home/ttamaki/New/results/plasmids_metadata.gb

