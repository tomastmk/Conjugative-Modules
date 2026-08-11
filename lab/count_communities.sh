output_path="/home/ttamaki/New/results_new/communities_count_test.tsv"

echo "metric    resolution  cut number_of_clusters" >> $output_path

for arquivo in /home/ttamaki/New/results/var_com/*.tsv
do 
    nome=$(basename "$arquivo" .tsv)
    metric=$(echo "$nome" | cut -d'_' -f1)
    r=$(echo "$nome" | grep -oP 'r\K[0-9.]+')
    c=$(echo "$nome" | grep -oP 'c\K[0-9.]+')
    n_clusters=$(tail -n 1 $arquivo | cut -f3)
    echo "$metric   $r  $c  $((${n_clusters%.*}+1))" >> $output_path
done
