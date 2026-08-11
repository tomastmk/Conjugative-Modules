#!/bin/bash

# Esse nome vai aparecer na execução do squeue
#SBATCH --job-name=clust

# Esse é o número de processos
#SBATCH --ntasks=1

# Essa é a memória que seu processo pode usar
#SBATCH --mem=80GB

# Esse é o número de cpus
#SBATCH --cpus-per-task=8

# Essa é a estimativa de tempo. NOTE que seu processo será encerrado após esse tempo limite
#SBATCH --time=70:00:00

# Esses são os arquivos de saída do StdOut e do StdErr, respectivamente
#SBATCH --output=/home/ttamaki/New/out/%j_clustering.out
#SBATCH --error=/home/ttamaki/New/err/%j_clustering.err

# Para exportar seu ambiente pixi
export PATH="/home/ttamaki/Conjugative-Modules/.pixi/bin:$PATH"

echo -e "\n## Job iniciado em $(date +'%d-%m-%Y as %T') ##\n"

# Muda o local de execução do job, fique atento se esse local faz sentido no seu projeto
cd "$SLURM_SUBMIT_DIR"

#pixi run diamond makedb --in /home/ttamaki/New/results/certo/pred_proteins.faa --db /home/ttamaki/New/results/certo/diamond_db.dmnd
pixi run diamond deepclust -d /home/ttamaki/New/results/diamond_db.dmnd -o /home/ttamaki/New/results2/deepclust_out.tsv --threads 8 --mutual-cover 80

echo -e "\n## Job finalizado em $(date +'%d-%m-%Y as %T') ##\n"
