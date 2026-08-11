#!/bin/bash

# Esse nome vai aparecer na execução do squeue
#SBATCH --job-name=com

# Esse é o número de processos
#SBATCH --ntasks=1

# Essa é a memória que seu processo pode usar
#SBATCH --mem=160GB

# Esse é o número de cpus
#SBATCH --cpus-per-task=12

# Essa é a estimativa de tempo. NOTE que seu processo será encerrado após esse tempo limite
#SBATCH --time=190:00:00

# Esses são os arquivos de saída do StdOut e do StdErr, respectivamente
#SBATCH --output=/home/ttamaki/New/out/%j_communities.out
#SBATCH --error=/home/ttamaki/New/err/%j_communities.err

# Para exportar seu ambiente pixi
export PATH="/home/ttamaki/Conjugative-Modules/.pixi/bin:$PATH"

echo -e "\n## Job iniciado em $(date +'%d-%m-%Y as %T') ##\n"

# Muda o local de execução do job, fique atento se esse local faz sentido no seu projeto
cd "$SLURM_SUBMIT_DIR"

pixi run python /home/ttamaki/New/src/find_communities.py \
    -i /home/ttamaki/New/results2/cooccurrence_out_graph.tsv \
    -o/home/ttamaki/New/results2/communities \
    -mi /home/ttamaki/New/results2/cooccurrence_out_member_id.tsv

echo -e "\n## Job finalizado em $(date +'%d-%m-%Y as %T') ##\n"
