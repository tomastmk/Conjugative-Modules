#!/bin/bash

# Esse nome vai aparecer na execução do squeue
#SBATCH --job-name=trees

# Esse é o número de processos
#SBATCH --ntasks=1

# Essa é a memória que seu processo pode usar
#SBATCH --mem=80GB

# Esse é o número de cpus
#SBATCH --cpus-per-task=8

# Essa é a estimativa de tempo. NOTE que seu processo será encerrado após esse tempo limite
#SBATCH --time=70:00:00

# Esses são os arquivos de saída do StdOut e do StdErr, respectivamente
#SBATCH --output=/home/ttamaki/Conjugative-Modules/out/%j_tree.out
#SBATCH --error=/home/ttamaki/Conjugative-Modules/err/%j_tree.err

# Para exportar seu ambiente pixi
export PATH="/home/ttamaki/Conjugative-Modules/.pixi/bin:$PATH"

echo -e "\n## Job iniciado em $(date +'%d-%m-%Y as %T') ##\n"

# Muda o local de execução do job, fique atento se esse local faz sentido no seu projeto
cd "$SLURM_SUBMIT_DIR"

pixi run /home/ttamaki/Conjugative-Modules/src/tree.sh \
  /home/ttamaki/Conjugative-Modules/results/enrichment/fubar/trees \
  /home/ttamaki/Conjugative-Modules/results/enrichment/fubar/protein_filtered

echo -e "\n## Job finalizado em $(date +'%d-%m-%Y as %T') ##\n"
