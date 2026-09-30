#!/bin/bash

# Esse nome vai aparecer na execução do squeue
#SBATCH --job-name=dp

# Esse é o número de processos
#SBATCH --ntasks=1

# Essa é a memória que seu processo pode usar
#SBATCH --mem=180GB

# Esse é o número de cpus
#SBATCH --cpus-per-task=8

# Essa é a estimativa de tempo. NOTE que seu processo será encerrado após esse tempo limite
#SBATCH --time=170:00:00

# Esses são os arquivos de saída do StdOut e do StdErr, respectivamente
#SBATCH --output=/home/ttamaki/Conjugative-Modules/out/%j_deep.out
#SBATCH --error=/home/ttamaki/Conjugative-Modules/err/%j_deep.err

# Para exportar seu ambiente pixi
export PATH="/home/ttamaki/DeepTMHMM-Academic-License-v1.0/.pixi/bin:$PATH"

echo -e "\n## Job iniciado em $(date +'%d-%m-%Y as %T') ##\n"

# Muda o local de execução do job, fique atento se esse local faz sentido no seu projeto
cd "$SLURM_SUBMIT_DIR"

source /home/ttamaki/DeepTMHMM-venv/bin/activate

pixi run  python3 /home/ttamaki/DeepTMHMM-Academic-License-v1.0/predict.py \
  --fasta /home/ttamaki/Conjugative-Modules/results/enrichment/alignments/protein_filtered/c/primeiras_proteinas.faa \
  --outputdir /home/ttamaki/Conjugative-Modules/results/enrichment/membrane


echo -e "\n## Job finalizado em $(date +'%d-%m-%Y as %T') ##\n"
