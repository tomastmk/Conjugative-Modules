"""
Este módulo recebe uma saída do meu programa de coocorrência
e reformata ele de maneira a separar os módulos de cluster em
módulos genes para cada replicon
"""


import argparse
import os
from collections import defaultdict

from tqdm import tqdm


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-cd", type=str, help="Communities directory path")
    parser.add_argument("-mi", type=str, help="Member Id path")
    args = parser.parse_args()

    # Dado um cluster, retorna os membros
    cluster_map_by_rep = defaultdict(lambda: defaultdict(list))

    # Dado um gene, retorna a qual cluster pertence
    member_to_cluster = {}

    with open(args.mi) as f:
        next(f)
        for line in f:
            m, c = line.strip().split("\t")
            c = int(c)

            replicon, gene = m.split("_")

            cluster_map_by_rep[c][replicon].append(int(gene))
            member_to_cluster[m] = c


    with os.scandir(args.cd) as entries:

        entries = list(entries)
        # Para cada arquivo de comunidade, cria um arquivo de módulo correspondente
        for entry in tqdm(entries, desc="Processing community files"):            
            modules = defaultdict(list)
            output_path = entry.path.replace("communities", "modules", 1)
            output_path = output_path.replace("com", "mod", 1)
            
            input_file = open(entry.path, "r", encoding="utf-8")
            next(input_file) # Pula header
                        
            output_file =  open(output_path, "w", encoding="utf-8")
            output_file.write("replicon\tgene\tcommunity_id\n")


            last_com_n = 0
            for line in input_file:

                print(line)
                replicon, gene_n, com_n = line.split("\t")
                com_n = int(com_n.strip("\n"))

                if com_n != last_com_n:
                    for plasmid, genes in modules.items():
                        for gene in genes:

                            line = f"{plasmid}\t{gene}\t{last_com_n}\n"

                            output_file.write(line)
                            
                    modules = defaultdict(list)
                    last_com_n = com_n

                gene = replicon+"_"+gene_n

                cluster_n: int = member_to_cluster[gene]
                cluster_genes: defaultdict[str,list] = cluster_map_by_rep[cluster_n]


                for rep in cluster_genes:
                    modules[rep].extend(cluster_genes[rep])

            for replicon, genes in modules.items():
                for gene in genes:

                    line = f"{replicon}\t{gene}\t{last_com_n}\n"

                    output_file.write(line)
    
            input_file.close()
            output_file.close()

if __name__ == "__main__":
    main()