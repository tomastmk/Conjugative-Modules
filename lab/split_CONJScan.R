#
# Script para separar CONJScan em MOB de T4SS
#

library(data.table)
library(stringr)


conjscan_path <- "/home/ttamaki/New/results/CONJScan/read_best_solution.tsv"
conjscan_data <- fread(conjscan_path)
conjscan <- conjscan_data[,.(replicon, hit_id, gene_name)]

setnames(conjscan, old = c("replicon","hit_id","gene_name"), new = c("plasmid","gene","conj_annot"), skip_absent = TRUE)
conjscan[,c("gene_n") := tstrsplit(conjscan$gene,"_")[2]]
conjscan



MOB_systems <- conjscan[
    str_detect(
        conj_annot,
        regex("MOB|VirD4|TraD|TrwB|TcpA|PcfC|T4CP", ignore_case = TRUE)
    ),
]
T4SS_systems <- conjscan[!conj_annot %in% MOB_systems$conj_annot & !conj_annot =="-"]

T4SS_systems <- T4SS_systems[,.(gene,plasmid,gene_n,conj_annot)]
MOB_systems <- MOB_systems[,.(gene,plasmid,gene_n,conj_annot)]
fwrite(MOB_systems,"/home/ttamaki/New/results/certo/CONJScan/systems_dt_MOB.tsv", sep="\t")
fwrite(T4SS_systems,"/home/ttamaki/New/results/certo/CONJScan/systems_dt_T4SS.tsv", sep = "\t")

