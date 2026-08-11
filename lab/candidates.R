library(data.table)
library(tidyr)
library(Biostrings)

module_path <- "results/modules/exp_mod_r7.1_c0.42.tsv"
modules <- fread(module_path, header = TRUE)
modules

memberID_path <- "/home/ttamaki/New/results/member_cluster_id_annotated.tsv"
member_ID <- fread(memberID_path, header = TRUE)

member_ID[cluster_n == 238938]
cluster_genes <- unique(member_ID[cluster_n == 238938]$gene)
faa <- readAAStringSet("results/pred_proteins.faa")
faa

selected_sequences <- faa[sub(" #.*", "", names(faa)) %in% cluster_genes]
writeXStringSet(selected_sequences, "results_mutual_coverage/proteinas_filtradas.faa")


faa_filtrado <- faa[sub(" .*", "", names(faa)) %in% unique(member_ID[cluster_n == 238938][,.(gene,interpro_description)]$gene)]
lista <- unique(member_ID[cluster_n == 238938][,.(gene,interpro_description)]$gene)
writeXStringSet(faa_filtrado, "proteinas_filtradas.faa")
id(faa)
as.character(faa[[1]])

faa
getwd()
writeXStringSet(aaa, "lab/proteinas_peptidase.faa")
