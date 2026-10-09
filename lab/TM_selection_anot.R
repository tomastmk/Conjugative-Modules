#
# Script for merging transmembrane proteins and positive selected proteins
#

library(data.table)

cluster_annot <- fread("results/clustering/member_cluster_table_annotated.tsv")

# Teste de enriquecimento
enriched <- fread("/home/ttamaki/Conjugative-Modules/results/enrichment/fisher/fisher_r7.6_c0.16_MPF.tsv")
enriched_filtered1 <- enriched[padj<0.05]
enriched_filtered2 <- enriched[odds_ratio>2]
enriched_filtered3<- enriched[odds_ratio>2][padj<0.05]
sum(a %in% enriched_filtered1$cluster_id)
sum(a %in% enriched_filtered2$cluster_id)
sum(a %in% enriched_filtered3$cluster_id)
sum(a %in% enriched$cluster_id)
unique(enriched_filtered3$cluster_id)

# Peptídeo Sinal
signalp_results <- fread("/home/ttamaki/Conjugative-Modules/results/enrichment/signalp_complete/prediction_results.txt")
SEC_proteins <- signalp_results[`LIPO(Sec/SPII)`>0.3]$`# ID`
SEC_clusters <- unique(cluster_annot[gene %in% SEC_proteins][,.(cluster_n)])
SEC_clusters[,SEC := TRUE]
SEC_clusters

# Seleção Positiva
fubar_results <- fread("/home/ttamaki/Conjugative-Modules/results/enrichment/fubar/fubar_adapted_output.tsv")
positive_selected <-unique(fubar_results[`Prob[alpha<beta]` > 0.5])
positive_selected_clusters <- unique(positive_selected$`print cluster_n`)
positive_selected_clusters

# Domínio Transmembrana
deep_TMHMM_results <- fread("/home/ttamaki/Conjugative-Modules/results/enrichment/membrane/transmembrane_proteins.tsv")
TM_proteins <- deep_TMHMM_results$Proteins
TM_clusters <- unique(cluster_annot[gene %in% TM_proteins][,.(cluster_n)])
TM_clusters[, TM := TRUE]
TM_clusters

### Intersec
intersec_clusters <- unique(cluster_annot[cluster_n %in% positive_selected_clusters][,.(cluster_n, interpro_description)])
intersec_clusters[,SEC := intersec_clusters$cluster_n %in% SEC_clusters$cluster_n]
intersec_clusters[,TM := intersec_clusters$cluster_n %in% TM_clusters$cluster_n]
intersec_clusters <- intersec_clusters[TM==TRUE | SEC==TRUE]
intersec_clusters

fwrite(intersec_clusters, "/home/ttamaki/Conjugative-Modules/results/enrichment/membrane_selected_proteins.tsv", sep = "\t")
fwrite(unique(intersec_clusters[,.(interpro_description)]), "/home/ttamaki/Conjugative-Modules/results/enrichment/membrane_selected_proteins_annotations.tsv", sep = "\t")
