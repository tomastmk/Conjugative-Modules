#
# Script for merging transmembrane proteins and positive selected proteins
#

library(data.table)

ab <- c("5598","5636","5764","6483","9323","10953","11973","13095","39561","42710")
cluster_annot <- fread("results/clustering/member_cluster_table_annotated.tsv")

filter <- cluster_annot[cluster_n %in% ab ]
unique(filter[,.(cluster_n,signature_description, interpro_description)])

cluster42710 <- fread("/home/ttamaki/Conjugative-Modules/results/clustering/member_cluster_table_annotated.tsv")
unique(cluster42710[cluster_n == 42710][signature_description == ''])

member_id <- fread('/home/ttamaki/Conjugative-Modules/results/clustering/member_cluster_table_annotated.tsv')
member_id
a <- member_id[grepl("TraN", interpro_description, ignore.case = FALSE)]
unique(a)
tem <- unique(a$gene)

writeLines(tem, "~/traN.txt")

enriched <- fread("/home/ttamaki/Conjugative-Modules/results/enrichment/fisher/fisher_r7.6_c0.16_MPF.tsv")
enriched_filtered1 <- enriched[padj<0.05]
enriched_filtered2 <- enriched[odds_ratio>2]
enriched_filtered3<- enriched[odds_ratio>2][padj<0.05]
sum(a %in% enriched_filtered1$cluster_id)
sum(a %in% enriched_filtered2$cluster_id)
sum(a %in% enriched_filtered3$cluster_id)
sum(a %in% enriched$cluster_id)
unique(enriched_filtered3$cluster_id)


enriched_filtered3
a[a %in% enriched_filtered3$cluster_id]

path_dir <- "/home/ttamaki/Conjugative-Modules/results/enrichment/fubar/selection_tsv"
path1 <- paste0(path_dir,"/cluster_5175.tsv")
data23 <- fread(path1)
data23[beta-alpha > 1]

temp[V1 %in% unique(member_id[cluster_n == 5175]$gene)]

temp <- fread("~/temp.txt")
temp

poposad <- unique(member_id[gene %in% temp$V1]$cluster_n)
poposad

temp %in% unique(member_id[cluster_n == 5175]$gene)


sum(unique(enriched_filtered3$cluster_id) %in% poposad)
