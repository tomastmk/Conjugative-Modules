library(data.table)

data <- fread("/home/ttamaki/Conjugative-Modules/results/enrichment/fisher/fisher_r7.6_c0.16_MPF.tsv")
pvalues <- data$padj
pvalues_less <- pvalues[pvalues<0.1]
hist(pvalues_less)
