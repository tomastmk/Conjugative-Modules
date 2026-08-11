library(data.table)

data <- fread("/home/ttamaki/New/results_mutual_coverage/fisher_r7.6_c0.16_MPF.tsv", header = TRUE)
data <- r2
data[padj<0.05,]
fwrite(filter7[padj<0.05], "/home/ttamaki/New/results_mutual_coverage/fisher_r7.6_c0.16_MPF_a0.05_filtered.tsv", sep = "\t")


## Filtering 
data <- fread("/home/ttamaki/New/results/performance/teste_fisher_r7.1_c0.42_MPF_a0.05.tsv", header = TRUE)
data
filter1 <- data[!grepl("Tra",data$`InterPro Description`),]
filter2 <- filter1[!grepl("relaxase",filter1$`InterPro Description`),]
filter3 <- filter2[!grepl("Type IV",filter2$`InterPro Description`),]
filter4 <- filter3[!grepl("Conjugal",filter3$`InterPro Description`),]
filter5 <- filter4[!grepl("Trb",filter4$`InterPro Description`),]
filter6 <- filter5[!grepl("Mob",filter5$`InterPro Description`),]
filter7 <- filter6[!grepl("Type-IV",filter6$`InterPro Description`),]


filter7$`InterPro Description`
data <- fread("/home/ttamaki/New/results/performance/fisher_r7.6_c0.16_MPF_filtered.tsv", header = TRUE)
data

