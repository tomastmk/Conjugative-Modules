library(data.table)
library(dplyr)

df <- fread("/home/ttamaki/Conj/data/processing/graphs/graph.csv",header = TRUE)
alpha <- fread("/home/ttamaki/Conj/data/processing/graphs/graph_filtered_alpha.tsv", header = TRUE)
df1 <- df[,c(1,2,7,8,9)]
df_num <- df[,c(7,8,9)]


df_final <- merge(
  df1,
  alpha,
  by = c("gene1", "gene2"),
  all = TRUE
)

df_final[is.na(df_final)] <- 0
df_final_num <- df_final[,c("ochiai","exp")]
df_final_num[,"alpha"] <- df_final_num[,"alpha"]/10

## analise sobre a matriz de covariancias
CP1 <- princomp(df_final_num)
summary(CP1, loadings=TRUE)
hist(CP1$scores[,1])


hist(df$alpha)
mean <- mean(df$alpha)
sigma <- sd(df$alpha)


max <- max(df$alpha)
min <- min(df$alpha)

df$norm <- (df$alpha-min)/(max-min)
df$norm2 <- (df$alpha-mean)/sigma
par(mfrow=c(2,1))

hist(df$alpha)
hist(df$norm)
hist(df$norm2)
results <- df[,c("gene1","gene2","norm")]
colnames(results) <- c("gene1","gene2","alpha")

output <- "aa.tsv"
write.table(results,
          file = output,
          sep = "\t",
          row.names = FALSE,
          quote = FALSE) 
