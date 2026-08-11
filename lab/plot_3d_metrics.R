library(data.table)
library(ggplot2)

path1 <-"/home/ttamaki/New/results/var_results/metricas_medias.tsv"
path2 <-"/home/ttamaki/New/results/var_results/metricas2_medias.tsv"
path3 <-"/home/ttamaki/New/results/var_results/metricas3_medias.tsv"
path4 <-"/home/ttamaki/New/results/var_results/metricas4_medias.tsv"
path5 <-"/home/ttamaki/New/results/var_results/metricas5_medias.tsv"


df1 <- fread(path1, header = TRUE, sep = "\t")
df2 <- fread(path2, header = TRUE, sep = "\t")
df3 <- fread(path3, header = TRUE, sep = "\t")
df4 <- fread(path4, header = TRUE, sep = "\t")
df5 <- fread(path5, header = TRUE, sep = "\t")

df <- rbindlist(list(df1,df2,df3,df4,df5))

df[, c("metrica", "tipo", "resolucao", "corte") :=
        tstrsplit(arquivo, "_", fixed = TRUE)]

df[, resolucao := as.numeric(sub("^r", "", resolucao))]
df[, corte := as.numeric(sub("\\.tsv$", "", sub("^c", "", corte)))]


df[,c("tipo","arquivo") := NULL]

setorder(df, -f1)
# nao continuo
ggplot(df[df$resolucao<10 & df$corte<0.5], aes(x = resolucao, y = corte, fill = f1)) +
  geom_tile() +
  facet_wrap(~ metrica) +
  scale_fill_viridis_c(option = "turbo")

par(mfrow=c(1,3))
hist(df[df$metrica =="exp"]$f1)
hist(df[df$metrica =="ochiai"]$f1)
hist(df[df$metrica =="jaccard"]$f1)



par(mfrow=c(1,3))
graph <- fread("/home/ttamaki/New/results/graphs/graph_filtered.tsv")

exp <- data.frame(exp = unlist(graph[graph$exp>0.05]$exp))
jaccard <- data.frame(jaccard = unlist(graph[graph$jaccard>0.05]$jaccard))
ochiai <- data.frame(ochiai = unlist(graph[graph$ochiai>0.05]$ochiai))



ggplot(exp, aes(x = exp)) +
  geom_histogram(
    bins = 30,
    fill = "#4C78A8",
    color = "white",
    linewidth = 0.3
  ) +
  labs(
    x = "Similaridade (exp)",
    y = "Frequência",
    title = expression("")
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

ggplot(jaccard, aes(x = jaccard)) +
  geom_histogram(
    bins = 30,
    fill = "#4C78A8",
    color = "white",
    linewidth = 0.3
  ) +
  labs(
    x = "Similaridade (jaccard)",
    y = "Frequência",
    title = expression("")
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

ggplot(ochiai, aes(x = ochiai)) +
  geom_histogram(
    bins = 30,
    fill = "#4C78A8",
    color = "white",
    linewidth = 0.3
  ) +
  labs(
    x = "Similaridade (ochiai)",
    y = "Frequência",
    title = expression("")
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )


# nao continuo
## f1
ggplot(df[df$resolucao <2 & df$corte<0.5], aes(x = resolucao, y = corte, fill = f1)) +
  geom_tile() +
  facet_wrap(~ metrica) +
  scale_fill_viridis_c(option = "G")



## precision
ggplot(df[df$resolucao <3 & df$corte<0.5], aes(x = resolucao, y = corte, fill = precision)) +
  geom_tile() +
  facet_wrap(~ metrica) +
  scale_fill_viridis_c(option = "G")

## recall
ggplot(df[df$resolucao < 3 & df$corte<0.5], aes(x = resolucao, y = corte, fill = recall)) +
  geom_tile() +
  facet_wrap(~ metrica) +
  scale_fill_viridis_c(option = "G")




# continuo
ggplot(df[df$resolucao < 2 & df$corte < 0.3, ],
       aes(x = resolucao, y = corte, z = f1)) +
  geom_contour_filled() +
  facet_wrap(~ metrica) +
  scale_fill_viridis_d(option = "G")

ggplot(df[df$resolucao < 1 & df$corte < 0.3, ],
       aes(x = resolucao, y = corte, z = f1)) +
  geom_contour_filled() +
  geom_contour(color = "black", alpha = 0.2) +
  facet_wrap(~ metrica)+
  scale_fill_viridis_d(option = "turbo")
