#
# Script for plotting f1-scores
#

library(data.table)
library(ggplot2)

data <- fread("/home/ttamaki/Conjugative-Modules/results/performance/f1_t4ss_conjscan.csv", header = FALSE)
setnames(data, c("File","TP","FP","FN"))
data
resultado <- data[
    ,
    lapply(.SD, sum),
    by = File,
    .SDcols = c("TP", "FP", "FN")
]
resultado
resultado[
    ,
    `:=`(
        metric = sub("_mod.*", "", File),
        resolution = as.numeric(sub(".*_r([0-9.]+)_c.*", "\\1", File)),
        cut = as.numeric(sub(".*_c([0-9.]+)\\.tsv$", "\\1", File))
    )
]
resultado <- resultado[,.(metric,resolution,cut,TP,FP,FN)]
resultado[,F1 := 2*TP/(2*TP+FP+FN)]
sort(resultado$F1)
resultado[F1 >= 0.65]
to_plot <- resultado[metric == "exp"]

data_filtered <- data[File == "exp_mod_r7.6_c0.16.tsv"]
data_filtered[, F1 := 2*TP/(2*TP+FP+FN)]
hist(data_filtered$F1)
expo <- data[metric == "exp" & resolution == 7.6 & cut == 0.16]
data_filtered
library(ggplot2)

data_filtered
a <- data_filtered$F1
a
p <- ggplot(data_filtered, aes(x = F1)) +

  geom_histogram(
    bins = 30,
    fill = "#183757"
  ) +

  theme_classic(base_size = 24) +

  labs(
    x = "F1-score (Similarity)",
    y = "Count"
  ) +

  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 60
    ),

    axis.title = element_text(
      face = "bold",
      size = 50
    ),

    axis.text = element_text(
      size = 40
    )
  )

p

ggsave(
    filename = "export.png",
    plot = p,
    width = 45, 
    height = 25, 
    units = "cm", 
    dpi = 300
    )
p <- ggplot(to_plot,
    aes(x = resolution, y = cut, fill = F1)) +
    geom_tile() +
    facet_wrap(~metric) +
    scale_fill_viridis_c(name = "F1") +
    labs(
        x = "Resolution",
        y = "Cutoff"
    ) +
    theme_bw()
p

p <- ggplot(
    to_plot,
    aes(x = resolution, y = cut, fill = F1)
) +
    geom_tile(
        color = "white",
        linewidth = 0.1
    ) +
    facet_wrap(~metric, nrow = 1) +
    scale_fill_gradientn(
    name = "F1",
    colours = c(
        "#FFF7BC",
        "#FEC44F",
        "#FE9929",
        "#EC7014",
        "#CC4C02",
        "#993404"
    ),
    values = scales::rescale(c(
        min(to_plot$F1, na.rm = TRUE),
        0.4,
        0.50,
        0.60,
        0.62,
        0.66,
        max(to_plot$F1, na.rm = TRUE)
    )),
    limits = range(to_plot$F1, na.rm = TRUE)
) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
        x = "Resolution",
        y = "Cutoff"
    ) +
    theme_minimal(base_size = 30) +
    theme(
        panel.grid = element_blank(),
        panel.border = element_rect(
            color = "grey70",
            fill = NA,
            linewidth = 0.5
        ),
        strip.background = element_rect(
            fill = "grey95",
            color = NA
        ),
        strip.text = element_text(
            face = "bold",
            size = 11
        ),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold", size = 15)
    )

p
1)

ggsave(
    filename = "~/send.png",
    plot = p,
    width = 45, 
    height = 30, 
    units = "cm", 
    dpi = 100
    )

head(resultado[order(resultado$F1, decreasing = TRUE),], n=10)
#############
# Ploting histogram of weights
grafo <- fread("/home/ttamaki/Conjugative-Modules/results/cooccurrence_out_graph.tsv")
grafo

jaccard <- grafo$jaccard
ochiai <- grafo$ochiai
expm <- grafo$exp

# Plot completo
library(ggplot2)
library(tidyr)

dados_hist_all <- data.frame(
  metrica = c(
    rep("Jaccard", length(as.numeric(jaccard))),
    rep("Otsuka-Ochiai", length(as.numeric(ochiai))),
    rep("Exp", length(as.numeric(expm)))
  ),
  valor = c(
    as.numeric(jaccard),
    as.numeric(ochiai),
    as.numeric(expm)
  )
) |>
  subset(!is.na(valor))

grafico_all <- ggplot(dados_hist_all, aes(x = valor, fill = metrica)) +
  geom_histogram(bins = 30, color = "white", na.rm = TRUE) +
  facet_wrap(~ metrica, nrow = 1, scales = "free") +
  guides(fill = "none") +
  labs(x = "Valor", y = "Frequência", title = "Distribuição dos pesos por métrica") +
  theme_minimal()

ggsave(
  "/home/ttamaki/Conjugative-Modules/figs/relatorio_pub/histogramas_all.png",
  plot = grafico_all,
  width = 12,
  height = 4,
  dpi = 300
)

# Plot somente valores maiores que 0.05
jaccard_filtro <- jaccard[jaccard > 0.05]
ochiai_filtro <- ochiai[ochiai > 0.05]
expm_filtro <- expm[expm > 0.05]

dados_hist_filtro <- data.frame(
  metrica = c(
    rep("Jaccard", length(as.numeric(jaccard_filtro))),
    rep("Otsuka-Ochiai", length(as.numeric(ochiai_filtro))),
    rep("Exp", length(as.numeric(expm_filtro)))
  ),
  valor = c(
    as.numeric(jaccard_filtro),
    as.numeric(ochiai_filtro),
    as.numeric(expm_filtro)
  )
) |>
  subset(!is.na(valor))

grafico_filtro <- ggplot(dados_hist_filtro, aes(x = valor, fill = metrica)) +
  geom_histogram(bins = 30, color = "white", na.rm = TRUE) +
  facet_wrap(~ metrica, nrow = 1, scales = "free") +
  guides(fill = "none") +
  labs(x = "Valor", y = "Frequência", title = "Distribuição dos pesos com valor > 0.05") +
  theme_minimal()

ggsave(
  "/home/ttamaki/Conjugative-Modules/figs/relatorio_pub/histogramas0.05.png",
  plot = grafico_filtro,
  width = 12,
  height = 4,
  dpi = 300
)

#############
# Ploting number of communities acording to the parameters
count_com <- fread("/home/ttamaki/New/results_new/communities_count_test.tsv", header = TRUE)
p <- ggplot(count_com,
    aes(x = resolution, y = cut, fill = number_of_clusters)) +
    geom_tile() +
    facet_wrap(~metric) +
    scale_fill_viridis_c(name = "Quantidade de comunidades") +
    labs(
        x = "Resolution",
        y = "Cutoff"
    ) +
    theme_bw()
p
ggsave(
    filename = "/home/ttamaki/New/figs/number_communities.png",
    plot = p,
    width = 30, 
    height = 15, 
    units = "cm", 
    dpi = 100
    )

p <- ggplot(count_com,
    aes(x = resolution, y = cut, fill = log(number_of_clusters))) +
    geom_tile() +
    facet_wrap(~metric) +
    scale_fill_viridis_c(name = "Quantidade de comunidades (log)") +
    labs(
        x = "Resolution",
        y = "Cutoff"
    ) +
    theme_bw()
p
ggsave(
    filename = "/home/ttamaki/New/figs/number_communities_log.png",
    plot = p,
    width = 30, 
    height = 15, 
    units = "cm", 
    dpi = 100
    )
