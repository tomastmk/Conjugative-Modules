library(data.table)
library(ggridges)
library(ggplot2)

# lista:
# resolucao | valor
df <- fread("/home/ttamaki/Conj/data/results/new/filtered_metrics.tsv", header = TRUE)
head(df)

### Ridge graph

setDT(df)
metrica_escolhida <- "jaccard"
res_keep <- sort(unique(df$resolucao))[
  seq(1, length(unique(df$resolucao)), by = 4)
]

df_sub <- df[
  metrica == metrica_escolhida &
  resolucao %in% res_keep
]

# ==========================================
# Plot sobreposto
# ==========================================

ggplot() +

  geom_density_ridges(
    data = transform(df_sub, metric = "F1"),
    aes(
      x = f1,
      y = factor(resolucao),
      fill = metric
    ),
    alpha = 0.3,
    scale = 1
  ) +

  geom_density_ridges(
    data = transform(df_sub, metric = "Recall"),
    aes(
      x = recall,
      y = factor(resolucao),
      fill = metric
    ),
    alpha = 1,
    scale = 1
  ) +

  geom_density_ridges(
    data = transform(df_sub, metric = "Precision"),
    aes(
      x = precision,
      y = factor(resolucao),
      fill = metric
    ),
    alpha = 0.3,
    scale = 1
  ) +

  scale_fill_manual(
    values = c(
      "F1" = "#D81B60",
      "Recall" = "#1E88E5",
      "Precision" = "#FFC107"
    )
  ) +

  labs(
    title = metrica_escolhida,
    x = "Valor",
    y = "Resolução",
    fill = "Métrica"
  ) +

  theme_minimal()


### Ridge plot delta len

ggplot(
  df_sub,
  aes(
    x = delta_len,
    y = factor(resolucao),
    fill = after_stat(x)
  )
) +
  geom_density_ridges_gradient(
    scale = 1,
    rel_min_height = 0.01
  ) +
  scale_fill_viridis_c() +
  coord_cartesian(xlim = c(-250, 20)) +
  theme_minimal()


### Ridge plot f1, recall e precision separados

plot_df <- melt(
  df_sub,
  id.vars = c("metrica", "resolucao"),
  measure.vars = c("f1", "recall", "precision"),
  variable.name = "metric",
  value.name = "value"
)

ggplot(
  plot_df,
  aes(
    x = value,
    y = factor(resolucao),
    fill = after_stat(density)
  )
) +
  geom_density_ridges_gradient(
    scale = 1,
    rel_min_height = 0.01
  ) +
  facet_wrap(~ metric, nrow = 1) +
  scale_fill_viridis_c() +
  labs(
    x = "Valor",
    y = "Resolução",
    fill = "Densidade"
  ) +
  theme_minimal()



## heatmap
metric <- "precision"

breaks <- seq(
  0,
  1,
  length.out = 50
)

heat_df <- df[, {

  h <- hist(
    get(metric),
    breaks = breaks,
    plot = FALSE
  )

  .(
    bin = h$mids,
    freq = h$counts
  )

}, by = .(metrica, resolucao)]

ggplot(
  heat_df,
  aes(
    x = resolucao,
    y = bin,
    fill = freq
  )
) +
  geom_tile() +
  scale_fill_viridis_c(option = "C") +
  facet_wrap(~ metrica) +
  labs(
    x = "Resolução",
    y = metric,
    fill = "Frequência"
  ) +
  theme_minimal()
    


### mediana

library(data.table)
library(ggplot2)

setDT(df)
# ==========================================
# Métrica escolhida
# ==========================================

metric <- "recall"
# ou:
# "recall"
# "precision"

# ==========================================
# Resumo por resolução
# ==========================================

summary_df <- df[, .(
  mean   = mean(get(metric), na.rm = TRUE),
  median = median(get(metric), na.rm = TRUE),
  q25    = quantile(get(metric), 0.25, na.rm = TRUE),
  q75    = quantile(get(metric), 0.75, na.rm = TRUE)
), by = .(metrica, resolucao)]

# ==========================================
# Plot
# ==========================================
ggplot(summary_df, aes(x = resolucao)) +

  # faixa interquartil
  geom_ribbon(
    aes(ymin = q25, ymax = q75),
    fill = "grey70",
    alpha = 0.4
  ) +

  # média
  geom_point(
    aes(y = mean, color = "Média"),
    size = 2
  ) +

  # mediana
  geom_point(
    aes(y = median, color = "Mediana"),
    size = 2
  ) +

  facet_wrap(~ metrica) +

  labs(
    x = "Resolução",
    y = metric,
    color = ""
  ) +

  theme_minimal()
