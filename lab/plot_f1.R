library(data.table)
library(ggplot2)



data <- fread("/home/ttamaki/New/results2/performance/f1_t4ss_conjscan.csv", header = FALSE)
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


p <- ggplot(resultado,
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
ggsave(
    filename = "figs/f1_best_module.png",
    plot = p,
    width = 15, 
    height = 15, 
    units = "cm", 
    dpi = 100
    )

head(resultado[order(resultado$F1, decreasing = TRUE),], n=10)
#############3

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
