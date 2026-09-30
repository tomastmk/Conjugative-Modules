#
# Script para fazer o plot e separar as proteínas com beta-alpha > 1 e P(beta > alpha) > 0.65
#


library(data.table)
library(ggplot2)
library(tidyr)
library(patchwork)

dir_path <- "/home/ttamaki/Conjugative-Modules/results/enrichment/fubar/selection_tsv"
output_path <- "/home/ttamaki/Conjugative-Modules/results/enrichment/fubar/beta-alpha_hist"


for (file in list.files(dir_path, full.names = FALSE)) {
    full_name <- file.path(dir_path, file)
    data <- fread(full_name)

    plot_data <- data |>
    pivot_longer(
        cols = c(`beta-alpha`, `Prob[alpha<beta]`),
        names_to = "variable",
        values_to = "value"
    )

    a <- ggplot(data, aes(x = `beta-alpha`)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "white") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    coord_cartesian(xlim = c(-10, 10)) +
    theme_classic()

    b <- ggplot(data, aes(x = `Prob[alpha<beta]`)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "white") +
    theme_classic()

    c <- a / b

    ggsave(
        filename = paste0(output_path, "/", file, "_histogram.png"),
        plot = c,
        width = 8,
        height = 6,
        dpi = 300
    )
}

data
y <- 0
total <- 0
files_selected <- character()

for (file in list.files(dir_path, full.names = FALSE)) {
    full_name <- file.path(dir_path, file)
    data <- fread(full_name)
    

    if (nrow(data[`Prob[alpha>beta]` > 0.65 & `beta-alpha` > 1]) > 1){
        files_selected <- c(files_selected, file)
    }
}


files_selected
