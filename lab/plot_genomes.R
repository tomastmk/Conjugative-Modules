#
# Script for plotting genomes
#


library(data.table)
library(gggenomes)
library(Biostrings)
library(dplyr)

### All predicted proteins
prediction_headers <- fread(
  cmd = "grep '^>' '/home/ttamaki/Conjugative-Modules/results/pred_proteins.faa'",
  header = FALSE
)$V1


prediction_headers <- sub("^>", "", prediction_headers)

### Genes start and end
all_genes_data <- data.table(header = prediction_headers)
all_genes_data[, c("gene", "start", "end", "strand","sla?") :=
    tstrsplit(header, " # ", fixed = TRUE)]
all_genes_data <- all_genes_data[,c("gene", "start", "end", "strand")]
all_genes_data[, c("replicon","gene_n") := tstrsplit(all_genes_data$gene, "_")]
all_genes_data <- all_genes_data[,.(replicon,gene,start,end,strand)]
all_genes_data$start <- as.numeric(all_genes_data$start)
all_genes_data$end <- as.numeric(all_genes_data$end)
all_genes_data


### My modules
my_modules <- fread("/home/ttamaki/Conjugative-Modules/results/modules/exp_mod_r7.6_c0.16.tsv")
setnames(my_modules,old = "gene", new = "gene_n")
my_modules[, gene := paste0(my_modules$replicon,"_", my_modules$gene_n)]
my_modules_genes <- my_modules %>%
  left_join(
    select(all_genes_data, gene, start, end, strand),
    by = "gene"
  )
my_modules_genes <- my_modules_genes[,.(replicon, gene, community_id, start, end, strand)]
my_modules_genes
a <- my_modules[replicon == "CP124978.1"]
a[order(gene_n, decreasing = FALSE)]

### Conjscan modules
T4SS_conjscan <- fread("/home/ttamaki/Conjugative-Modules/results/CONJScan/systems_dt_T4SS.tsv", header=TRUE)
T4SS_conjscan_genes <- left_join(T4SS_conjscan, all_genes_data, by = "gene")
setnames(T4SS_conjscan_genes, old = "plasmid", new = "replicon")
T4SS_conjscan_genes[,community_id := "CONJScan"]
T4SS_conjscan_genes <- T4SS_conjscan_genes[,.(replicon, gene, community_id, start, end, strand)]
T4SS_conjscan_genes

T4SS_conjscan_genes

### Preparing plot
# Selecting Plasmid


reorder_conjscan <- function(x) {
  others <- setdiff(unique(x), "CONJScan")
  levels <- append(others, "CONJScan", after = 2)
  factor(x, levels = levels)
}

plot_modules <- function(plasmid) {

  ## Communities containing at least one CONJScan gene
  conjugative_communities <- unique(
    my_modules_genes[
      gene %in% T4SS_conjscan_genes[replicon == plasmid]$gene
    ]$community_id
  )

  ## Data to plot
  modules_to_plot <- rbind(
    my_modules_genes[community_id %in% conjugative_communities],
    T4SS_conjscan_genes
  )[replicon == plasmid]

  modules_to_plot[, `:=`(
    start = as.numeric(start),
    end   = as.numeric(end)
  )]

  ## Coordinates
  min_start <- min(modules_to_plot$start)
  max_end   <- max(modules_to_plot$end)
  plot_length <- max_end - min_start

  modules_to_plot[, `:=`(
    start = start - min_start,
    end   = end - min_start
  )]

  ## Shared genes
  shared <- modules_to_plot[, if (.N > 1) .SD, by = gene]

  shared_genes <- modules_to_plot[
    gene %in% unique(shared$gene) &
      community_id != "CONJScan"
  ]

  modules_to_plot[, shared := "Other"]

  modules_to_plot[
    shared_genes,
    on = "gene",
    shared := i.community_id
  ]

  ## Sequences
  s0 <- tibble(
    seq_id = c("ref", unique(modules_to_plot$community_id)),
    length = plot_length
  )

  ## Genes
  g0 <- tibble(
    seq_id = modules_to_plot$community_id,
    start  = modules_to_plot$start,
    end    = modules_to_plot$end,
    shared = modules_to_plot$shared
  )

  ## Reference plasmid
  fill_genes <- all_genes_data[
    replicon == plasmid &
      start >= min_start &
      end <= max_end
  ]

  g0 <- bind_rows(
    g0,
    tibble(
      seq_id = "ref",
      start = fill_genes$start - min_start,
      end   = fill_genes$end - min_start,
      shared = "Other"
    )
  )

  ## Links
  shared_conjscan <- shared[community_id == "CONJScan"]
  shared_modules  <- shared[community_id != "CONJScan"]

  l0 <- tibble(
    seq_id  = shared_conjscan$community_id,
    start   = shared_conjscan$start,
    end     = shared_conjscan$end,
    seq_id2 = shared_modules$community_id,
    start2  = shared_conjscan$start,
    end2    = shared_conjscan$end
  )

  ## Ordering
  order_tracks <- function(df) {
    df %>%
      mutate(seq_id = reorder_conjscan(seq_id)) %>%
      arrange(seq_id) %>%
      mutate(seq_id = as.character(seq_id))
  }

  s0 <- order_tracks(s0)
  g0 <- order_tracks(g0)
  l0 <- order_tracks(l0)

  ## Colors
  cols <- setNames(
    hcl.colors(length(unique(modules_to_plot$shared)), "Dark 3"),
    sort(unique(modules_to_plot$shared))
  )
  cols["Other"] <- "grey80"

  ## Plot
  p <- gggenomes(
    genes = g0,
    seqs = s0,
    links = l0
  ) +
    geom_seq(size = 1) +
    geom_gene(size = 5, aes(fill = shared)) +
    geom_link(aes(fill = seq_id2)) +
    scale_fill_manual(values = cols) +
    labs(
      title = paste("Plasmid:", plasmid),
      subtitle = "Plasmids and CONJScan module",
      x = "Position",
      y = NULL,
      fill = "Community"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )
  ggsave(
    filename = paste0("/home/ttamaki/Conjugative-Modules/figs/relatorio_pub/modulos", plasmid, "_modules_plot.png"), 
    plot = p, 
    width = 8, 
    height = 12, 
    units = "in", 
    dpi = 300
  )
  
}


plasmids <- sample(my_modules_genes[gene %in% T4SS_conjscan_genes$gene]$replicon, 100)


for (plasmid in plasmids){
  plot_modules(plasmid)
}



library(ggplot2)
library(data.table)
df <- data.frame(
  distance = 1:6,
  similarity = exp(1 - 1:6)
)
df
p <- ggplot(df, aes(x = distance, y = similarity)) +

  geom_line(linewidth = 1.5) +

  geom_point(size = 4) +

  labs(
    x = "Gene distance, d(x,y)",
    y = "Similarity",
    title = expression(S(x,y) == e^{1-d(x,y)})
  ) +

  theme_classic(base_size = 24) +

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
  ) +

  scale_x_continuous(
    breaks = 1:8
  ) +

  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)
  )

p

ggsave(
  filename = "export2.png",
  plot = p,
  width = 45,
  height = 30,
  units = "cm",
  dpi = 300
)
