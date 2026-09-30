#
# Script for reading metadata
#

library(data.table)-
library(ggplot2)

metadata_path <- "/home/ttamaki/shared/sequence_metadata.tsv"
data <- fread(metadata_path, sep = "\t", fill = TRUE, header = TRUE)
filtered_data <- data[,.(sequence_accession,scientific_name,tax_id)]
unique(filtered_data$scientific_name)
write.table(
    unique(filtered_data$tax_id),
    file = "/home/ttamaki/Conjugative-Modules/results/metadata/list_of_species.txt",
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
)

species <- fread("/home/ttamaki/Conjugative-Modules/results/metadata/taxonomy.tsv", header = TRUE)
species_filtered <- species[, .(Taxid, `Domain/Realm name`, `Kingdom name`, `Phylum name`, `Class name`, `Order name`, `Family name`, `Genus name`, `Species name`)]
species_filtered
filtered_data
merged_data <- merge(filtered_data, species_filtered, by.x = "tax_id", by.y = "Taxid")
merged_data

write.table(
    merged_data,
    file = "/home/ttamaki/Conjugative-Modules/results/metadata/merged_metadata.tsv",
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
)

merged_data <- fread("/home/ttamaki/Conjugative-Modules/results/metadata/merged_metadata.tsv", header = TRUE)
Domain_count <- table(merged_data$`Domain/Realm name`)
Domain_count
Kingdom_count <- table(merged_data$`Kingdom name`)
Kingdom_count
Phylum_count <- table(merged_data$`Phylum name`)
Phylum_count
Class_count <- table(merged_data$`Class name`)
Class_count
Order_count <- table(merged_data$`Order name`)
Order_count
Family_count <- table(merged_data$`Family name`)
Family_count[order(Family_count, decreasing = FALSE)]

merged_data[`Family name` %in% c("Borreliaceae","Bacillaceae","Enterobacteriaceae")]$`sequence_accession`

Genus_count <- table(merged_data$`Genus name`)
Genus_count
Species_count <- table(merged_data$`Species name`)
Species_count
length(Family_count)
different_counts <- list(
    Domain_count = length(Domain_count),
    Kingdom_count = length(Kingdom_count),
    Phylum_count = length(Phylum_count),
    Class_count = length(Class_count),
    Order_count = length(Order_count),
    Family_count = length(Family_count),
    Genus_count = length(Genus_count),
    Species_count = length(Species_count)
)
different_counts

df <- data.frame(
  rank = c(
    "Species",
    "Genus",
    "Family",
    "Order",
    "Class",
    "Phylum"
  ),
  n = c(
    length(Species_count),
    length(Genus_count),
    length(Family_count),
    length(Order_count),
    length(Class_count),
    length(Phylum_count)
  )
)

df$rank <- factor(
  df$rank,
  levels = rev(df$rank)
)

p <- ggplot(df, aes(x = n, y = rank)) +
  geom_col(width = 0.65) +
  geom_text(
    aes(label = format(n, big.mark = ",")),
    hjust = -0.15,
    size = 7,
    fontface = "bold"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 24) +
  theme(
    axis.text.y = element_text(
      size = 24,
      face = "bold"
    ),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    plot.margin = margin(10, 40, 10, 10)
  )

  p

library(dplyr)

df_phylum <- data.frame(
  phylum = names(Phylum_count),
  n = as.numeric(Phylum_count)
)

df_phylum <- df_phylum %>%
  arrange(desc(n))

# Manter os 8 maiores
df_phylum_plot <- df_phylum %>%
  slice_head(n = 8) %>%
  bind_rows(
    df_phylum %>%
      slice(-(1:8)) %>%
      summarise(
        phylum = "Other phyla",
        n = sum(n)
      )
  )
 df_phylum_plot <- df_phylum_plot[df_phylum_plot$phylum != "Pseudomonadota", ]
df_phylum_plot$phylum <- factor(
  df_phylum_plot$phylum,
  levels = rev(df_phylum_plot$phylum)
)

p_phylum <- ggplot(
  df_phylum_plot,
  aes(x = n, y = phylum)
) +
  geom_col(width = 0.65) +
  geom_text(
    aes(label = format(n, big.mark = ",")),
    hjust = -0.15,
    size = 6,
    fontface = "bold"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    x = "Number of plasmids",
    y = NULL
  ) +
  theme_classic(base_size = 22) +
  theme(
    axis.text.y = element_text(
      size = 18,
      face = "bold"
    ),
    axis.text.x = element_text(size = 17),
    axis.title.x = element_text(
      size = 21,
      face = "bold"
    )
  )

p_phylum
