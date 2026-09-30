#
# Script for split conjugative and non conjugative enriched cluster
#

library(data.table)

data <- fread("/home/ttamaki/Conjugative-Modules/results/enrichment/fisher/fisher_r7.6_c0.16_MPF.tsv", header = TRUE)
data_filtered <- data[padj < 0.05 & odds_ratio > 2]
length(unique(data_filtered$cluster_id))


select <- c("TraH","TraH-2","TraG","Type II/IV secretion", "TraK", "TraQ-like", "TraN-like", "TraD","TraP","TraI", "MobA/VirD2-like", "TraA", "TrbF", "TraC-like", "TraD/TraG", "TraU", "Type IV secretory", "TraG", "TraI", "TraX", "TraS", "MobL","TraG", "MobI", "TraK-like" ,"Tra", "relaxase", "Type IV secretion", "Conjugal", "Trb", "Mob", "Type-IV", "conjugative", "TcpE", "Type 4", "relaxation", "TcpC", "Pilus", "Conjugative" )
pattern <- paste0( "(?<![[:alnum:]_])(", paste(select, collapse = "|"), ")(?![[:alnum:]_])" )
known_elements <- data_filtered[ grepl( pattern, data_filtered$`InterPro Description`, ignore.case = TRUE, perl = TRUE ), ]
unknown_elements <- data_filtered[cluster_id %in% known_elements$cluster_id == FALSE]
data_filtered$
known_elements <- known_elements[order(odds_ratio, decreasing = TRUE)]
unknown_elements <- unknown_elements[order(odds_ratio, decreasing = TRUE)]

unique(data_filtered$`cluster_id`))
known <- unique(known_elements$`cluster_id`)
unknown <- unique(unknown_elements$`cluster_id`)

all_data <- fread("/home/ttamaki/Conjugative-Modules/results/member_cluster_table.tsv")
nrow(all_data[cluster_n %in% known])
nrow(all_data[cluster_n %in% unknown])


56 mil plasmideos
53332 161 genes 158 925 cluster sem singletons 905 cluster enriquecidos dos quais 495 sao conjugativos e 410 nao 