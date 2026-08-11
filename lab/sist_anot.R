library(data.table)

df <- fread("/home/ttamaki/New/results_mutual_coverage/cooccurrence_out_member_id.tsv", header = TRUE)

annot <- fread(
  "/home/ttamaki/New/results/interproscan/genbank_interproscan_annotations.tsv",
  sep = "\t",
  header = FALSE,
  quote = ""
)

colnames(annot) <- c(
  "gene",
  "sequence_md5",
  "sequence_length",
  "analysis",
  "signature_accession",
  "signature_description",
  "start",
  "stop",
  "score",
  "status",
  "date",
  "interpro_accession",
  "interpro_description",
  "go_annotations",
  "pathways"
)

# =========================================================
# mantém apenas colunas relevantes
# =========================================================

annot <- annot[, .(
  gene,
  analysis,
  signature_accession,
  signature_description,
  interpro_accession,
  interpro_description
)]

# =========================================================
# left join
# =========================================================

memberId_annot <- merge(
  df,
  annot,
  by = "gene",
  all.x = TRUE
)

memberId_annot[, gene_n := tstrsplit(gene, "_", fixed = TRUE)[2]]

memberId_annot <- memberId_annot[order(cluster_n, gene_n)]
memberId_annot[, c("gene_n") := NULL]


fwrite(memberId_annot,
  file = "/home/ttamaki/New/results_mutual_coverage/member_cluster_id_annotated.tsv",
  sep="\t")


com <- fread("/home/ttamaki/Conj/data/processing/old_com/communities_filtered/jaccard_com_2.0.tsv", header = TRUE)
com[,gene := paste0(id,"_",gene_n)]

# Equivalente
community_annotations <- function(memberId_annot, annot_grep, random = TRUE, cluster_member = ""){

  relaxases <- memberId_annot[
    grepl(annot_grep, signature_description, ignore.case = TRUE) |
    grepl(annot_grep, interpro_description, ignore.case = TRUE)
  ]$gene

  if (random){
    n <- sample.int(length(relaxases), 1)
    cluster_member <- relaxases[n]
  }
  cluster_id <- unique(memberId_annot[memberId_annot$gene==cluster_member]$cluster_n)
  cluster <- memberId_annot[memberId_annot$cluster_n == cluster_id]$gene
  community_id <- com[com$gene %in% cluster]$community_n
  centroids <- com[com$community_n == community_id]$gene

  clusters_annot <- list()

  for (i in seq(from = 1, to = length(centroids))){
    centroid <- centroids[i]
    centroid_cluster_id <- unique(memberId_annot[gene == centroid]$cluster_n)
    x <- unique(
      memberId_annot[cluster_n == centroid_cluster_id]$interpro_description
    )

    x[x == "-"] <- NA

    centroid_cluster <- list(unique(x))
    clusters_annot[i] <- centroid_cluster
  }
  return(unique(clusters_annot))
}

a <- community_annotations(memberId_annot, "relaxase")
print(a)

###
relaxase_community_ids <- memberId_annot[memberId_annot$gene %in% relaxase_centroids]$cluster_n
relaxase_community <- memberId_annot[memberId_annot$cluster_n %in% relaxase_community_ids]
relaxase_community_annot <- relaxase_community$interpro_description


# separa por qualquer caractere que não seja letra ou número
palavras <- unlist(
  strsplit(
    tolower(relaxase_community_annot),
    "[ ,]+"
  )
)

palavras <- trimws(
  unlist(
    strsplit(relaxase_community_annot, ",")
  )
)
# remove vazios
palavras <- palavras[palavras != ""]


remover <- c(
  "domain",
  "family",
  "protein",
  "type",
  "n",
  "of",
  "like"
)

palavras <- palavras[!palavras %in% remover]
# conta frequências
freq <- sort(table(palavras), decreasing = TRUE)
head(freq,20)
table(palavras)

library(tibble)
library(dplyr)
data_frame(df) |> arrange(-f1)

unique(memberId_annot[cluster_n == 141633]$gene)

