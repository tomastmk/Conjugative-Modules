library(data.table)

df <- fread("/home/ttamaki/Conj/data/results/new/jaccard_mod_2.0.tsv", header = TRUE)
setDT(df)

df[, id := rleid(replicon)]

df[, gene_id := paste0(replicon, "_", gene)]
    
df[, c("replicon") := NULL]
df[, c("module_id") := NULL]

print(df)

annot <- fread(
  "/home/ttamaki/Conj/data/results_interproscan/genbank_interproscan_annotations.tsv",
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

resultado <- merge(
  df,
  annot,
  by.x = "gene_id",
  by.y = "gene",
  all.x = TRUE
)
resultado <- resultado[order(id,gene)]
df <- df[order(id,gene)]
resultado[, c("gene") := NULL]

print(resultado)

write.table(resultado, "/home/ttamaki/Conj/data/results/new/annotated_jaccard_mod_2.0.tsv", sep = "\t", row.names = FALSE)
