library(Matrix)


args <- commandArgs(trailingOnly = TRUE)
INPUT_PATH <- args[1] #"/home/ttamaki/Conjugative-Modules/l_data/processing/cluster_sequence.tsv"
OUTPUT_PATH <- args[2] #"/home/ttamaki/Conjugative-Modules/l_data/processing/alpha_contingency_table.ssv"
print(INPUT_PATH)
print(OUTPUT_PATH)


# 1. Ler arquivo (já ignorando header)
lines <- readLines(INPUT_PATH)[-1]

# 2. Separar por TAB
split_data <- strsplit(lines, "\t")

# 3. Nome das linhas
row_names <- sapply(split_data, `[`, 1)

# 4. Separar as colunas por vírgula
col_list <- lapply(split_data, function(x) {
  as.integer(strsplit(x[2], ",")[[1]])
})

# 5. Construir índices
i <- rep(seq_along(col_list), lengths(col_list))
j <- unlist(col_list)
j <- j + 1
# 6. Criar matriz esparsa
sparse_mat <- sparseMatrix(i = i, j = j, x = 1)

# 7. Nomear linhas
rownames(sparse_mat) <- row_names

# 8. Cutoff
occurrence_mat <- sparse_mat[, colSums(sparse_mat) > 20]
occurrence_mat <- occurrence_mat > 1

OUTPUT_PATH <- "/home/ttamaki/Conjugative-Modules/l_data/processing/alpha_contingency_table.ssv"

p <- ncol(occurrence_mat)
N <- nrow(occurrence_mat)

m <- colSums(occurrence_mat)

cat("g1 g2 a b c d", "\n", file = OUTPUT_PATH, append = TRUE)

for (A in 1:(p-1)){
  print(A)
  for (B in (A+1):p){

    m_A <- m[A]
    m_B <- m[B]

    if(m_A<5 | m_B <5) next

    a <- t(occurrence_mat[,A]) %*% occurrence_mat[,B]
    
    if(a==0) next

    b <- m_B-a
    c <- m_A-a
    d <- N-m_B-c

    data <- c(A,B,a,b,c,d)
    cat(data, "\n", file = OUTPUT_PATH, append = TRUE)
    }
}



