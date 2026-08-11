library(CooccurrenceAffinity)
library(future.apply)

# =========================================================
# Paralelização
# =========================================================

n_workers <- 8
plan(multicore, workers = n_workers)

# =========================================================
# Argumentos
# =========================================================

args <- commandArgs(trailingOnly = TRUE)

contingency_table_path <- args[1]
output_path <- args[2]

# =========================================================
# Arquivos
# =========================================================

contingency_table <- file(contingency_table_path, "r")
output <- file(output_path, "w")

# =========================================================
# Configurações
# =========================================================

chunk_size <- 1000000

writeLines(
  "p1 p2 a b c d alpha pval LowerLimitCI UpperLimitCI",
  output
)

# =========================================================
# Pula header
# =========================================================

readLines(contingency_table, n = 1)

# =========================================================
# Progresso
# =========================================================

chunk_n <- 0
total_lines <- 0

# =========================================================
# Loop principal
# =========================================================

repeat {

  lines <- readLines(contingency_table, n = chunk_size)

  if (length(lines) == 0) {
    break
  }

  chunk_n <- chunk_n + 1
  total_lines <- total_lines + length(lines)

  cat(
    sprintf(
      "\rChunk %d | linhas processadas: %d",
      chunk_n,
      total_lines
    )
  )

  # =======================================================
  # Divide linhas entre workers
  # =======================================================

  split_lines <- split(
    lines,
    cut(
      seq_along(lines),
      n_workers,
      labels = FALSE
    )
  )

  # =======================================================
  # Processa em paralelo
  # =======================================================

  results <- future_lapply(split_lines, function(chunk) {

    out <- character(length(chunk))

    for (i in seq_along(chunk)) {

      vals <- scan(
        text = chunk[[i]],
        what = numeric(),
        quiet = TRUE
      )

      p1 <- vals[1]
      p2 <- vals[2]

      a  <- vals[3]
      b  <- vals[4]
      c_ <- vals[5]
      d  <- vals[6]

      alpha <- ML.Alpha(
        a,
        c(a + c_, a + b, a + b + c_ + d)
      )

      out[i] <- paste(
        c(
          p1,
          p2,
          a,
          b,
          c_,
          d,
          alpha$est,
          alpha$pval,
          alpha$MedianIntrvl
        ),
        collapse = " "
      )
    }

    out
  })

  # =======================================================
  # Salva
  # =======================================================

  writeLines(
    unlist(results, use.names = FALSE),
    output
  )
}

cat("\nFinalizado.\n")

# =========================================================
# Fecha arquivos
# =========================================================

close(contingency_table)
close(output)