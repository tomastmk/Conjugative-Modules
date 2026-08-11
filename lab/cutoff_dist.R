
getwd()

dados <- read.table("/home/ttamaki/Conj/data/processing/graphs/graph.csv",header=TRUE, sep=",")
head(dados)
jaccard = dados[,c(1,2,7)]
ochiai = dados[,c(1,2,8)]
exp = dados[,c(1,2,9)]


find_cutoff <- function(values, plot = TRUE) {
  # Ordenar valores
  sorted_values <- sort(values)
  len <- length(values)
  
  distances <- numeric(len)
  
  # Reta entre primeiro e último ponto (y = mx + n)
  m <- (sorted_values[len] - sorted_values[1]) / (len - 1)
  n <- sorted_values[1] - m * 1
  
  for (idx in seq_len(len)) {
    y_value <- sorted_values[idx]
    
    # intercepto da perpendicular (y = -1/m * x + b)
    b <- y_value - (-1 / m) * (idx - 1)
    
    # interseção entre reta principal e perpendicular
    intersect_y <- (n + b * m^2) /
      (1 + m^2)
    
    intersect_x <- (intersect_y - n) / m
    
    # distância euclidiana
    distance <- sqrt((intersect_x - (idx - 1))^2 + (intersect_y - y_value)^2)
    
    distances[idx] <- distance
  }
  

  cutoff_index <- which.max(distances)
  return(sorted_values[cutoff_index])
}

jaccard_cutoff <- find_cutoff(jaccard[,"jaccard"],plot=FALSE)
ochiai_cutoff <- find_cutoff(ochiai[,"ochiai"],plot=FALSE)
exp_cutoff <- find_cutoff(exp[,"exp"],plot=FALSE)

jaccard_filtered <- jaccard[jaccard["jaccard"]>jaccard_cutoff,]
ochiai_filtered <- ochiai[ochiai["ochiai"]>ochiai_cutoff,]
exp_filtered <- exp[exp["exp"]>exp_cutoff,]

hist(jaccard[jaccard["jaccard"]>0,]$jaccard)

print(
  sprintf(
    "Em jaccard restaram %.2f%% dos dados",
    nrow(jaccard_filtered) / nrow(jaccard) * 100
  )
)

print(
  sprintf(
    "Em exp restaram %.2f%% dos dados",
    nrow(exp_filtered) / nrow(exp) * 100
  )
)

print(
  sprintf(
    "Em ochiai restaram %.2f%% dos dados",
    nrow(ochiai_filtered) / nrow(ochiai) * 100
  )
)


par(mfrow=c(2,3)) # Divide a tela em 3

boxplot(ochiai["ochiai"], ylim=c(0,1))
boxplot(exp["exp"], ylim=c(0,1))
boxplot(jaccard["jaccard"], ylim=c(0,1))
boxplot(ochiai_filtered["ochiai"], ylim=c(0,1))
boxplot(exp_filtered["exp"], ylim=c(0,1))
boxplot(jaccard_filtered["jaccard"], ylim=c(0,1))

  
library(dplyr)
library(tidyr)

merged <- ochiai_filtered %>%
  full_join(exp_filtered, by = c("gene1", "gene2")) %>%
  full_join(jaccard_filtered, by = c("gene1", "gene2")) %>%
  mutate(across(where(is.numeric), ~replace_na(., 0)))

output <- "/home/ttamaki/Conj/data/graph_filtered.tsv"
write.table(merged,
          file = output,
          sep = "\t",
          row.names = FALSE,
          quote = FALSE) 
