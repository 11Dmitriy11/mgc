args <- commandArgs(trailingOnly = TRUE)

# Подключение библиотек
library(dplyr)
library(tidyr)
library(igraph)
library(mcmcRanking)
library(BioNet)
library(xgboost)

# Функции
min_max_standardize <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

fisher_combine <- function(pvals) {
  X2 <- -2 * sum(log(pvals))
  df <- 2 * length(pvals)
  1 - pchisq(X2, df)
}

find_k_closest_gene_paths <- function(graph, gene_names, k) {
  k_closest_paths <- list()
  
  for (v in V(graph)) {
    paths <- list()
    path_lengths <- numeric()
    
    for (gene in gene_names) {
      if (gene %in% V(graph)$name) {
        path <- shortest_paths(graph, from = v, to = gene, output = "vpath")$vpath[[1]]
        path_lengths[gene] <- length(path) - 1
        paths[[gene]] <- path
      }
    }
    
    sorted_genes <- names(sort(path_lengths))
    k_shortest_paths <- lapply(sorted_genes[1:min(k, length(sorted_genes))], function(g) paths[[g]])
    k_closest_paths[[V(graph)$name[v]]] <- k_shortest_paths
  }
  
  k_closest_paths
}

# Загрузка данных
edges <- read.table(args[1], sep = '\t', header = FALSE, col.names = c('u', 'v'))
g <- graph_from_edgelist(as.matrix(edges))
g <- igraph::simplify(g)

true_genes <- read.table(args[2], sep = '\t', header = FALSE)$V1

# Проверка компоненты связности
components <- components(g)
largest_comp_id <- which.max(components$csize)
general_subgraph <- induced_subgraph(g, which(components$membership == largest_comp_id))

# Параметры
l <- 3  # задаем количество раундов, можно вынести как аргумент
set.seed(123)

# Случайные p-value для демонстрации
random <- runif(vcount(general_subgraph))
names(random) <- V(general_subgraph)$name

# Обработка
shortest_paths_to_genes <- find_k_closest_gene_paths(general_subgraph, true_genes, 3)

# Сбор путей (по количеству вершин)
get_length <- function(path) length(path)
first_elements <- sapply(shortest_paths_to_genes, function(x) if (length(x) >= 1) get_length(x[[1]]) else NA)
second_elements <- sapply(shortest_paths_to_genes, function(x) if (length(x) >= 2) get_length(x[[2]]) else NA)
third_elements <- sapply(shortest_paths_to_genes, function(x) if (length(x) >= 3) get_length(x[[3]]) else NA)

# p-value и объединение
d <- data.frame(name = V(general_subgraph)$name)
for (i in 1:l) {
  V(general_subgraph)$pval <- random[V(general_subgraph)$name]
  d[, paste0("likelihood", i)] <- V(general_subgraph)$pval
}

if (l == 1) {
  d$Product <- d$likelihood1
} else {
  d$Product <- apply(d[, paste0("likelihood", 1:l)], 1, fisher_combine)
}

V(general_subgraph)$pval <- d$Product
fdr <- quantile(V(general_subgraph)$pval, 0.05)
general_subgraph <- set_likelihood(graph = general_subgraph, fdr = as.numeric(fdr))

# MCMC выбор
z <- mcmc_sample(
  graph = general_subgraph,
  times = 100,
  niter = 10000,
  exp_lh = 1 / 2 ^ (l:0)
)

p <- get_frequency(z, prob = TRUE)
d_f <- data.frame(prob = p, names = names(p))
d_f$y <- ifelse(d_f$names %in% true_genes, 1, 0)
d_f$pval <- V(general_subgraph)$pval[match(d_f$names, V(general_subgraph)$name)]
d_f$first <- first_elements[match(d_f$names, names(first_elements))]
d_f$second <- second_elements[match(d_f$names, names(second_elements))]
d_f$third <- third_elements[match(d_f$names, names(third_elements))]

# Подготовка данных для XGBoost
features <- d_f %>% select(prob, pval, first, second, third) 
label <- d_f$y
dtrain <- xgb.DMatrix(data = as.matrix(features), label = label)

# Обучение XGBoost
model_xgb <- xgboost(data = dtrain, max_depth = 3, eta = 0.1, nrounds = 100, objective = "binary:logistic", verbose = 0)
d_f$pred <- predict(model_xgb, as.matrix(features))

# Запись результатов
write.table(d_f, file = 'test.txt', sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)


