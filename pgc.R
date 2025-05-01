args <- commandArgs(trailingOnly = TRUE)
library(dplyr)
library(readxl)
library(tidyr)
library(mcmcRanking)
library(igraph)
library(BioNet)
library(ROCR)
library(xgboost)

   
min_max_standardize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}  
'%!in%' <- function(x,y)!('%in%'(x,y))

fisher_combine <- function(pvals) {
  X2 <- -2 * sum(log(pvals))
  df <- 2 * length(pvals)
  combined_p <- 1 - pchisq(X2, df)
  return(combined_p)
}

get_random_connected_subgraph <- function(graph, subgraph_size) {
   V(graph)$likelihood <- 1
    z <-
      mcmc_sample(
        graph = graph,
        times = 1,
        niter = 1e4,
        subgraph_order =  subgraph_size
      )
    p <- get_frequency(z)
    selected_vertices <- names(p[p>0])
subgraph <- induced_subgraph(graph, selected_vertices)
return(subgraph)
}

get_connected_subgraph <- function(graph, target_size) {
  
  if (vcount(graph) < target_size) {
    stop("The graph does not have enough vertices.")
  }
  start_vertex <- sample(V(graph), 1)
  
  bfs_result <- igraph::bfs(graph, root = start_vertex, dist = TRUE)
  
  if (length(bfs_result$order) < target_size) {
    stop("The graph is not large enough to contain a subgraph with the target size.")
  }
  
  selected_vertices <- bfs_result$order[1:target_size]
  
  subgraph <- induced_subgraph(graph, selected_vertices)
  
  return(subgraph)
}

find_k_closest_gene_paths <- function(graph, gene_names, k) {
  k_closest_paths <- list()
  
  for (v in V(graph)) {
    paths <- list()
    path_lengths <- numeric()
    
    for (gene in gene_names) {
      if (gene %in% V(graph)$name) {
        path <- shortest_paths(graph, from = v, to = which(V(graph)$name == gene), output = "vpath")$vpath[[1]]
        path_length <- length(path) - 1  
        
        paths[[gene]] <- path
        path_lengths[gene] <- path_length
      }
    }
    
    sorted_genes <- names(sort(path_lengths, decreasing = FALSE))
    k_shortest_paths <- lapply(sorted_genes[1:min(k, length(sorted_genes))], function(gene) paths[[gene]])
    
    k_closest_paths[[V(graph)$name[v]]] <- k_shortest_paths
  }
  
  return(k_closest_paths)
}

extract_4_element <- function(sublist) {
  return(sublist[4])
}
extract_5_element <- function(sublist) {
  return(sublist[5])
}
extract_3_element <- function(sublist) {
  return(sublist[3])
}
extract_2_element <- function(sublist) {
  return(sublist[2])
}
extract_1_element <- function(sublist) {
  return(sublist[1])
}
count_vertices <- function(sublist) {
  return(length(sublist[[1]]))
}

edges <- read.table(args[1], sep='\t',header=F, col.names=c('u', 'v'))
g <- graph_from_edgelist(as.matrix(edges))
g <- largestComp(g)
g <- simplify(g)

true_genes <- read.table(args[2], sep='\t',header=F)
true_genes <- true_genes$V1
shortest_paths_to_genes <- find_k_closest_gene_paths(as.undirected(general_subgraph), true_genes, 3)

first_elements <- sapply(lapply(shortest_paths_to_genes, extract_1_element),count_vertices)
second_elements <- sapply(lapply(shortest_paths_to_genes, extract_2_element),count_vertices)
third_elements <- sapply(lapply(shortest_paths_to_genes, extract_3_element),count_vertices)

for (i in 1:l) {
  V(general_subgraph)$likelihood <- V(general_subgraph)$pval
if (i == 1) {
  V(general_subgraph)$pval <- random[V(general_subgraph)$name]
  d <-  data.frame(name = V(general_subgraph)$name, likelihood = V(general_subgraph)$likelihood)
}
else { 
  d[,paste("likelihood",i,sep="")] <- V(general_subgraph)$likelihood
}
}

if (l == 1) {
  d$Product <- d$likelihood
  V(general_subgraph)$pval <- d$Product
  fdr <- quantile(V(general_subgraph)$pval,0.05)
  general_subgraph <- set_likelihood(graph = general_subgraph, fdr = as.numeric(fdr))
} else {
  d$Product <- apply(d[, 2:(l+1)], 1, fisher_combine)
  V(general_subgraph)$pval <- d$Product
  fdr <- quantile(V(general_subgraph)$pval,0.05)
  general_subgraph <- set_likelihood(graph = general_subgraph, fdr = as.numeric(fdr))
}

V(general_subgraph)$pval <- d$Product
fdr <- quantile(V(general_subgraph)$pval,0.05)
general_subgraph <- set_likelihood(graph = general_subgraph, fdr = as.numeric(fdr))


z <-
  mcmc_sample(
    graph = general_subgraph,
    times = 1e2,
    niter = 1e4,
    exp_lh = 1 / 2 ^ (depth:0)
  )

p <- get_frequency(z, prob = TRUE)
d_f <- data.frame(prob = p, names = names(p))

d_f <- d_f %>% mutate(y = case_when(names %in% names(subgraph_vertices) ~ 1,
                                    TRUE ~ 0))

named_vector <- setNames(V(general_subgraph)$pval, V(general_subgraph)$name)
d_f$pval <-  named_vector[row.names(d_f)]

d_f$first <- first_elements[row.names(d_f)]
d_f$second <- second_elements[row.names(d_f)]
d_f$third <- third_elements[row.names(d_f)]

d_f$y <- ifelse(row.names(d_f) %in% true_genes == 1, 0)
row.names(RES) <- 1:nrow(RES)

label <- RES$y

tr <- xgb.DMatrix(data = as.matrix(RES[,-c(2,3,4,ncol(RES))]), label = as.numeric(RES$y))
model_xgb <- xgboost(data = tr, label = label,max_depth = 3, eta = 0.1, nrounds = 100, objective = "binary:logistic")
res <- predict(model_xgb, as.matrix(RES[,-c(2,3,4,ncol(RES))]),  type = "prob")
RES$pred <- res







