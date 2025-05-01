args = commandArgs(trailingOnly = TRUE)

library(tidyr)
library(dplyr)
library(tools)

CDBL_PCA <- read.table(args[1], sep = '\t', header = TRUE)
ECCE <- read.table(args[2], sep = '\t', header = TRUE)
pheno <- read.table(args[3], sep = '\t', header = FALSE)$V1

CDBL_PCA <- CDBL_PCA %>% dplyr::select(all_of(pheno))
ECCE <- ECCE %>% dplyr::select(all_of(pheno))
print(head(ECCE))
print(nrow(ECCE))
print(head(CDBL_PCA))
print(nrow(CDBL_PCA))

pca_test <- prcomp(CDBL_PCA, center = TRUE, scale. = TRUE)

# Определяем оптимальное количество компонент (например, ≥90% объяснённой дисперсии)
explained_var <- pca_test$sdev^2 / sum(pca_test$sdev^2)
cum_explained <- cumsum(explained_var)
optimal_pc <- which(cum_explained >= 0.4)[1]

# Сообщаем пользователю
cat('Selected', optimal_pc, 'PC, which explain', round(cum_explained[optimal_pc] * 100, 2), '% variance.\n')

# Берём оптимальное количество компонент
test_proj <- pca_test$x[, 1:optimal_pc, drop = FALSE]
print(head(test_proj))
print(nrow(test_proj))
eigenvecs <- pca_test$rotation[, 1:optimal_pc, drop = FALSE]

test_means <- colMeans(test_proj)
test_vars <- apply(test_proj, 2, var)

ECCE_centered <- scale(ECCE, center = pca_test$center, scale = FALSE)
ECCE_proj <- as.matrix(ECCE_centered) %*% eigenvecs

print(head(ECCE_proj))
print(nrow(ECCE_proj))

cov_mat <- diag(test_vars)
d2 <- mahalanobis(ECCE_proj, center = test_means, cov = cov_mat)

threshold <- qchisq(0.95, df = 2)

selected_controls <- ECCE[d2 <= threshold, ]
selected_indices <- which(d2 <= threshold)

test_file_name <- basename(args[1])
test_file_base <- file_path_sans_ext(test_file_name)  
output_file <- file.path(getwd(), paste0(test_file_base, '_selected_indices.txt'))

write.table(selected_indices, file = output_file, sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
