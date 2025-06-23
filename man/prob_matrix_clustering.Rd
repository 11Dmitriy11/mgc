\name{prob_matrix_clustering}
\alias{prob_matrix_clustering}
\title{Probabilistic Matrix Clustering}
\description{
Selects a control cohort using Mahalanobis distance in SVD space.
}
\usage{
prob_matrix_clustering(X, test_idx, n_ctrl = NULL, apriori = NULL,
                       rank_output = FALSE)
}
\arguments{
  \item{X}{Numeric matrix of observations \(\times\) features.}
  \item{test_idx}{Row indices of the test (case) group.}
  \item{n_ctrl}{Number of controls to sample (optional).}
  \item{apriori}{Optional prior weights for each row.}
  \item{rank_output}{Logical; return ranked data.frame if TRUE.}
}
\value{Vector of indices or data.frame with distances/probabilities.}
\examples{
set.seed(1)
X <- matrix(rnorm(1000), ncol = 10)
prob_matrix_clustering(X, 1:10, n_ctrl = 20)
}
