#' Internal utilities
#' @keywords internal
#' @param_doc Common arguments
NULL

#' @title irnt_transform
#' @description Inverse-rank normal transform (Blom).
#' @keywords internal
irnt_transform <- function(x) {
  r <- rank(x, ties.method = "average")
  qnorm((r - 0.5) / length(r))
}
