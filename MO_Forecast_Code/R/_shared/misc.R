# ==============================================================================
# File: misc.R
# ==============================================================================
# Purpose
#   Small utility functions used across multiple pipeline stages.
#
# Function index
#   `%||%`(x, y)
#     Null-coalescing operator. Returns x if not NULL, else y.
#
#   softmax_row(x)
#     Numerically stable softmax for a single numeric vector.
#     x: numeric vector of logits. Returns: numeric vector of probabilities.
#
#   row_logsumexp(mat_NK)
#     Row-wise log-sum-exp for an N x K matrix.
#     mat_NK: numeric matrix. Returns: numeric vector of length N.
#
#   inv_logit(x)
#     Numerically stable inverse logit (sigmoid).
#     x: numeric vector. Returns: numeric vector of probabilities.
# ==============================================================================

`%||%` <- function(x, y) if (!is.null(x)) x else y


softmax_row <- function(x){
  z <- x - max(x, na.rm = TRUE)
  ez <- exp(z)
  ez / sum(ez)
}

row_logsumexp <- function(mat_NK) {
  m <- apply(mat_NK, 1, max)
  m + log(rowSums(exp(mat_NK - m)))
}

inv_logit <- function(x) {
  ifelse(x > 0, 1/(1 + exp(-x)), exp(x)/(exp(x) + 1))
}
