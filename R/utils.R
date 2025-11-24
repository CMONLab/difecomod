#' Utility Functions for Statistical Calculations
#'
#' Internal helper functions used throughout the package.
#'
#' @name utils
#' @keywords internal
NULL

#' Get Mode of Distribution
#'
#' Calculates the mode (most frequent value) of a vector.
#'
#' @param v Numeric or character vector
#' @return The most frequent value in the vector
#' @keywords internal
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#' Reshuffle Matrix Columns
#'
#' Randomly reshuffles the elements within each column of a matrix.
#' This function is used in the random null model generation.
#'
#' @param Xocc Binary occurrence matrix (taxa x samples)
#' @return Matrix with the same dimensions but shuffled column elements
#' @keywords internal
reshuffle_matrix_columns <- function(Xocc) {
  randXocc <- apply(Xocc, 2, sample)
  return(randXocc)
}
