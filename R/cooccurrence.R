#' Extract Pairwise Co-occurrence Frequencies
#'
#' Computes the frequency at which each pair of taxa co-occur across samples
#' in a binary occurrence matrix. This is a core function for differential
#' co-occurrence analysis.
#'
#' @param Xocc Binary occurrence matrix (samples x taxa). Columns represent taxa,
#'  rows represent samples. Matrix elements are 1 for presence and 0 for absence.
#'
#' @return A list containing:
#'   \item{O}{Numeric vector of co-occurrence frequencies for all taxa pairs}
#'   \item{I}{List of index pairs corresponding to each co-occurrence frequency}
#'
#' @details
#' For each pair of taxa (i,j), the co-occurrence frequency is calculated as:
#' \deqn{O_{ij} = \frac{\sum_{k=1}^{n} X_{ki} \cdot X_{kj}}{n}}
#' where n is the number of samples, and X_ki indicates presence/absence of taxon i in sample k.
#'
#' @examples
#' # Create example binary occurrence matrix
#' set.seed(123)
#' occurrence_matrix <- matrix(rbinom(40, 1, 0.5), nrow = 8, ncol = 5)
#' rownames(occurrence_matrix) <- paste0("Sample_", 1:8)
#' colnames(occurrence_matrix) <- paste0("Taxa_", 1:5)
#' 
#' # Extract co-occurrence frequencies
#' result <- extract_cooccurrence_list(occurrence_matrix)
#' print(head(result$O))  # Co-occurrence frequencies
#' print(head(result$I))  # Taxa pair indices
#'
#' @importFrom utils head
#' @export
extract_cooccurrence_list <- function(Xocc) {
  # Validate input
  if (!is.matrix(Xocc)) {
    stop("Input must be a matrix")
  }
  if (any(is.na(Xocc))) {
    warning("Matrix contains NA values which will be treated as 0")
    Xocc[is.na(Xocc)] <- 0
  }
  if (!all(Xocc %in% c(0, 1))) {
    warning("Matrix should contain only 0s and 1s for binary occurrence data")
  }
  
  number_I <- (ncol(Xocc) * (ncol(Xocc) - 1)) / 2
  O <- rep(0, number_I)
  links_I <- vector(mode = 'list', length = number_I)
  idx <- 0
  
  for (i in 1:(ncol(Xocc) - 1)) {
    for (j in (i + 1):ncol(Xocc)) {
      idx <- idx + 1
      # Fraction of samples where taxa pair (i,j) co-occur
      O[idx] <- sum(Xocc[, i] * Xocc[, j], na.rm = TRUE) / nrow(Xocc)
      links_I[[idx]] <- c(i, j)
    }
  }
  
  return(list(O = O, I = links_I))
}

#' Create Occurrence Matrix from Abundance Data
#'
#' Converts a relative abundance matrix to a binary occurrence matrix using
#' a specified abundance threshold.
#'
#' @param abundance_matrix Numeric matrix of relative abundances (taxa x samples)
#' @param detection_threshold Numeric threshold (abundance percentage) for determining taxa presence/absence (default: 0.01%)
#'
#' @return Binary occurrence matrix (samples x taxa)
#'
#' @examples
#' # Create example abundance data
#' set.seed(456)
#' abundance_data <- matrix(runif(40, 0, 0.01), nrow = 5, ncol = 8)
#' rownames(abundance_data) <- paste0("Taxa_", 1:5)
#' colnames(abundance_data) <- paste0("Sample_", 1:8)
#' 
#' # Convert to occurrence matrix
#' occurrence_table <- calculate_cooccurrence_matrix(abundance_data, detection_threshold = 0.1)
#'
#' @export
calculate_cooccurrence_matrix <- function(abundance_matrix, detection_threshold = 0.01) {
  # Validate inputs
  if (!is.matrix(abundance_matrix)) {
    stop("abundance_matrix must be a matrix")
  }
  if (detection_threshold < 0 || detection_threshold > 100) {
    stop("threshold must be between 0 and 100%")
  }

  # Convert to percentage
  perc_table <- abundance_matrix * 100
  
  # Create binary occurrence matrix
  occurrence_matrix <- t((perc_table >= (detection_threshold)) * 1)
  
  return(occurrence_matrix)
}