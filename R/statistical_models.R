#' Differential Co-occurrence Statistical Models
#'
#' Functions for generating null distributions and statistical testing
#' in differential co-occurrence analysis.
#'
#' @name statistical_models
NULL

#' Random Model for Differential Co-occurrence Analysis
#'
#' Generates null distribution by randomly shuffling occurrence patterns
#' within each taxon while preserving overall occurrence frequencies.
#'
#' @param X1occ Binary occurrence matrix for condition 1 (samples x taxa)
#' @param X2occ Binary occurrence matrix for condition 2 (samples x taxa)
#' @param n Number of random column shuffling (default: 5000)
#'
#' @return List containing:
#'   \item{Values}{Matrix of null co-occurrence models (n_shufflings x n_cooccurrences)}
#'   \item{Stats}{Data frame with mean, median, SD, and 5\% and 95\% percentile values of null co-occurrence models}
#'
#' @details
#' The random model tests whether observed differential co-occurrence patterns
#' could arise by chance given the individual occurrence frequencies of each taxon.
#' This is achieved by randomly shuffling the occurrence patterns within each
#' taxon column, breaking any co-occurrence structure while preserving marginal
#' occurrence rates.
#'
#' @examples
#' # Create example data
#' set.seed(42)
#' X1 <- matrix(rbinom(40, 1, 0.6), nrow = 8, ncol = 5)
#' X2 <- matrix(rbinom(35, 1, 0.4), nrow = 7, ncol = 5)
#' colnames(X1) <- colnames(X2) <- paste0("Taxa_", 1:5)
#' 
#' # Generate null model (use small n for example)
#' null_random <- delta_cooccurrence_random_model(X1, X2, n = 100)
#' head(null_random$Stats)
#'
#' @importFrom stats quantile sd median
#' @export
delta_cooccurrence_random_model <- function(X1occ, X2occ, n = 5000) {
  # Validate inputs
  if (!is.matrix(X1occ) || !is.matrix(X2occ)) {
    stop("Both X1occ and X2occ must be matrices")
  }
  if (ncol(X1occ) != ncol(X2occ)) {
    stop("Both matrices must have the same number of taxa (columns)")
  }
  if (n < 100) {
    warning("Small number of permutations (n < 100) may lead to unstable results")
  }
  
  null_model_delta_cooccurrence <- c()
  
  for (i in 1:n) {
    null_O1 <- extract_cooccurrence_list(reshuffle_matrix_columns(X1occ))$O
    null_O2 <- extract_cooccurrence_list(reshuffle_matrix_columns(X2occ))$O
    null_delta <- null_O1 - null_O2
    null_model_delta_cooccurrence <- cbind(null_model_delta_cooccurrence, null_delta)
  }
  
  df_null_model_delta_cooccurrence_stats <- data.frame(
    Mean = apply(null_model_delta_cooccurrence, 1, mean),
    Median = apply(null_model_delta_cooccurrence, 1, median),
    STD = apply(null_model_delta_cooccurrence, 1, sd),
    Min = apply(null_model_delta_cooccurrence, 1, min),
    Max = apply(null_model_delta_cooccurrence, 1, max),
    C95th_perc = apply(null_model_delta_cooccurrence, 1, function(x) {quantile(x, c(0.05, 0.95))})[1,],
    C5th_perc = apply(null_model_delta_cooccurrence, 1, function(x) {quantile(x, c(0.05, 0.95))})[2,]
  )
  
  return(list(Values = null_model_delta_cooccurrence,
              Stats = df_null_model_delta_cooccurrence_stats))
}

#' Permutation Model for Differential Co-occurrence Analysis
#'
#' Generates null distribution by randomly permuting sample condition labels
#' while preserving group sizes and overall sample structure.
#'
#' @param X1occ Binary occurrence matrix for condition 1
#' @param X2occ Binary occurrence matrix for condition 2
#' @param n Number of group permutations (default: 5000)
#'
#' @return List containing:
#'   \item{Values}{Matrix of null co-occurrence models (n_permutations x n_cooccurrences)}
#'   \item{Stats}{Data frame with mean, median, SD, and 5\% and 95\% percentile values of null co-occurrence models}
#'
#' @details
#' The permutation model tests whether differential co-occurrence patterns
#' could arise simply from random assignment of samples to conditions,
#' preserving the actual group sizes but randomizing which samples belong
#' to which condition.
#'
#' @examples
#' set.seed(42)
#' X1 <- matrix(rbinom(40, 1, 0.6), nrow = 8, ncol = 5)
#' X2 <- matrix(rbinom(35, 1, 0.4), nrow = 7, ncol = 5)
#' 
#' null_perm <- delta_cooccurrence_permutation_model(X1, X2, n = 100)
#' head(null_perm$Stats)
#'
#' @importFrom stats quantile sd median
#' @export
delta_cooccurrence_permutation_model <- function(X1occ, X2occ, n = 5000) {
  # Validate inputs
  if (!is.matrix(X1occ) || !is.matrix(X2occ)) {
    stop("Both X1occ and X2occ must be matrices")
  }
  if (ncol(X1occ) != ncol(X2occ)) {
    stop("Both matrices must have the same number of taxa (columns)")
  }
  
  null_model_delta_cooccurrence <- c()
  n1 <- dim(X1occ)[1]
  n2 <- dim(X2occ)[1]
  Xocc <- rbind(X1occ, X2occ)
  
  for (i in 1:n) {
    Xocc <- Xocc[sample(1:(n1 + n2)), ]
    null_O1 <- extract_cooccurrence_list(Xocc[1:n1, ])$O
    null_O2 <- extract_cooccurrence_list(Xocc[(n1 + 1):(n1 + n2), ])$O
    null_delta <- null_O1 - null_O2
    null_model_delta_cooccurrence <- cbind(null_model_delta_cooccurrence, null_delta)
  }
  
  df_null_model_delta_cooccurrence_stats <- data.frame(
    Mean = apply(null_model_delta_cooccurrence, 1, mean),
    Median = apply(null_model_delta_cooccurrence, 1, median),
    STD = apply(null_model_delta_cooccurrence, 1, sd),
    Min = apply(null_model_delta_cooccurrence, 1, min),
    Max = apply(null_model_delta_cooccurrence, 1, max),
    C95th_perc = apply(null_model_delta_cooccurrence, 1, function(x) {quantile(x, c(0.05, 0.95))})[1,],
    C5th_perc = apply(null_model_delta_cooccurrence, 1, function(x) {quantile(x, c(0.05, 0.95))})[2,]
  )
  
  return(list(Values = null_model_delta_cooccurrence,
              Stats = df_null_model_delta_cooccurrence_stats))
}

#' Validate Binary Occurrence Matrix
#'
#' Performs quality control checks on binary occurrence matrices used in
#' differential co-occurrence analysis.
#'
#' @param occurrence_matrix Binary matrix (samples x taxa)
#' @param min_prevalence Minimum prevalence threshold (default: 0.05)
#' @param max_prevalence Maximum prevalence threshold (default: 0.95)
#' @param min_samples Minimum number of samples required (default: 5)
#'
#' @return List with validation results and filtered matrix
#'
#' @examples
#' set.seed(123)
#' test_matrix <- matrix(rbinom(50, 1, 0.3), nrow = 10, ncol = 5)
#' colnames(test_matrix) <- paste0("Taxa_", 1:5)
#' 
#' validation <- validate_occurrence_matrix(test_matrix)
#' print(validation$summary)
#'
#' @export
validate_occurrence_matrix <- function(occurrence_matrix, min_prevalence = 0.1, 
                                       max_prevalence = 0.9, min_samples = 5) {
  
  # Basic checks
  if (!is.matrix(occurrence_matrix)) {
    stop("Input must be a matrix")
  }
  if (nrow(occurrence_matrix) < min_samples) {
    stop(sprintf("Matrix must have at least %d samples (rows)", min_samples))
  }
  
  # Check for binary data
  unique_vals <- unique(as.vector(occurrence_matrix))
  if (!all(unique_vals %in% c(0, 1, NA))) {
    warning("Matrix contains non-binary values. Converting to binary using threshold 0.5")
    occurrence_matrix <- (occurrence_matrix > 0.5) * 1
  }
  
  # Calculate prevalences
  prevalences <- colMeans(occurrence_matrix, na.rm = TRUE)
  
  # Filter taxa based on prevalence
  valid_taxa <- (prevalences >= min_prevalence) & (prevalences <= max_prevalence)
  
  summary_stats <- list(
    n_samples = nrow(occurrence_matrix),
    n_taxa_original = ncol(occurrence_matrix),
    n_taxa_filtered = sum(valid_taxa),
    prevalence_range = range(prevalences),
    n_missing_values = sum(is.na(occurrence_matrix))
  )
  
  filtered_matrix <- occurrence_matrix[, valid_taxa, drop = FALSE]
  
  return(list(
    filtered_matrix = filtered_matrix,
    valid_taxa = valid_taxa,
    prevalences = prevalences,
    summary = summary_stats
  ))
}

#' Compute Thresholds from Random Model
#'
#' This function calculates upper and lower statistical thresholds for each variable
#' based on the mean and standard deviation of the null distribution from the random model.
#' You can control the width of the threshold with the \code{n_sigma} parameter.
#'
#' @param values A numerical matrix (n_permutations x n_cooccurrences) containing the null co-occurrence distributions 
#'   values (e.g., \code{random_model$Values}).
#' @param n_sigma Number of standard deviations from the mean used as threshold (default: 1).
#' @return A list with two vectors: \code{upper} (upper thresholds) and \code{lower} (lower thresholds), each of length nrow(values).
#' @examples
#' \dontrun{
#' # Assuming you have run delta_cooccurrence_random_model:
#' # random_model <- delta_cooccurrence_random_model(X1, X2, n = 5000)
#' upper_lower <- compute_random_model_thresholds(random_model$Values, n_sigma = 1)
#' upper_thr_rand <- upper_lower$upper
#' lower_thr_rand <- upper_lower$lower
#' }
#' @export
compute_random_model_thresholds <- function(values, n_sigma = 1) {
  upper <- apply(values, 1, function(x) mean(x, na.rm = TRUE) + n_sigma * sd(x, na.rm = TRUE))
  lower <- apply(values, 1, function(x) mean(x, na.rm = TRUE) - n_sigma * sd(x, na.rm = TRUE))
  return(list(upper = upper, lower = lower))
}

#' Compute Thresholds from Permutation Model
#'
#' This function calculates upper and lower statistical thresholds for each variable
#' based on quantiles of the null distribution from the permutation model.
#'
#' @param values A numerical matrix (n_permutations x n_cooccurrences) containing the null co-occurrence distributions
#'   values (e.g., \code{perm_model$Values}).
#' @param conf Confidence level for the lower/upper bounds (default: 0.05, equivalent to 5th and 95th percentile).
#' @return A list with two vectors: \code{upper} (upper thresholds) and \code{lower} (lower thresholds), each of length nrow(values).
#' @examples
#' \dontrun{
#' # Assuming you have run delta_cooccurrence_permutation_model:
#' # perm_model <- delta_cooccurrence_permutation_model(X1, X2, n = 5000)
#' upper_lower <- compute_permutation_model_thresholds(perm_model$Values, conf = 0.05)
#' upper_thr_perm <- upper_lower$upper
#' lower_thr_perm <- upper_lower$lower
#' }
#' @export
compute_permutation_model_thresholds <- function(values, conf = 0.05) {
  qmat <- apply(values, 1, function(x) quantile(x, c(conf, 1 - conf), na.rm = TRUE))
  lower <- qmat[1, ]
  upper <- qmat[2, ]
  return(list(upper = upper, lower = lower))
}