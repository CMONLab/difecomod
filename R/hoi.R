#' Higher-Order Interaction Analysis
#'
#' Functions for analyzing three-way and higher-order microbial interactions.
#'
#' @name hoi
NULL

#' Extract Higher-Order Interaction (HOI) List
#'
#' Computes all possible three-way interactions (triplets) from an occurrence matrix
#' and calculates their co-occurrence frequencies.
#'
#' @param Xocc Binary occurrence matrix (samples x taxa). Rows represent samples,
#'   columns represent features (e.g., microbial taxa).
#'
#' @return A list containing:
#'   \item{O}{Numeric vector of co-occurrence frequencies for each triplet}
#'   \item{I}{List of integer vectors, each containing indices of three features forming a triplet}
#'
#' @details
#' The function calculates the fraction of samples where three features co-occur.
#' For n features, this generates n*(n-1)*(n-2)/6 unique triplets.
#'
#' @examples
#' \dontrun{
#' # Create binary occurrence matrix
#' set.seed(123)
#' occurrence_matrix <- matrix(rbinom(40, 1, 0.5), nrow = 8, ncol = 5)
#' rownames(occurrence_matrix) <- paste0("Sample_", 1:8)
#' colnames(occurrence_matrix) <- paste0("Taxa_", 1:5)
#' hoi_result <- extract_hoi_list(occurrence_matrix)
#' print(head(hoi_result$O))
#' }
#'
#' @export
extract_hoi_list <- function(Xocc) {
  # Validate input
  if (!is.matrix(Xocc)) {
    stop("Xocc must be a matrix")
  }
  if (ncol(Xocc) < 3) {
    stop("Xocc must have at least 3 columns to compute triplets")
  }
  
  # Calculate number of possible triplets
  number_hoi <- (ncol(Xocc) - 1) * (ncol(Xocc) - 2) * (ncol(Xocc) - 3) / 6
  
  # Initialize output vectors
  O <- rep(0, number_hoi)
  links_hoi <- vector(mode = 'list', length = number_hoi)
  
  idx <- 0
  
  # Iterate through all possible triplets
  for (i in 1:(ncol(Xocc) - 2)) {
    for (j in (i + 1):(ncol(Xocc) - 1)) {
      for (k in (j + 1):ncol(Xocc)) {
        idx <- idx + 1
        
        # Calculate fraction of samples where all three features co-occur
        O[idx] <- sum(Xocc[, i] * Xocc[, j] * Xocc[, k], na.rm = TRUE) / nrow(Xocc)
        
        # Store triplet indices
        links_hoi[[idx]] <- c(i, j, k)
      }
    }
  }
  
  return(list(O = O, I = links_hoi))
}


#' Extract Higher-Order Interaction (HOI) Sublist
#'
#' Computes co-occurrence frequencies for a specified subset of triplets (e.g., triangles
#' in a network).
#'
#' @param Xocc Binary occurrence matrix (samples x taxa).
#' @param sub_list Matrix or data frame where each row contains three feature indices
#'   representing a triplet to analyze.
#'
#' @return A list containing:
#'   \item{O}{Numeric vector of co-occurrence frequencies for specified triplets}
#'   \item{I}{List of remapped triplet indices relative to unique features in sub_list}
#'
#' @details
#' This function is useful when analyzing specific network motifs (e.g., triangular structures)
#' rather than all possible triplets.
#'
#' @examples
#' \dontrun{
#' # Create binary occurrence matrix
#' set.seed(123)
#' occurrence_matrix <- matrix(rbinom(40, 1, 0.5), nrow = 8, ncol = 5)
#' rownames(occurrence_matrix) <- paste0("Sample_", 1:8)
#' colnames(occurrence_matrix) <- paste0("Taxa_", 1:5)
#' # Define specific triplets to analyze
#' triangles <- matrix(c(1,2,3, 2,3,4, 1,3,5), nrow = 3, byrow = TRUE)
#' hoi_sub <- extract_hoi_sublist(occurrence_matrix, sub_list = triangles)
#' print(head(hoi_sub$O))
#' }
#'
#' @export
extract_hoi_sublist <- function(Xocc, sub_list = NULL) {
  # Validate inputs
  if (!is.matrix(Xocc) && !is.data.frame(Xocc)) {
    stop("Xocc must be a matrix or data frame")
  }
  if (is.null(sub_list)) {
    stop("sub_list must be provided")
  }
  if (!is.matrix(sub_list) && !is.data.frame(sub_list)) {
    sub_list <- as.matrix(sub_list)
  }
  if (ncol(sub_list) != 3) {
    stop("sub_list must have exactly 3 columns")
  }
  
  # Convert to matrix if needed
  if (is.data.frame(Xocc)) {
    Xocc <- as.matrix(Xocc)
  }
  
  number_hoi <- nrow(sub_list)
  O <- rep(0, number_hoi)
  links_hoi <- vector(mode = 'list', length = number_hoi)
  
  # Get unique indices from triangles list
  new_ind_hoi <- sort(unique(as.vector(sub_list)))
  
  # Calculate co-occurrence for each triplet
  for (i in 1:number_hoi) {
    O[i] <- sum(
      Xocc[, sub_list[i, 1]] * 
        Xocc[, sub_list[i, 2]] * 
        Xocc[, sub_list[i, 3]], 
      na.rm = TRUE
    ) / nrow(Xocc)
    
    # Remap indices to positions in the unique feature set
    links_hoi[[i]] <- c(
      which(new_ind_hoi == sub_list[i, 1]),
      which(new_ind_hoi == sub_list[i, 2]),
      which(new_ind_hoi == sub_list[i, 3])
    )
  }
  
  return(list(O = O, I = links_hoi))
}


#' Random Model for Differential Higher-Order Interactions
#'
#' Generates null distribution for differential HOI analysis using random shuffling.
#'
#' @param X1occ Binary occurrence matrix for condition 1
#' @param X2occ Binary occurrence matrix for condition 2  
#' @param n Number of random column shuffling (default: 5000)
#' @param n_sigma Number of standard deviations from the mean used as threshold (default: 1)
#' @param triangle_list Optional matrix of network triplets on which to test HOI. 
#'   If NULL, all possible triplets are tested.
#' @param values Logical, whether to return full null distribution values (default: FALSE)
#'
#' @return List containing:
#'   \item{Stats}{Data frame with Mean, Median, STD, Min, Max, low_thr, up_thr}
#'   \item{Values}{Matrix of null distribution values (only if values = TRUE)}
#'
#' @details
#' This model preserves occurrence patterns within each condition but breaks
#' feature associations by independently shuffling columns.
#'
#' @examples
#' set.seed(123)
#' X1 <- matrix(rbinom(30, 1, 0.7), nrow = 6, ncol = 5)
#' X2 <- matrix(rbinom(35, 1, 0.5), nrow = 7, ncol = 5)
#' 
#' # Test all triplets (use smaller n for example speed)
#' null_model <- delta_hoi_random_model(X1, X2, n = 100)
#' 
#' # Test specific triplets
#' triangles <- matrix(c(1,2,3, 2,3,4), nrow = 2, byrow = TRUE)
#' null_filtered <- delta_hoi_random_model(X1, X2, n = 100, triangle_list = triangles)
#'
#' @importFrom stats quantile sd median
#' @export
delta_hoi_random_model <- function(X1occ, X2occ, n = 5000, n_sigma = 1, 
                                   triangle_list = NULL, values = FALSE) {
  
  # Validate inputs
  if (!is.matrix(X1occ) || !is.matrix(X2occ)) {
    stop("Inputs must be matrices")
  }
  if (ncol(X1occ) != ncol(X2occ)) {
    stop("Both matrices must have the same number of taxa (columns)")
  }
  if (ncol(X1occ) < 3) {
    stop("Matrices must have at least 3 columns to compute triplets")
  }
  
  null_model_delta_cooccurrence <- NULL
  
  for (i in 1:n) {
    if (!is.null(triangle_list)) {
      null_O1 <- extract_hoi_sublist(
        reshuffle_matrix_columns(X1occ), 
        sub_list = triangle_list
      )$O
      null_O2 <- extract_hoi_sublist(
        reshuffle_matrix_columns(X2occ), 
        sub_list = triangle_list
      )$O
    } else {
      null_O1 <- extract_hoi_list(reshuffle_matrix_columns(X1occ))$O
      null_O2 <- extract_hoi_list(reshuffle_matrix_columns(X2occ))$O
    }
    
    null_delta <- null_O1 - null_O2
    null_model_delta_cooccurrence <- cbind(null_model_delta_cooccurrence, null_delta)
  }
  
  # Calculate summary statistics
  df_null_model_delta_cooccurrence_stats <- data.frame(
    Mean = apply(null_model_delta_cooccurrence, 1, mean),
    Median = apply(null_model_delta_cooccurrence, 1, median),
    STD = apply(null_model_delta_cooccurrence, 1, sd),
    Min = apply(null_model_delta_cooccurrence, 1, min),
    Max = apply(null_model_delta_cooccurrence, 1, max),
    low_thr = apply(null_model_delta_cooccurrence, 1, function(x) {
      mean(x) - n_sigma * sd(x)
    }),
    up_thr = apply(null_model_delta_cooccurrence, 1, function(x) {
      mean(x) + n_sigma * sd(x)
    })
  )
  
  if (values == TRUE) {
    return(list(
      Values = null_model_delta_cooccurrence,
      Stats = df_null_model_delta_cooccurrence_stats
    ))
  } else {
    return(list(Stats = df_null_model_delta_cooccurrence_stats))
  }
}


#' Permutation Model for Differential Higher-Order Interactions
#'
#' Generates null distribution using sample label permutation for HOI analysis.
#'
#' @param X1occ Binary occurrence matrix for condition 1
#' @param X2occ Binary occurrence matrix for condition 2
#' @param n Number of group permutations (default: 5000)
#' @param conf Confidence level for the lower/upper bounds (default: 0.05, 
#'   equivalent to 5th and 95th percentile)
#' @param triangle_list Optional matrix of network triplets on which to test HOI.
#'   If NULL, all possible triplets are tested.
#' @param values Logical, whether to return full null distribution values (default: FALSE)
#'
#' @return List containing:
#'   \item{Stats}{Data frame with Mean, Median, STD, Min, Max, low_thr, up_thr}
#'   \item{Values}{Matrix of null distribution values (only if values = TRUE)}
#'
#' @details
#' This model tests whether differential HOI patterns exceed chance expectations
#' under the null hypothesis that condition labels are arbitrary.
#'
#' @examples
#' set.seed(123)
#' X1 <- matrix(rbinom(30, 1, 0.7), nrow = 6, ncol = 5)
#' X2 <- matrix(rbinom(35, 1, 0.5), nrow = 7, ncol = 5)
#' 
#' # Test all triplets (use smaller n for example speed)
#' perm_model <- delta_hoi_permutation_model(X1, X2, n = 100)
#'
#' @importFrom stats quantile sd median
#' @export
delta_hoi_permutation_model <- function(X1occ, X2occ, n = 5000, conf = 0.05, 
                                        triangle_list = NULL, values = FALSE) {
  
  # Validate inputs
  if (!is.matrix(X1occ) || !is.matrix(X2occ)) {
    stop("Inputs must be matrices")
  }
  if (ncol(X1occ) != ncol(X2occ)) {
    stop("Both matrices must have the same number of taxa (columns)")
  }
  if (ncol(X1occ) < 3) {
    stop("Matrices must have at least 3 columns to compute triplets")
  }
  
  null_model_delta_cooccurrence <- NULL
  n1 <- nrow(X1occ)
  n2 <- nrow(X2occ)
  Xocc <- rbind(X1occ, X2occ)
  
  for (i in 1:n) {
    # Randomly permute sample labels
    Xocc <- Xocc[sample(1:(n1 + n2)), ]
    
    if (!is.null(triangle_list)) {
      null_O1 <- extract_hoi_sublist(Xocc[1:n1, ], sub_list = triangle_list)$O
      null_O2 <- extract_hoi_sublist(Xocc[(n1 + 1):(n1 + n2), ], sub_list = triangle_list)$O
    } else {
      null_O1 <- extract_hoi_list(Xocc[1:n1, ])$O
      null_O2 <- extract_hoi_list(Xocc[(n1 + 1):(n1 + n2), ])$O
    }
    
    null_delta <- null_O1 - null_O2
    null_model_delta_cooccurrence <- cbind(null_model_delta_cooccurrence, null_delta)
  }
  
  # Calculate summary statistics with quantile-based thresholds
  df_null_model_delta_cooccurrence_stats <- data.frame(
    Mean = apply(null_model_delta_cooccurrence, 1, mean),
    Median = apply(null_model_delta_cooccurrence, 1, median),
    STD = apply(null_model_delta_cooccurrence, 1, sd),
    Min = apply(null_model_delta_cooccurrence, 1, min),
    Max = apply(null_model_delta_cooccurrence, 1, max),
    low_thr = apply(null_model_delta_cooccurrence, 1, function(x) {
      quantile(x, c(conf, 1 - conf))[1]
    }),
    up_thr = apply(null_model_delta_cooccurrence, 1, function(x) {
      quantile(x, c(conf, 1 - conf))[2]
    })
  )
  
  if (values == TRUE) {
    return(list(
      Values = null_model_delta_cooccurrence,
      Stats = df_null_model_delta_cooccurrence_stats
    ))
  } else {
    return(list(Stats = df_null_model_delta_cooccurrence_stats))
  }
}
