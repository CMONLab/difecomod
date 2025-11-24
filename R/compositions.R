#' Functions to Handle Compositional Data
#'
#' Internal helper functions used throughout the package.
#'
#' @name compositions
#' @keywords internal
NULL

#' Extract Core Microbiota from Abundance Data 
#'
#' Filters taxa from a relative abundance matrix based on prevalence across samples
#'
#' @param abundance_matrix Numeric matrix of relative abundances (taxa x samples)
#' @param prevalence_threshold Numeric threshold (abundance percentage) for determining taxa prevalence (default: 0.1%)
#' @param min_prevalence Minimum prevalence (fraction of samples) required to retain taxa (default: 0.1)
#'
#' @return Numeric matrix of relative abundances (taxa x samples) with filtered taxa
#'
#' @examples
#' # Create example count abundance data
#' set.seed(456)
#' abundance_matrix <- matrix(runif(40, 0, 0.1), nrow = 5, ncol = 8)
#' colnames(abundance_matrix) <- paste0("Sample_", 1:8)
#' rownames(abundance_matrix) <- paste0("Taxa_", 1:5)
#' 
#' # Convert to occurrence matrix
#' core_abundance_matrix <- core_extraction(abundance_matrix, prevalence_threshold = 5, min_prevalence = 0.5)
#'
#' @export
core_extraction <- function(abundance_matrix, prevalence_threshold = 1, min_prevalence = 0.1) {
  # Validate inputs
  if (!is.matrix(abundance_matrix)) {
    stop("abundance_matrix must be a matrix")
  }
  if (prevalence_threshold < 0 || prevalence_threshold > 100) {
    stop("threshold must be between 0 and 100 %")
  }
  if (min_prevalence < 0 || min_prevalence > 1) {
    stop("min_prevalence must be between 0 and 1")
  }
  
  # Convert to percentage
  perc_table <- abundance_matrix * 100
  
  # Filter taxa based on prevalence
  n_samples <- ncol(abundance_matrix)
  prevalence <- rowSums(perc_table >= prevalence_threshold)
  selected_taxa <- prevalence >= round(min_prevalence * n_samples)
  
  if (sum(selected_taxa) == 0) {
    stop("No taxa meet the prevalence criteria. Consider lowering min_prevalence.")
  }
  
  # Filter abundance matrix with selected taxa
  core_abundance_matrix <- abundance_matrix[selected_taxa, ]
  
  message(sprintf("Filtered to %d taxa (from %d) based on prevalence >= %0.1f%%", 
                  sum(selected_taxa), nrow(abundance_matrix), min_prevalence * 100))
  
  return(core_abundance_matrix)
}

#' Centered Log-Ratio (CLR) transformation
#'
#'transform compositional data in count table by (i) imputing zeros and (ii) taking the logarithm 
#'of each part divided by the geometric mean of all parts in the composition
#'
#' @param count_matrix Numeric matrix of counts (taxa x samples)
#' @param zero_imputation_method Zero count imputation criterion (see cmultRepl function): "GBM" (default), "CZM", or "pseudocount"
#'
#' @return Numeric matrix of clr-transformed abundances (taxa x samples)
#'
#' @importFrom zCompositions cmultRepl
#'
#' @examples
#' # Create example count data
#' set.seed(456)
#' count_data <- matrix(floor(runif(40, 0, 10)), nrow = 5, ncol = 8)
#' colnames(count_data) <- paste0("Sample_", 1:8)
#' rownames(count_data) <- paste0("Taxa_", 1:5)
#'  
#' # Convert to occurrence matrix
#' clr_table <- clr_transformation(count_data, zero_imputation_method="GBM")
#'
#' @export
clr_transformation <- function(count_matrix, zero_imputation_method="GBM") {
  
  # Imputation of zero counts
  if (zero_imputation_method == "pseudocount"){
    #count_matrix[count_matrix == 0] <- 1 # add pseudocount only to null values
    count_matrix <- count_matrix + 1
    clr_table <- count_matrix
  }else{
    clr_table <- zCompositions::cmultRepl(count_matrix, method = zero_imputation_method);
  }
  clr_table <- apply(log(clr_table), 2, function(x){x - mean(x)});
  return(clr_table)
}
