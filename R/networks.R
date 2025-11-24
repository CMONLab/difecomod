#' Network Construction and Analysis
#'
#' Functions for constructing and analyzing differential co-occurrence networks.
#'
#' @name networks
NULL

#' Construct Differential Co-occurrence Network
#'
#' Creates a network representation of significant differential co-occurrences
#' based on random and permutation model thresholds.
#'
#' @param delta_I Numeric vector of observed differential co-occurrence values
#' @param links_I List of taxa pair indices corresponding to delta_I
#' @param nodes_I Character vector of node (taxa) names
#' @param stat_thr_rand list of 2 containing $ lower and $ upper statistical thresholds (from delta_cooccurrence_random_model)
#' @param stat_thr_perm list of 2 containing $ lower and $ upper statistical thresholds (from delta_cooccurrence_permutation_model)
#'
#' @return Adjacency matrix representing the differential co-occurrence network.
#'   Values: +1 = over-represented in condition 1, -1 = under-represented in condition 1,
#'   0 = not significant
#'
#' The sign indicates the direction of the differential co-occurrence.
#'
#' @examples
#' # This example requires results from the statistical models
#' # See package vignette for complete workflow
#' 
#' @export
construct_network_delta <- function(delta_I = NULL, links_I = NULL, nodes_I = NULL, 
                                    stat_thr_rand = NULL, stat_thr_perm = NULL) {
  
  # Validate inputs
  if (any(sapply(list(delta_I, links_I, nodes_I, stat_thr_rand, stat_thr_perm), is.null))) {
    stop("All parameters must be provided")
  }
  
  # Determine significance while keeping permutation model sign
  sign_rand <- (delta_I > stat_thr_rand$upper)* 1 + (delta_I < stat_thr_rand$lower)* 1;
  sign_perm <- (delta_I > stat_thr_perm$upper)* 1 + (delta_I < stat_thr_perm$lower)* -1;
  
  if (sum((sign_perm * sign(delta_I)) < 0)>0){print("warning: permutation direction inconsistent!")}
  Fd <- sign_rand * sign_perm; 
  
  if(sum(Fd)==0){
    print("No significant interactions!")
    return(0)
  }
  Net <- matrix(data=0, nrow=length(nodes_I), ncol=length(nodes_I));
  
  colnames(Net) <- nodes_I;
  rownames(Net) <- nodes_I;
  
  for (l in which(Fd!=0)){
    # keep trace of Net for pheatmap visualization
    Net[links_I[[l]], links_I[[l]]] <- Fd[l];
  }
  Net <- Net[rowSums(abs(Net))!=0, colSums(abs(Net))!=0];
  diag(Net) <- 0;
  
  return(Net)
}


#' Extract Consistent-Sign Triangles from Network
#'
#' Identifies all triangular motifs (3-node cliques) in a network and filters
#' to retain only those where all three edges have the same sign.
#'
#' @param adj_matrix Numeric adjacency matrix representing the network. 
#'   Positive values indicate positive associations, negative values indicate 
#'   negative associations. The matrix should be square and symmetric.
#' @param filter_sign Logical, whether to filter for consistent-sign triangles only.
#'   If TRUE (default), returns only triangles where all edges are positive OR all negative.
#'   If FALSE, returns all triangles regardless of edge signs.
#'
#' @return A list containing:
#'   \item{triangles_list}{Matrix where each row contains three node indices forming a triangle}
#'   \item{triangle_taxa_indices}{Numeric vector of unique taxa indices involved in triangles}
#'   \item{triangle_taxa}{Character vector of unique taxa names involved in triangles}
#'   \item{n_triangles}{Total number of triangles found}
#'   \item{n_positive}{Number of triangles with all positive edges}
#'   \item{n_negative}{Number of triangles with all negative edges}
#'
#' @details
#' This function:
#' 1. Converts the adjacency matrix to an igraph object
#' 2. Identifies all triangular motifs using igraph::triangles()
#' 3. Sorts node indices within each triangle for consistency
#' 4. Optionally filters to keep only triangles with consistent edge signs
#' 5. Returns triangle indices and associated feature names
#'
#' Triangles with consistent signs represent coherent ecological modules:
#' - All positive: mutualistic/cooperative associations
#' - All negative: competitive/antagonistic associations
#'
#' @examples
#' \dontrun{
#' # Create example network matrix
#' set.seed(123)
#' network_mat <- matrix(rnorm(25), 5, 5)
#' network_mat[network_mat > 0] <- 1
#' network_mat[network_mat < 0] <- -1
#' diag(network_mat) <- 0
#' rownames(network_mat) <- colnames(network_mat) <- paste0("Taxa", 1:5)
#' 
#' # Extract consistent-sign triangles
#' triangles_result <- extract_network_triangles(network_mat, filter_sign = TRUE)
#' print(triangles_result$n_triangles)
#' print(triangles_result$triangle_taxa)
#' }
#'
#' @export
extract_network_triangles <- function(adj_matrix, filter_sign = TRUE) {
  
  # Validate input
  if (!is.matrix(adj_matrix)) {
    stop("adj_matrix must be a matrix")
  }
  if (nrow(adj_matrix) != ncol(adj_matrix)) {
    stop("adj_matrix must be square")
  }
  if (nrow(adj_matrix) < 3) {
    stop("Network must have at least 3 nodes to form triangles")
  }
  
  # Create igraph object from absolute values (to find structural triangles)
  net_igraph <- igraph::graph_from_adjacency_matrix(
    abs(adj_matrix), 
    mode = "undirected"
  )
  
  # Extract all triangles from the network
  triangles_raw <- igraph::triangles(net_igraph)
  
  # Check if any triangles exist
  if (length(triangles_raw) == 0) {
    message("No triangles found in the network")
    return(list(
      triangles_list = matrix(nrow = 0, ncol = 3),
      triangle_taxa_indices = 0,
      triangle_taxa = character(0),
      n_triangles = 0,
      n_positive = 0,
      n_negative = 0
    ))
  }
  
  # Convert to matrix format (each row is one triangle)
  triangles_list <- t(matrix(triangles_raw, nrow = 3))
  
  # Sort node indices within each triangle for consistency
  triangles_list <- t(apply(triangles_list, 1, sort))
  
  # Initialize counters
  n_positive <- 0
  n_negative <- 0
  
  # Filter for consistent-sign triangles if requested
  if (filter_sign) {
    # Helper function to get edge signs for a triangle
    get_edge_signs <- function(triangle) {
      i <- triangle[1]; j <- triangle[2]; k <- triangle[3]
      c(adj_matrix[i, j], adj_matrix[j, k], adj_matrix[k, i])
    }
    
    # Check each triangle for sign consistency
    is_consistent <- apply(triangles_list, 1, function(tri) {
      signs <- get_edge_signs(tri)
      all(signs > 0) || all(signs < 0)
    })
    
    # Count positive and negative triangles
    triangle_types <- apply(triangles_list, 1, function(tri) {
      signs <- get_edge_signs(tri)
      if (all(signs > 0)) return("positive")
      if (all(signs < 0)) return("negative")
      return("mixed")
    })
    
    n_positive <- sum(triangle_types == "positive")
    n_negative <- sum(triangle_types == "negative")
    
    # Filter to keep only consistent triangles
    triangles_list <- triangles_list[is_consistent, , drop = FALSE]
    
    message(sprintf(
      "Found %d total triangles: %d positive, %d negative, %d mixed",
      length(is_consistent), n_positive, n_negative, 
      sum(!is_consistent)
    ))
    message(sprintf("Retained %d consistent-sign triangles", nrow(triangles_list)))
  } else {
    message(sprintf("Found %d triangles (no filtering applied)", nrow(triangles_list)))
  }
  
  # Extract unique feature names involved in triangles
  if (nrow(triangles_list) > 0) {
    unique_indices <- sort(unique(as.vector(triangles_list)))
    triangle_taxa <- rownames(adj_matrix)[unique_indices]
  } else {
    unique_indices <- character(0)
    triangle_taxa <- character(0)
  }
  
  # Return results
  return(list(
    triangles_list = triangles_list,
    triangle_taxa_indices = unique_indices,
    triangle_taxa = triangle_taxa,
    n_triangles = nrow(triangles_list),
    n_positive = n_positive,
    n_negative = n_negative
  ))
}


#' Construct Higher-Order Interaction Diagram
#'
#' Creates a matrix representation of significant triplet interactions for 
#' visualization, combining results from both random and permutation model thresholds.
#'
#' @param delta_I Numeric vector of observed differential HOI values (condition 1 - condition 2)
#' @param links_I List of triplet indices, where each element contains three feature indices
#' @param nodes_I Character vector of node (taxa) names involved in the triplets
#' @param stat_thr_rand List of 2 containing $lower and $upper statistical thresholds 
#'   (from delta_hoi_random_model)
#' @param stat_thr_perm List of 2 containing $lower and $upper statistical thresholds 
#'   (from delta_hoi_permutation_model)
#'
#' @return Matrix where rows are features and columns are significant triplets.
#'   Values: +1 = over-represented in condition 1 (healthy), 
#'           -1 = under-represented in condition 1 (enriched in condition 2/disease),
#'           0 = not significant.
#'   Returns 0 if no significant HOI detected.
#'
#' The sign indicates the direction of differential HOI. This function creates an
#' UpSet-style visualization matrix where each column represents one significant triplet
#' and rows show which taxa are involved.
#'
#' @examples
#' # This example requires results from the statistical models
#' # See package vignette for complete workflow
#' 
#' @export
construct_hoi_diagram <- function(delta_I = NULL, links_I = NULL, nodes_I = NULL, 
                                  stat_thr_rand = NULL, stat_thr_perm = NULL) {
  
  # Validate inputs
  if (any(sapply(list(delta_I, links_I, nodes_I, stat_thr_rand, stat_thr_perm), is.null))) {
    stop("All parameters must be provided")
  }
  
  # Determine significance while keeping permutation model sign
  sign_rand <- (delta_I > stat_thr_rand$upper) * 1 + (delta_I < stat_thr_rand$lower) * 1
  sign_perm <- (delta_I > stat_thr_perm$upper) * 1 + (delta_I < stat_thr_perm$lower) * -1
  
  # Check for inconsistencies between permutation direction and observed values
  if (sum((sign_perm * sign(delta_I)) < 0) > 0) {
    print("warning: permutation direction inconsistent!")
  }
  
  # Combine both tests: only significant if both models agree
  Fd <- sign_rand * sign_perm
  
  # Check if any significant HOI exist
  if (sum(Fd) == 0) {
    print("No significant HOI!")
    return(0)
  }
  
  # Create UpSet-style matrix
  Upset <- matrix(data = 0, nrow = length(nodes_I), ncol = sum(Fd != 0))
  colnames(Upset) <- seq(1:ncol(Upset))
  rownames(Upset) <- nodes_I
  
  # # Fill matrix with significant triplet information
  # for (l in which(Fd != 0)) {
  #   Upset[links_I[[l]], which(which(Fd != 0) == l)] <- Fd[l]
  # }
  idx <-0;
  for (l in which(Fd!=0)){
    idx <- idx + 1;
    Upset[links_I[[l]], idx] <- Fd[l];
  }  
  
  return(Upset)
}

#' Validate Higher Order Co-occurrences Against Pairwise Co-occurrences
#'
#' Filters HOI diagram by keeping only those triplets where the sign of the triplet
#' matches all three corresponding pairwise signs in the given network matrix.
#'
#' @param hoi_diagram Matrix where rows are taxa, columns are triplets.
#'   Values are +1 (enriched in condition 1), -1 (enriched in condition 2), or 0.
#'   Row names must be taxa names.
#' @param network_matrix Square adjacency matrix with pairwise differential signs.
#'   Row and column names must be taxa names matching those in hoi_diagram.
#'
#' @return Filtered HOI diagram containing only triplets where all three pairwise
#'   associations have the same sign as the triplet. Returns 0 if none remain.
#'
#' @details
#' For each triplet column in the diagram:
#' 1. Identifies the three taxa with non-zero values
#' 2. Extracts their triplet sign (+1 or -1)
#' 3. Checks all three pairwise signs in network_matrix
#' 4. Keeps triplet only if all pairwise signs match the triplet sign
#'
#' This ensures biological coherence: higher-order patterns align with
#' underlying pairwise associations.
#'
#' @examples
#' \dontrun{
#' hoi_validated <- validate_hoi_diagram(
#'   hoi_diagram = hoi_diagram_raw,
#'   network_matrix = network_matrix
#' )
#' }
#'
#' @export
validate_hoi_diagram <- function(hoi_diagram, network_matrix) {
  
  # Input validation
  if (!is.matrix(hoi_diagram) || ncol(hoi_diagram) == 0) {
    return(0)
  }
  
  if (!is.matrix(network_matrix)) {
    stop("network_matrix must be a matrix")
  }
  
  if (is.null(rownames(hoi_diagram)) || is.null(rownames(network_matrix))) {
    stop("Both hoi_diagram and network_matrix must have row names (taxa names)")
  }
  
  coherent_cols <- rep(FALSE, ncol(hoi_diagram))
  
  # Check each triplet column
  for (col in seq_len(ncol(hoi_diagram))) {
    
    # Get taxa indices with non-zero values
    taxa_idx <- which(hoi_diagram[, col] != 0)
    
    if (length(taxa_idx) != 3) {
      next  # Skip if not exactly 3 taxa
    }
    
    # Get triplet sign and taxa names
    triplet_sign <- sign(hoi_diagram[taxa_idx[1], col])
    taxa_names <- rownames(hoi_diagram)[taxa_idx]
    
    # Get pairwise signs from network_matrix
    sign_12 <- sign(network_matrix[taxa_names[1], taxa_names[2]])
    sign_23 <- sign(network_matrix[taxa_names[2], taxa_names[3]])
    sign_13 <- sign(network_matrix[taxa_names[1], taxa_names[3]])
    
    # Check strict coherence: all three pairs must match triplet sign
    if (sign_12 == triplet_sign && sign_23 == triplet_sign && sign_13 == triplet_sign) {
      coherent_cols[col] <- TRUE
    }
  }
  
  # Filter to coherent triplets only
  hoi_diagram_filtered <- hoi_diagram[, coherent_cols, drop = FALSE]
  
  if (ncol(hoi_diagram_filtered) == 0) {
    return(0)
  }
  
  return(hoi_diagram_filtered)
}