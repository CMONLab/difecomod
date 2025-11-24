#' Visualization Functions for Differential Co-occurrence Analysis
#'
#' Functions for creating plots and visualizations of differential co-occurrence
#' networks and results.
#'
#' @name visualization
NULL

#' Visualize Core Microbiota from Abundance Data 
#'
#' Plot taxa prevalence across samples
#'
#' @param abundance_matrix Numeric matrix of relative abundances (taxa x samples)
#' @param prevalence_threshold Numeric threshold (abundance percentage) for determining taxa prevalence (default: 0.1%)
#' @param min_prevalence Minimum prevalence (fraction of samples) required to retain taxa (default: 0.1)
#'
#' @return ggplot object
#'
#' @examples
#' # Example requires a completed analysis
#' # See package vignette for full workflow
#'
#' @importFrom ggpubr ggbarplot
#'
#' @export
core_visualization <- function(abundance_matrix, prevalence_threshold = 1, min_prevalence = 0.1) {
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
  
  # Number of samples is number of columns!
  n_samples <- ncol(abundance_matrix)
  
  # Calculate prevalence per taxon (number of samples where abundance >= threshold)
  prevalence <- rowSums(perc_table >= prevalence_threshold)
  
  # Create a tidy data frame for plotting
  df.taxa <- data.frame(taxa = rownames(abundance_matrix), prevalence = prevalence)
  df.taxa <- df.taxa[df.taxa$prevalence > 4, ] # keep only taxa present in at least 5 samples
  df.taxa$taxa <- factor(df.taxa$taxa, levels = df.taxa$taxa[order(df.taxa$prevalence, decreasing = TRUE)])
  
  # Set the cutoff line for core
  cutoff <- round(n_samples * min_prevalence)
  df.taxa$color <- ifelse(df.taxa$prevalence >= cutoff, "darkgreen", "darkgray")
  
  # Build the plot
  x_pos <- nrow(df.taxa)
  y_pos <- cutoff + max(1, 0.015 * max(df.taxa$prevalence))
  
  p <- ggpubr::ggbarplot(
    df.taxa, x = 'taxa', y = 'prevalence',
    fill = "color",
    palette = c("darkgrey", "darkgreen"),
    x.text.angle = 90,
    ylab = sprintf("Numb. of samples w.t. abund. >=%s%%", prevalence_threshold),
    xlab = FALSE,
    title = sprintf('Core Microbiota (n=%s patients)', n_samples)
  ) +
    ggpubr::font("ylab", size = 14, face = "bold") +
    ggpubr::font("xy.text", size = 9, face = "bold") +
    ggplot2::geom_hline(yintercept = cutoff, linetype = 2, color = "black") +
    ggplot2::annotate(
      "text",
      x = x_pos + 0.2,
      y = y_pos,
      label = sprintf('%s%% population', min_prevalence * 100),
      hjust = + 2,
      vjust = -0.5,
      size = 4.5,
      fontface = "italic"
    ) +
    ggplot2::theme(legend.position = "none")
  
  return(p)
}

#' Microbiota Heatmap
#'
#' Plot the microbiota clr-transformed abundance profiles across samples
#'
#' @param scaled_clr_table Numeric matrix of standardized clr-transformed abundance profiles (taxa x samples)
#' @param covariates Numeric matrix of sample(patients) covariates (samples x covariate)
#'
#' @return ggplot object
#'
#' @examples
#' # Example requires a completed analysis
#' # See package vignette for full workflow
#'
#'@importFrom pheatmap pheatmap
#'
#' @export
microbiota_heatmap <- function(scaled_clr_table=NULL, covariates=NULL) {
  # Validate inputs
  if (is.null(scaled_clr_table)) {
    stop("scaled_clr_table must be provided")
  }

  # Use pheatmap for clustering and plotting
  phm <- pheatmap::pheatmap(
    mat = scaled_clr_table,
    breaks = seq(-max(abs(min(scaled_clr_table)),abs(max(scaled_clr_table))), max(abs(min(scaled_clr_table)),abs(max(scaled_clr_table))),length.out=(100)),
    color=grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(99),
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    clustering_distance_rows = "euclidean",
    treeheight_row = 0,
    clustering_distance_cols = "euclidean",
    clustering_method = "ward.D2",
    annotation_col = covariates,
    scale = 'none',
    fontsize_row = 6,
    show_colnames = FALSE,
  )

  return(phm)
}

#' Plot Co-occurrence Network
#'
#' Creates a network visualization of differential co-occurrences using ggplot2 and network packages.
#'
#' @param network_matrix Adjacency matrix from construct_network_delta
#' @param layout Network layout algorithm (default: "fruchtermanreingold")
#' @param node_size_var Variable for node sizing (default: "degree")
#' @param edge_colors Colors for positive and negative edges (default: c("green", "red"))
#' @param node_labels Logical, whether to show node labels (default: TRUE)
#' @param title Plot title
#'
#' @return ggplot object
#'
#' @examples
#' # Example requires a completed analysis
#' # See package vignette for full workflow
#'
#' @importFrom network network
#' @importFrom graphics par
#' @export
plot_cooccurrence_network <- function(network_matrix, layout = "fruchtermanreingold",
                                      node_size_var = "degree", 
                                      edge_colors = c("green", "red"),
                                      node_labels = TRUE, title = "Differential Co-occurrence Network") {
  
  # Validate inputs
  if (!is.matrix(network_matrix)) {
    stop("network_matrix must be a matrix")
  }
  if (nrow(network_matrix) != ncol(network_matrix)) {
    stop("network_matrix must be square")
  }
  if (length(edge_colors) != 2) {
    stop("edge_colors must be a vector of length 2")
  }
  
  # Remove diagonal for network creation
  net_for_plotting <- network_matrix
  diag(net_for_plotting) <- 0
  
  # Create network object
  if (requireNamespace("network", quietly = TRUE) && requireNamespace("GGally", quietly = TRUE)) {
    
    net <- network::network(net_for_plotting, directed = FALSE, ignore.eval = FALSE, names.eval = "weights")
    network::set.edge.attribute(net, "size", abs(network::get.edge.attribute(net, "weights")))
    network::set.vertex.attribute(net, "degree", as.numeric(rowSums(abs(network_matrix))))
    network::set.edge.attribute(net, "color", ifelse(network::get.edge.attribute(net, "weights") > 0, edge_colors[1], edge_colors[2]))
    ##################
    # ISSUE: potential issue when network has edges of one type only
    # edge_weights <- network::get.edge.attribute(net, "weights")
    # colors_used <- unique(ifelse(edge_weights > 0, edge_colors[1], edge_colors[2]))
    # # If only positive or negative weights are found
    # if (length(colors_used) == 1) {
    #   network::set.edge.attribute(net, "color", colors_used)
    # }else{
    #   network::set.edge.attribute(net, "color", 
    #                             ifelse(network::get.edge.attribute(net, "weights") > 0, edge_colors[1], edge_colors[2]))
    # }
    ##################
    # Create plot
    p <- GGally::ggnet2(net, 
                        label = node_labels, 
                        node.size = node_size_var, 
                        edge.size = "size", 
                        edge.color = "color",
                        layout.exp = 0.2) +
      ggplot2::ggtitle(title) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    
    return(p)
    
  } else {
    stop("Packages 'network' and 'GGally' are required for network plotting")
  }
}
