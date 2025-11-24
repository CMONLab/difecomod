---
title: "difecomod: Differential Co-occurrence Analysis for Ecological Modules"
author: "Jacopo Iacovacci"
date: "2025-11-24"
output:
  rmarkdown::html_vignette:
    toc: true
    toc_depth: 2
vignette: >
  %\VignetteIndexEntry{Differential Co-occurrence Analysis for Ecological Modules}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



## Introduction

The `difecomod` package implements differential co-occurrence analysis methodology for detecting ecosystem blocks or modules that change between host conditions in microbiome data. This approach is particularly useful for understanding how microbial community interactions are altered by disease states, treatments, or environmental changes.

### Methodology Overview

Given the modularity of microbial network interactions, we hypothesize that different ecosystem blocks or modules coexist within the microbiota and that changes in the host environment (e.g., disease conditions) might favor certain blocks of interacting species over others.

The methodology consists of several key steps:

1.  **Data Processing**: Convert relative abundance data to binary occurrence matrices
2.  **Co-occurrence Calculation**: Compute frequencies of taxa co-occurrence in each condition\
3.  **Statistical Testing**: Use null random and permutation models for co-occurrences to assess significance
4.  **Network Construction**: Build microbial networks of significant differential co-occurrences
5.  **Higher-Order Interactions**: Analyze three-taxa co-occurrences in data

## Installation


``` r
# Install devtools if not already installed
if (!require("devtools")) install.packages("devtools")

# Install difecomod from GitHub
devtools::install_github("CMONLab/difecomod")

# Load the package
library(difecomod)
```


``` r
library(difecomod)
#> Error in library(difecomod): there is no package called 'difecomod'
library(SummarizedExperiment)
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: 'MatrixGenerics'
#> The following objects are masked from 'package:matrixStats':
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> Loading required package: GenomicRanges
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> 
#> Attaching package: 'BiocGenerics'
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:base':
#> 
#>     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
#>     as.data.frame, basename, cbind, colnames, dirname, do.call,
#>     duplicated, eval, evalq, get, grep, grepl, intersect, is.unsorted,
#>     lapply, mapply, match, mget, order, paste, pmax, pmax.int, pmin,
#>     pmin.int, rank, rbind, rownames, sapply, setdiff, sort, table,
#>     tapply, union, unique, unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: 'S4Vectors'
#> The following object is masked from 'package:utils':
#> 
#>     findMatches
#> The following objects are masked from 'package:base':
#> 
#>     I, expand.grid, unname
#> Loading required package: IRanges
#> Loading required package: GenomeInfoDb
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: 'Biobase'
#> The following object is masked from 'package:MatrixGenerics':
#> 
#>     rowMedians
#> The following objects are masked from 'package:matrixStats':
#> 
#>     anyMissing, rowMedians
library(ggplot2)
library(pheatmap)
library(igraph)
#> 
#> Attaching package: 'igraph'
#> The following object is masked from 'package:GenomicRanges':
#> 
#>     union
#> The following object is masked from 'package:IRanges':
#> 
#>     union
#> The following object is masked from 'package:S4Vectors':
#> 
#>     union
#> The following objects are masked from 'package:BiocGenerics':
#> 
#>     normalize, path, union
#> The following objects are masked from 'package:stats':
#> 
#>     decompose, spectrum
#> The following object is masked from 'package:base':
#> 
#>     union
library(tictoc)
#> Warning: package 'tictoc' was built under R version 4.3.2
#> 
#> Attaching package: 'tictoc'
#> The following object is masked from 'package:SummarizedExperiment':
#> 
#>     shift
#> The following object is masked from 'package:GenomicRanges':
#> 
#>     shift
#> The following object is masked from 'package:IRanges':
#> 
#>     shift
set.seed(1234)  # For reproducible examples
```

## Basic Usage

### Using Gupta Data Set as Example Data

For demonstration, we'll use microbiome data from faecal samples in two host conditions: colorectal cancer (CRC) versus healthy individuals:


``` r
# Load JAMS processed data from package
# The data file is located in the package's extdata directory
data_path <- system.file("extdata", "CRC_cohorts_JAMS_processed", "Gupta.RData",
                         package = "difecomod")

# Check if data file exists
if (data_path == "" || !file.exists(data_path)) {
  stop("Data file not found. Please ensure the package was installed correctly with data files.")
}
#> Error: Data file not found. Please ensure the package was installed correctly with data files.

load(data_path)
#> Warning in readChar(con, 5L, useBytes = TRUE): cannot open compressed file '',
#> probable reason 'No such file or directory'
#> Error in readChar(con, 5L, useBytes = TRUE): cannot open the connection

# Extract abundance table and metadata
count_table <- assay(expvec$LKT, "BaseCounts");
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'assay': object 'expvec' not found
pheno <- as.data.frame(colData(expvec$LKT));
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'as.data.frame': error in evaluating the argument 'x' in selecting a method for function 'colData': object 'expvec' not found
pheno$condition <- pheno$disease;
#> Error: object 'pheno' not found

# Calculate total read counts
tot_counts <- colSums(count_table);
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'colSums': object 'count_table' not found

# Calculate percentage tables
abundance_table <- scale(count_table, center = FALSE, scale = tot_counts);
#> Error: object 'count_table' not found

# Visualize core microbiota
core_plot <- core_visualization(  
  abundance_table,
  prevalence_threshold = 1,
  min_prevalence = 0.1
);
#> Error in core_visualization(abundance_table, prevalence_threshold = 1, : could not find function "core_visualization"

print(core_plot)
#> Error: object 'core_plot' not found

# Extract core microbiota taxa
abundance_table <- core_extraction(
  abundance_table,
  prevalence_threshold = 1,
  min_prevalence = 0.1
);
#> Error in core_extraction(abundance_table, prevalence_threshold = 1, min_prevalence = 0.1): could not find function "core_extraction"

# Separate abundance tables for healthy and disease 
samples_healthy <- rownames(pheno[pheno$condition=='healthy',]);
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'rownames': object 'pheno' not found
samples_disease <- rownames(pheno[pheno$condition=='CRC',]);
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'rownames': object 'pheno' not found

abundance_healthy <- abundance_table[, samples_healthy];
#> Error: object 'abundance_table' not found
abundance_disease <- abundance_table[, samples_disease];
#> Error: object 'abundance_table' not found

# Display data summaries
cat("Healthy condition - Samples:", ncol(abundance_healthy), "Taxa:", nrow(abundance_healthy), "\n")
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'ncol': object 'abundance_healthy' not found
cat("Disease condition - Samples:", ncol(abundance_disease), "Taxa:", nrow(abundance_disease), "\n")
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'ncol': object 'abundance_disease' not found
```

## Step-by-Step Differential Co-occurrence Analysis

For more control over the analysis process, you can run each step individually:

### Step 1: Convert to Occurrence Matrices


``` r
# Convert abundance to binary occurrence matrices
occ_healthy <- calculate_cooccurrence_matrix(
  abundance_healthy,
  detection_threshold = 0.01
)
#> Error in calculate_cooccurrence_matrix(abundance_healthy, detection_threshold = 0.01): could not find function "calculate_cooccurrence_matrix"

occ_disease <- calculate_cooccurrence_matrix(
  abundance_disease,
  detection_threshold = 0.01
)
#> Error in calculate_cooccurrence_matrix(abundance_disease, detection_threshold = 0.01): could not find function "calculate_cooccurrence_matrix"

cat("Occurrence matrices created:\n")
#> Occurrence matrices created:

## Ensure same taxa in both matrices
#common_taxa <- sort(intersect(colnames(occ_healthy), colnames(occ_disease)))
#occ_healthy <- occ_healthy[, common_taxa]
#occ_disease <- occ_disease[, common_taxa]

#cat("Common taxa for analysis:", length(common_taxa), "\n")
#print(common_taxa)
```

### Step 2: Calculate Co-occurrence Frequencies


``` r
# Extract co-occurrence frequencies for each condition
cooccur_healthy <- extract_cooccurrence_list(occ_healthy)
#> Error in extract_cooccurrence_list(occ_healthy): could not find function "extract_cooccurrence_list"
cooccur_disease <- extract_cooccurrence_list(occ_disease)
#> Error in extract_cooccurrence_list(occ_disease): could not find function "extract_cooccurrence_list"

# Calculate delta co-occurrence
delta_cooccur <- cooccur_healthy$O - cooccur_disease$O
#> Error: object 'cooccur_healthy' not found

cat("Co-occurrence statistics:\n")
#> Co-occurrence statistics:
cat("Number of taxa pairs:", length(delta_cooccur), "\n")
#> Error: object 'delta_cooccur' not found
cat("Delta co-occurrence (Healthy - Disease) frequency range:", 
    sprintf("[%.4f, %.4f]", min(delta_cooccur), max(delta_cooccur)), "\n")
#> Error: object 'delta_cooccur' not found
```

### Step 3: Statistical Testing


``` r
tic()
# Random delta co-occurrence model - tests against chance co-occurrence
cat("Running random delta co-occurrence model...\n")
#> Running random delta co-occurrence model...
random_model <- delta_cooccurrence_random_model(
  occ_healthy, occ_disease, 
  n = 8000
)
#> Error in delta_cooccurrence_random_model(occ_healthy, occ_disease, n = 8000): could not find function "delta_cooccurrence_random_model"
toc()
#> 0.008 sec elapsed
tic()
# Permutation delta co-occurrence model - tests against random sample grouping  
cat("Running permutation delta co-occurrence model...\n")
#> Running permutation delta co-occurrence model...
perm_model <- delta_cooccurrence_permutation_model(
  occ_healthy, occ_disease,
  n = 8000
)
#> Error in delta_cooccurrence_permutation_model(occ_healthy, occ_disease, : could not find function "delta_cooccurrence_permutation_model"
toc()
#> 0.009 sec elapsed
# Calculate statistical thresholds
upper_lower_rand <- compute_random_model_thresholds(random_model$Values, n_sigma = 1)
#> Error in compute_random_model_thresholds(random_model$Values, n_sigma = 1): could not find function "compute_random_model_thresholds"
upper_lower_perm <- compute_permutation_model_thresholds(perm_model$Values, conf = 0.05)
#> Error in compute_permutation_model_thresholds(perm_model$Values, conf = 0.05): could not find function "compute_permutation_model_thresholds"
```

### Step 4: Network Construction


``` r
# Construct differential co-occurrence network
network_matrix <- construct_network_delta(
  delta_I = delta_cooccur,
  links_I = cooccur_healthy$I, 
  nodes_I = colnames(occ_healthy), #common_taxa,
  stat_thr_rand = upper_lower_rand,
  stat_thr_perm = upper_lower_perm
)
#> Error in construct_network_delta(delta_I = delta_cooccur, links_I = cooccur_healthy$I, : could not find function "construct_network_delta"

print(sort(colSums(abs(network_matrix))))
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'sort': error in evaluating the argument 'x' in selecting a method for function 'colSums': object 'network_matrix' not found

cat("Network construction results:\n")
#> Network construction results:
if (length(network_matrix) > 1) {
  n_significant <- sum(abs(network_matrix) > 0) / 2  # Symmetric matrix
  n_positive <- sum(network_matrix > 0) / 2
  n_negative <- sum(network_matrix < 0) / 2
  
  cat("Significant interactions:", n_significant, "\n")
  cat("Over-represented in Healthy VS Disease:", n_positive, "\n") 
  cat("Under-represented in Healthy VS Disease:", n_negative, "\n")
  cat("Connected taxa:", nrow(network_matrix), "\n")
  
  network_plot <- plot_cooccurrence_network(network_matrix, title = "Differential Co-occurrence Network")
  print(network_plot)
  # To export in pdf
  # ggsave(sprintf("differential_cooccurrence_network.pdf"), plot = network_plot, device = NULL,
  #       path = "./", scale = 0.9, width = 20, height = 16, dpi = 300,
  #       limitsize = TRUE);
  
} else {
  cat("No significant interactions found.\n")
}
#> Error: object 'network_matrix' not found
```

## Abundance Data Visualization

Even when differential co-occurrence patterns are subtle, we can visualize the abundance patterns:


``` r
#  Select taxa in Differential Co-occurrence Network
dcn_taxa <- names(rowSums(abs(network_matrix))[rowSums(abs(network_matrix))>0]);
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'rowSums': object 'network_matrix' not found

# Centered-log ratio transformation of core count abundances
clr_table <- clr_transformation(
  count_table[rownames(abundance_table), ],
  zero_imputation_method = "GBM"
);
#> Error in clr_transformation(count_table[rownames(abundance_table), ], : could not find function "clr_transformation"

# Standardization of clr-table
z <- t(scale(t(clr_table), center = TRUE, scale = TRUE))
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 't': error in evaluating the argument 'x' in selecting a method for function 't': object 'clr_table' not found

# Extract covariate data
col_covariates <- data.frame(pheno[, c("condition", "gender")])
#> Error: object 'pheno' not found
col_covariates$gender <- factor(col_covariates$gender, levels=c("male", "female"))
#> Error: object 'col_covariates' not found
col_covariates$condition <- factor(col_covariates$condition, levels=c("healthy", "CRC"))
#> Error: object 'col_covariates' not found

heatmap_plot <- microbiota_heatmap(z, col_covariates)
#> Error in microbiota_heatmap(z, col_covariates): could not find function "microbiota_heatmap"

# To export in pdf
# ggsave(sprintf("TEST_differential_cooccurrence_heatmap.pdf"), plot = heatmap_plot, device = NULL,
#       path = "./", scale = 0.9, width = 8, height = 8, dpi = 300,
#       limitsize = TRUE);

```

## Higher-Order Interactions Analysis

The package supports analysis of three-taxa co-occurrences:


``` r
# Extract Higher-Order Interactions (HOI) Analysis
# ================================================
# Higher-order interactions capture complex relationships between three taxa
# that co-occur together more or less frequently than expected.

cat("Analyzing higher-order interactions...\n\n")
#> Analyzing higher-order interactions...

# Extract all triplet interactions from both conditions
hoi_healthy <- extract_hoi_list(occ_healthy)
#> Error in extract_hoi_list(occ_healthy): could not find function "extract_hoi_list"
hoi_disease <- extract_hoi_list(occ_disease)
#> Error in extract_hoi_list(occ_disease): could not find function "extract_hoi_list"

# Calculate delta-HOI for all possible triplets
delta_hoi_all <- hoi_healthy$O - hoi_disease$O
#> Error: object 'hoi_healthy' not found

cat("Number of all possible triplet interactions:", length(hoi_healthy$O), "\n")
#> Error: object 'hoi_healthy' not found

# Identify network triangles (3-node cliques)
# ===================================================
# Look for triangles in the co-occurrence network (motifs where 
# three taxa form a connected cluster) with edges having consistent signs).

# Extract network triangles with consistent edge values 
cat("Identifying network triangles with consistent edge signs...\n")
#> Identifying network triangles with consistent edge signs...

triangle_results <- extract_network_triangles(
  adj_matrix = network_matrix,
  filter_sign = TRUE
)
#> Error in extract_network_triangles(adj_matrix = network_matrix, filter_sign = TRUE): could not find function "extract_network_triangles"

#INetF <- graph_from_adjacency_matrix(abs(network_matrix), mode="undirected");
#triangles_nodes <- count_triangles(INetF);
#triangles_list <- t(matrix(triangles(INetF), nrow = 3));
#triangles_list <- t(apply(triangles_list, 1, sort));
  
triangles_list <- triangle_results$triangles_list
#> Error: object 'triangle_results' not found
triangle_taxa_indices <- triangle_results$triangle_taxa_indices
#> Error: object 'triangle_results' not found
triangle_taxa <- triangle_results$triangle_taxa
#> Error: object 'triangle_results' not found

table(c(triangles_list[,1],triangles_list[,2],triangles_list[,3]))
#> Error: object 'triangles_list' not found

cat("Network triangles identified:", nrow(triangles_list), "\n")
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'nrow': object 'triangles_list' not found

# # Get unique taxa involved in triangles
cat("Unique taxa in triangles:", length(triangle_taxa), "\n")
#> Error: object 'triangle_taxa' not found
cat("Triangle taxa:\n")
#> Triangle taxa:
print(head(triangle_taxa, 10))
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'head': object 'triangle_taxa' not found
cat("  ... (showing first 10)\n\n")
#>   ... (showing first 10)

# Extract HOI specifically for network triangles with remapped indices
hoi_healthy_tri <- extract_hoi_sublist(occ_healthy, sub_list = triangles_list)
#> Error in extract_hoi_sublist(occ_healthy, sub_list = triangles_list): could not find function "extract_hoi_sublist"
hoi_disease_tri <- extract_hoi_sublist(occ_disease, sub_list = triangles_list)
#> Error in extract_hoi_sublist(occ_disease, sub_list = triangles_list): could not find function "extract_hoi_sublist"
delta_hoi_triangles <- hoi_healthy_tri$O - hoi_disease_tri$O
#> Error: object 'hoi_healthy_tri' not found

cat("Differential HOI statistics for triangles (Healthy - Disease):\n")
#> Differential HOI statistics for triangles (Healthy - Disease):
cat("  Mean change:", round(mean(delta_hoi_triangles), 3), "\n")
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'mean': object 'delta_hoi_triangles' not found
cat("  Range:", round(min(delta_hoi_triangles), 3), "to", 
    round(max(delta_hoi_triangles), 3), "\n\n")
#> Error: object 'delta_hoi_triangles' not found

# Generate null models for statistical testing
# =====================================================
#     Random Model: Shuffles feature columns independently within each condition
#     Tests whether differential HOI exceeds what's expected if taxa associations
#     were random but marginal frequencies were preserved.

tic()
cat("Computing random model null distribution...\n")
#> Computing random model null distribution...
hoi_random <- delta_hoi_random_model(
  X1occ = occ_healthy, 
  X2occ = occ_disease,
  n = 8000, 
  n_sigma = 1,
  triangle_list = triangles_list,
  values = TRUE
)
#> Error in delta_hoi_random_model(X1occ = occ_healthy, X2occ = occ_disease, : could not find function "delta_hoi_random_model"
toc()
#> 0.043 sec elapsed

#     Permutation Model: Permutes sample labels between conditions
#     Tests whether differential HOI exceeds what's expected if condition
#     assignment were arbitrary.
tic()
cat("Computing permutation model null distribution...\n")
#> Computing permutation model null distribution...
hoi_perm <- delta_hoi_permutation_model(
  X1occ = occ_healthy,
  X2occ = occ_disease,
  n = 8000,
  conf = 0.05,
  triangle_list = triangles_list,
  values = TRUE
)
#> Error in delta_hoi_permutation_model(X1occ = occ_healthy, X2occ = occ_disease, : could not find function "delta_hoi_permutation_model"
toc()
#> 0.434 sec elapsed
```

# Visualize Higher-Order Interactions


``` r
# This creates a matrix showing which triplets are significantly enriched/depleted
# based on both null models. The diagram visualizes the overlap between
# statistical tests and highlights robust HOI changes.

cat("Constructing HOI diagram...\n")
#> Constructing HOI diagram...

# Extract threshold statistics from null models
# The construct_HOI_diagram function requires thresholds in list format
stat_thr_rand <- list(
  lower = hoi_random$Stats$low_thr,
  upper = hoi_random$Stats$up_thr
)
#> Error: object 'hoi_random' not found

stat_thr_perm <- list(
  lower = hoi_perm$Stats$low_thr,
  upper = hoi_perm$Stats$up_thr
)
#> Error: object 'hoi_perm' not found

# Construct the HOI diagram
hoi_diagram_raw <- construct_hoi_diagram(
  delta_I = delta_hoi_triangles,
  links_I = hoi_healthy_tri$I,
  nodes_I = triangle_taxa,
  stat_thr_rand = stat_thr_rand,
  stat_thr_perm = stat_thr_perm
)
#> Error in construct_hoi_diagram(delta_I = delta_hoi_triangles, links_I = hoi_healthy_tri$I, : could not find function "construct_hoi_diagram"

hoi_diagram <- validate_hoi_diagram(
  hoi_diagram = hoi_diagram_raw,
  network_matrix = network_matrix
)
#> Error in validate_hoi_diagram(hoi_diagram = hoi_diagram_raw, network_matrix = network_matrix): could not find function "validate_hoi_diagram"

# Visualize if significant triplets exist
if (!is.null(dim(hoi_diagram)) && is.matrix(hoi_diagram)) {
  cat("HOI diagram created with", nrow(hoi_diagram), "taxa in", 
      ncol(hoi_diagram), "significant triplets\n")
  cat("Generating visualization...\n\n")
  
  # Create heatmap visualization
  # Red = depleted in healthy (enriched in disease)
  # Green = enriched in healthy (depleted in disease)
  heatmap_hoi <- pheatmap::pheatmap(
    hoi_diagram,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_colnames = FALSE,
    color = c("red", "black", "black", "green"),
    breaks = c(-1, -0.5, 0, 0.5, 1),
    main = "Higher-Order Differential Co-occurrences",
    fontsize = 10
  )

  # # To export in pdf
  # ggsave(sprintf("TEST_differential_cooccurrence_hoi.pdf"), plot = heatmap_hoi, device = NULL,
  #        path = "./", scale = 0.9, width = 8, height = 8, dpi = 300, limitsize = TRUE);

  cat("HOI diagram visualization displayed above.\n")
  cat("  - Green: HOI triplets more frequent in healthy condition\n")
  cat("  - Red: HOI triplets more frequent in disease condition\n")
  cat("  - Black: Non-significant changes\n\n")
  
} else if (identical(hoi_diagram, 0)) {
  # Function returned 0 (no significant HOI)
  cat("No significant HOI changes detected.\n")
  cat("This may indicate:\n")
  cat("  - Similar interaction patterns between conditions\n")
  cat("  - Insufficient statistical power (increase n permutations)\n")
  cat("  - Very few network triangles in the data\n\n")
  
} else {
  cat("HOI diagram could not be created.\n")
}
#> Error: object 'hoi_diagram' not found

# HOI Summary Statistics
# ======================
cat("\nSummary of significant HOI changes:\n")
#> 
#> Summary of significant HOI changes:

# Count significant triplets by each model
sig_random <- which(
  delta_hoi_triangles > hoi_random$Stats$up_thr | 
  delta_hoi_triangles < hoi_random$Stats$low_thr
)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'which': object 'delta_hoi_triangles' not found

sig_perm <- which(
  delta_hoi_triangles > hoi_perm$Stats$up_thr | 
  delta_hoi_triangles < hoi_perm$Stats$low_thr
)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'which': object 'delta_hoi_triangles' not found

sig_both <- intersect(sig_random, sig_perm)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'intersect': object 'sig_random' not found

cat("  Significant by random model:", length(sig_random), "triplets\n")
#> Error: object 'sig_random' not found
cat("  Significant by permutation model:", length(sig_perm), "triplets\n")
#> Error: object 'sig_perm' not found
cat("  Significant by both models:", length(sig_both), "triplets\n")
#> Error: object 'sig_both' not found

if (length(sig_both) > 0) {
  cat("    - Enriched in healthy:", 
      sum(delta_hoi_triangles[sig_both] > 0), "triplets\n")
  cat("    - Enriched in disease:", 
      sum(delta_hoi_triangles[sig_both] < 0), "triplets\n")
} else {
  cat("  No triplets passed both statistical tests.\n")
}
#> Error: object 'sig_both' not found

cat("\nHOI visualization section completed.\n")
#> 
#> HOI visualization section completed.
```

## Real Data Considerations

When working with real microbiome data:

1.  **Sample Size**: Ensure adequate sample sizes (\>20 samples per condition recommended)

2.  **Taxa Filtering**: Use appropriate prevalence thresholds to extract core genera (10% recommended)

3.  **Abundance Threshold**: Set meaningful detection thresholds (0.01-0.1% relative abundance)

4.  **Permutation/Reshuffling Number**: Use sufficient number of permutations/reshufflings for robust statistical inference (\>=5000)

5.  **Biological Interpretation**: Validate significant interactions using biological knowledge and literature

## Session Information


``` r
sessionInfo()
#> R version 4.3.1 (2023-06-16)
#> Platform: x86_64-apple-darwin20 (64-bit)
#> Running under: macOS Monterey 12.6.7
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRblas.0.dylib 
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
#> 
#> locale:
#> [1] C/UTF-8/C/C/C/C
#> 
#> time zone: Europe/London
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] tictoc_1.2.1                igraph_1.5.0               
#>  [3] pheatmap_1.0.12             ggplot2_3.4.2              
#>  [5] SummarizedExperiment_1.30.2 Biobase_2.60.0             
#>  [7] GenomicRanges_1.52.0        GenomeInfoDb_1.36.1        
#>  [9] IRanges_2.34.1              S4Vectors_0.38.1           
#> [11] BiocGenerics_0.46.0         MatrixGenerics_1.12.2      
#> [13] matrixStats_1.0.0          
#> 
#> loaded via a namespace (and not attached):
#>  [1] Matrix_1.6-0            gtable_0.3.3            dplyr_1.1.2            
#>  [4] compiler_4.3.1          crayon_1.5.2            tidyselect_1.2.0       
#>  [7] bitops_1.0-7            scales_1.2.1            lattice_0.21-8         
#> [10] R6_2.5.1                XVector_0.40.0          generics_0.1.3         
#> [13] S4Arrays_1.0.4          knitr_1.50              tibble_3.2.1           
#> [16] DelayedArray_0.26.6     munsell_0.5.0           GenomeInfoDbData_1.2.10
#> [19] RColorBrewer_1.1-3      pillar_1.9.0            rlang_1.1.6            
#> [22] utf8_1.2.3              xfun_0.52               cli_3.6.5              
#> [25] withr_3.0.2             magrittr_2.0.3          zlibbioc_1.46.0        
#> [28] grid_4.3.1              lifecycle_1.0.4         vctrs_0.6.3            
#> [31] evaluate_1.0.4          glue_1.6.2              RCurl_1.98-1.12        
#> [34] fansi_1.0.4             colorspace_2.1-0        pkgconfig_2.0.3        
#> [37] tools_4.3.1
```
