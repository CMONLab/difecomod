# difecomod: Differential Co-occurrence Analysis for Ecological Modules

**difecomod** implements the differential co-occurrence analysis methodology for inferring ecological network modules that change between environmental conditions in data sets describing presence/absence of taxa.
The workflow is designed for the analysis of clinical microbiome data sets, and is particularly useful to  understand how microbial community interactions are altered by disease states, treatments, or environmental changes.


## Installation

You can install the package difecomod from GitHub:

```r
# Install devtools if not already installed
if (!require("devtools")) install.packages("devtools")

# Install difecomod (dependencies will be installed automatically)
devtools::install_github("CMONLab/difecomod")
```

**Note**: All required dependencies (igraph, network, pheatmap, ggplot2, zCompositions, RColorBrewer, and others) will be automatically installed. For the vignette examples, you may also need to install: `SummarizedExperiment` and `tictoc`.


## Getting Started

For a comprehensive tutorial with real microbiome data examples, see the package vignette:

```r
library(difecomod)

# View the introductory vignette
vignette("introduction_difecomod", package = "difecomod")
```

The vignette demonstrates:
- Loading and processing microbiome data
- Converting abundance to occurrence matrices
- Running statistical tests (random and permutation models)
- Constructing differential co-occurrence networks
- Analyzing higher-order interactions (triplets)
- Visualizing results with heatmaps and network plots


## Workflow overview

```
┌─────────────────────────────────────────────────────────────────────┐
│                    DIFECOMOD ANALYSIS WORKFLOW                      │
└─────────────────────────────────────────────────────────────────────┘

                           INPUT DATA
                                │
                ┌───────────────┴───────────────┐
                │                               │
           Condition 1                     Condition 2
      (e.g., Healthy samples)         (e.g., Disease samples)
      Abundance Matrix (N1 × M)       Abundance Matrix (N2 × M)
                │                               │
                └───────────────┬───────────────┘
                                ↓
  ┌───────────────────────────────────────────────────────────────────────────┐
  │ STEP 1: CORE TAXA EXTRACTION                                              │
  │  • Filter by prevalence threshold (e.g., abundance >=1% in >=10% samples) │
  │  • Reduce to core microbiota (M → m taxa)                                 │
  └───────────────────────────────────────────────────────────────────────────┘
                                ↓
                ┌───────────────┴───────────────┐
                │                               │
  ┌─────────────────────────┐       ┌─────────────────────────┐
  │ STEP 2: BINARIZATION    │       │ STEP 2: BINARIZATION    │
  │  Abundance → Occurrence │       │  Abundance → Occurrence │
  │  (threshold: 0.01%)     │       │  (threshold: 0.01%)     │
  │  Binary Matrix X1       │       │  Binary Matrix X2       │
  │  (N1 × m)               │       │  (N2 × m)               │
  └───────────┬─────────────┘       └───────────┬─────────────┘
              │                                 │
              └────────────┬────────────────────┘
                           ↓
  ┌─────────────────────────────────────────────────────────────────────┐
  │ STEP 3: CO-OCCURRENCE COMPUTATION                                   │
  │  • Calculate pairwise co-occurrence frequencies                     │
  │  • O1[i,j] = freq(taxon_i & taxon_j in condition 1)                 │
  │  • O2[i,j] = freq(taxon_i & taxon_j in condition 2)                 │
  │  • ΔO[i,j] = O1[i,j] - O2[i,j]                                      │
  └─────────────────────────────────────────────────────────────────────┘
                                ↓
  ┌─────────────────────────────────────────────────────────────────────┐
  │ STEP 4: NULL MODEL GENERATION (n = 5000-8000 iterations)            │
  ├─────────────────────────────────────────────────────────────────────┤
  │  ┌──────────────────────────┐    ┌──────────────────────────┐       │
  │  │  RANDOM MODEL            │    │  PERMUTATION MODEL       │       │
  │  │  • Shuffle columns       │    │  • Shuffle sample labels │       │
  │  │    independently         │    │    between conditions    │       │
  │  │  • Test: random          │    │  • Test: arbitrary       │       │
  │  │    co-occurrence         │    │    grouping              │       │
  │  │  • Preserves marginal    │    │  • Preserves total       │       │
  │  │    frequencies           │    │    sample structure      │       │
  │  └──────────────────────────┘    └──────────────────────────┘       │
  └─────────────────────────────────────────────────────────────────────┘
                                ↓
  ┌─────────────────────────────────────────────────────────────────────┐
  │ STEP 5: STATISTICAL THRESHOLDS                                      │
  │  • Random model: μ ± σ (or n_sigma × σ)                             │
  │  • Permutation model: 5th & 95th percentiles (conf = 0.05)          │
  │  • Significant if BOTH criteria exceeded                            │
  └─────────────────────────────────────────────────────────────────────┘
                                ↓
  ┌─────────────────────────────────────────────────────────────────────┐
  │ STEP 6: NETWORK CONSTRUCTION                                        │
  │  • Nodes = Taxa with significant differential co-occurrences        │
  │  • Edges = Significant ΔO values                                    │
  │  • Edge sign:                                                       │
  │    (+) = Over-represented in Condition 1                            │
  │    (−) = Under-represented in Condition 1                           │
  └─────────────────────────────────────────────────────────────────────┘
                                ↓
  ┌─────────────────────────────────────────────────────────────────────┐
  │ STEP 7: HIGHER-ORDER INTERACTIONS (HOI)                             │
  │  • Identify triangles in network (3-node cliques)                   │
  │  • Calculate three-way co-occurrence: O[i,j,k]                      │
  │  • Generate null models for triplets                                │
  │  • Test significance of Δ-HOI                                       │
  │  • Construct HOI diagram                                            │
  └─────────────────────────────────────────────────────────────────────┘
                                ↓
  ┌─────────────────────────────────────────────────────────────────────┐
  │ STEP 8: VISUALIZATION                                               │
  ├─────────────────────────────────────────────────────────────────────┤
  │  ┌──────────────────────┐  ┌──────────────────────┐                 │
  │  │ Network Plot         │  │ Abundance Heatmap    │                 │
  │  │ • Differential       │  │ • CLR-transformed    │                 │
  │  │   co-occurrence      │  │ • Sample clustering  │                 │
  │  │   network            │  │ • Condition labels   │                 │
  │  └──────────────────────┘  └──────────────────────┘                 │
  │  ┌──────────────────────────────────────────────┐                   │
  │  │ HOI Heatmap                                  │                   │
  │  │ • Triplet interactions                       │                   │
  │  │ • Green = enriched in Condition 1            │                   │
  │  │ • Red = enriched in Condition 2              │                   │
  │  └──────────────────────────────────────────────┘                   │
  └─────────────────────────────────────────────────────────────────────┘
                                ↓
                          OUTPUT RESULTS
                ┌───────────────┴───────────────┐
                │                               │
        Network Matrix                   HOI Diagram
      (Significant ΔO)              (Significant Δ-HOI)
```

**Key Concepts:**
- **ΔO** = Differential co-occurrence (Condition 1 - Condition 2)
- **Dual statistical criteria** = Both random AND permutation models must be significant
- **HOI** = Higher-Order Interactions (three-taxa co-occurrences)
- **m** = Number of core taxa after filtering
- **N1, N2** = Number of samples in each condition

## Best Practices

When working with real microbiome data:

1. **Sample Size**: Ensure adequate sample sizes (>20 samples per condition recommended).
2. **Core Extraction**: Use appropriate prevalence thresholds to extract core genera (e.g., abundance greater or equal to 1%  in at least 10% of samples).
3. **Binarization**: Set meaningful detection thresholds (default 0.01% relative abundance).
4. **Null Models Generation**: Use sufficient number of permutations/shufflings for robust statistical inference (≥5000 recommended).
5. **Biological Interpretation**: Validate significant interactions using biological knowledge and literature.


## Reference

If you use difecomod in your research, please cite:

**Differential Co-occurrence Analysis: a Method to Extract Ecological Modules from Clinical Microbiome Data**. J. Iacovacci, N. Cannon, J.A. McCulloch, T. Rancati, G. Trinchieri (2025).

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

- **Author**: Jacopo Iacovacci
- **Email**: jacopo.iacovacci@istitutotumori.mi.it
- **Issues**: [GitHub Issues](https://github.com/CMONLab/difecomod/issues)