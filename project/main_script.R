####--LIBRARIES--####

library(dplyr)
library(tidyverse)
library(ggplot2)
library(spacexr) #devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
library(nicheDE) #devtools::install_github("kaishumason/NicheDE")
library(spatialGE) #devtools::install_github("fridleylab/spatialGE")
library(Giotto) #pak::pkg_install("drieslab/Giotto")


####--DIRECTORIES--####

dir <- "./project/material/GSE165098"
sample.ids <- c("LPS1", "LPS2", "LPS3", "saline1", "saline2", "saline3")
sample.path <- file.path(dir, sample.ids)


####--SAMPLE 1--####

# Obtain data
s1 <- createGiottoVisiumObject(
  visium_dir = sample.path[1], #LPS1
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = createGiottoInstructions(
    save_dir = "./project/outcomes/LPS1",
    save_plot = T, show_plot = F))

s1.matrix <- s1@expression[["cell"]][["rna"]][["raw"]]@exprMat # Extract expression matrix

# Visualize spots
s1.spots <- spatPlot(s1, point_alpha = 0.6, 
                     default_save_name = "spots") # Image not properly aligned
s1.spots

# Quality control 
## Select over tissue genes
in.tissue <- pDataDT(s1)[in_tissue == 1]$cell_ID

if (sum((in.tissue == pDataDT(s1)$cell_ID) == F) == 0) {
  cat("All genes are over tissue")
  rm(in.tissue)
} else {
  s1 <- subsetGiotto(s1, cell_ids = in.tissue)
  rm(in.tissue)
  } # All genes are over tissue

## Identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", s1@feat_ID$rna)
if (sum(is_mito == T) != 0) {
  s1@feat_ID$is_mito <- is_mito
  s1.mito <- subsetGiotto(s1, feat_ids = s1@feat_ID$rna[is_mito == T])
  cat(sum(is_mito == T), "mitochondrial genes found:", s1.mito@feat_ID$rna)
  rm(is_mito)
  } else {
    cat('No mitochondrial genes found')
    rm(is_mito)
    } # 13 mitochondrial genes found

### Mito genes proportions
s1.expr.sum <- apply(s1.matrix, 2, sum) # Sum of gene expression values per cell
s1.mito.sum <- apply(s1.matrix[s1.mito@feat_ID$rna, ], 2, sum) # Sum of mito gene expression values per cell

s1.mito.prop <- which(s1.mito.sum/s1.expr.sum > 0.28) # No cells exceed 0.28 mito genes proportion
rm(s1.expr.sum, s1.mito.sum)

## Plots
s1.detected <- filterDistributions(s1, detection = "cells", nr_bins = 150,
                                   method = "threshold", expression_threshold = 1, 
                                   default_save_name = "detected_genes")
s1.detected # Detected genes per cell, per sample


s1.libsize <- filterDistributions(s1, detection = "cells", nr_bins = 150,
                                  method = "sum", 
                                  default_save_name = "library_size")
s1.libsize # Library size, per sample


s1.thresholds <- filterCombinations(s1, expression_thresholds = c(1, 2, 3),
                                    feat_det_in_min_cells = c(50, 50, 100, 100),
                                    min_det_feats_per_cell = c(100, 250, 500, 1000), 
                                    show_plot = F,
                                    default_save_name = "thresholds")
s1.thresholds # Threshold evaluation, per sample


####--SAMPLE 2--####