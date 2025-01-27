####--LIBRARIES--####

library(dplyr)
library(tidyverse)
library(ggplot2)
library(spacexr) #devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
library(nicheDE) #devtools::install_github("kaishumason/NicheDE")
library(spatialGE) #devtools::install_github("fridleylab/spatialGE")
library(Giotto) #pak::pkg_install("drieslab/Giotto")

####--DIRECTORIES--####

# Create save instructions
instr <- createGiottoInstructions(save_dir = "./project/outcomes",
                                  save_plot = T, show_plot = F)


# Create sample dirs
dir <- "./project/material/GSE165098"
sample.ids <- c("LPS1", "LPS2", "LPS3", "saline1", "saline2", "saline3")
sample.path <- file.path(dir, sample.ids)

####--SAMPLE 1--####

# Obtain data
s1 <- createGiottoVisiumObject(visium_dir = sample.path[1], #LPS1
                               gene_column_index = 2, # Use gene symbols
                               expr_data = "filter", # Use filtered data
                               png_name = "tissue_lowres_image.png", # Lowres
                               instructions = instr)

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
  cat(sum(is_mito == T), "mitochondrial genes found:", which(is_mito == T))
  rm(is_mito)
  } else {
    cat('No mitochondrial genes found')
    rm(is_mito)
    } # 13 mitochondrial genes found

## Histograms
### Detected genes per cell
s1.detected <- filterDistributions(s1, detection = "cells", nr_bins = 150,
                                   method = "threshold",
                                   default_save_name = "detected_genes")
s1.detected

### Library size
s1.libsize <- filterDistributions(s1, detection = "cells", nr_bins = 150,
                                  method = "sum",
                                  default_save_name = "library_size")
s1.libsize

## Tresholds
s1.thresholds <- filterCombinations(s1, show_plot = F, 
                                    default_save_name = "thresholds")
s1.thresholds