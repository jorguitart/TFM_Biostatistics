####----LOAD LIBRARIES----####

library(dplyr)
library(tidyverse)
library(ggplot2)
library(spacexr) #devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
library(nicheDE) #devtools::install_github("kaishumason/NicheDE")
library(spatialGE) #devtools::install_github("fridleylab/spatialGE")
library(Giotto) #pak::pkg_install("drieslab/Giotto")
library(scran)

####--DIRECTORIES--####

# Create save instructions
instr <- createGiottoInstructions(save_dir = "./project/outcomes",
                                  save_plot = TRUE)


# Create sample dirs
dir <- "./project/material/GSE165098"
sample.ids <- c("LPS1", "LPS2", "LPS3", "saline1", "saline2", "saline3")
sample.path <- file.path(dir, sample.ids)

####--SAMPLE 1--####

# Obtain data
s1 <- createGiottoVisiumObject(visium_dir = sample.path[1], # LPS1
                               expr_data = "filter",
                               png_name = "tissue_lowres_image.png") # Not well aligned

# Visualize image
spatPlot(s1, point_alpha = 0.6)

# Quality control 
## Select over tissue genes

if (sum((in.tissue == pDataDT(s1)$cell_ID) == F) == 0) {
  cat("All genes are over tissue")
} else {
  s1 <- subsetGiotto(s1, cell_ids = in.tissue)
  }

## Identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", s1@feat_ID)
if (sum(is_mito == T) != 0) {
  cat(sum(is_mito == T), "mitochondrial genes found:", which(is_mito == T))
  } else {
    cat('No mitochondrial genes found')
    } # 13 mitochondrial genes

## Histograms
filterDistributions(s1, detection = "genes", nr_bins = 100) # Genes detected per spot (threshold = 0)


## Tresholds
filterCombinations(s1)
