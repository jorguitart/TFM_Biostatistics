####----LOAD LIBRARIES----####

library(dplyr)
library(tidyverse)
library(ggplot2)
library(spacexr) #devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
library(nicheDE) #devtools::install_github("kaishumason/NicheDE")
library(spatialGE) #devtools::install_github("fridleylab/spatialGE")
library(Giotto) #remotes::install_github("drieslab/Giotto@main")

####--DIRECTORIES--####

# Create save instructions
instr <- createGiottoInstructions(save_dir = "./project/outcomes",
                                  save_plot = TRUE)

# Create sample dirs
dir <- "./project/material/GSE165098"
sample.ids <- c("LPS1", "LPS2", "LPS3", "saline1", "saline2", "saline3")
h5.path <- file.path(dir, sample.ids, "filtered_feature_bc_matrix.h5")
pos.path <- file.path(dir, sample.ids, "spatial/tissue_positions_list.csv")
png.path <- file.path(dir, sample.ids, "spatial/tissue_lowres_image.png")


####--SAMPLE 1--####

# Obtain data
s1 <- createGiottoVisiumObject(h5_visium_path = h5.path[1], 
                               h5_tissue_positions_path = pos.path[1], 
                               h5_image_png_path = png.path[1], 
                               instructions = instr) #LPS1

# Visualize image
s1 <- updateGiottoImage(s1, image_name = "image", # Aligning
                        xmax_adj = 690, xmin_adj = 690,
                        ymax_adj = 820, ymin_adj = 520)

spatPlot(s1, show_image = T, point_alpha = 0.6)

# Quality control 
## Select over tissue genes
in.tissue <- pDataDT(s1)[in_tissue == 1]$cell_ID
s1 <- subsetGiotto(s1, cell_ids = in.tissue) 
s1@parameters$"0_subset" # No genes/cells deleted

## Identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", s1@gene_ID)
if (sum(is_mito == T) != 0) {
  cat(sum(is_mito == T), "mitochondrial genes found:", which(is_mito == T))
  } else {
    cat('No mitochondrial genes found')
    } # 13 mitochondrial genes

## Histograms
filterDistributions(s1, detection = "cells")
filterDistributions(s1, detection = "genes")

## Tresholds
filterCombinations(s1)
