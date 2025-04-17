####--LIBRARIES--####

library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(Giotto) #pak::pkg_install("drieslab/Giotto")
# library(spacexr) #devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
# library(nicheDE) #devtools::install_github("kaishumason/NicheDE")
# library(spatialGE) #devtools::install_github("fridleylab/spatialGE")

####--SAMPLE 01--####

# Normalization
load("./project/material/filtered_samples/s01_filtered.R") # Load data
s01 <- normalizeGiotto(s01.filtered); s01; rm(s01.filtered)

## Add statistics
s01 <- addFeatStatistics(s01); s01 <- addCellStatistics(s01)


# Visualization

## Find mito genes percentatge per spot
mito_genes <- grep("MT-", s01@feat_ID, value = T)
s01 <- addFeatsPerc(s01, feats = mito_genes, vector_name = "mito"); rm(mito_genes)

## Create plot
p1 <- spatPlot2D(s01, cell_color = "nr_feats", color_as_factor = F)
p2 <- spatPlot2D(s01, cell_color = "mito", color_as_factor = F)
p01 <- ggarrange(p1, p2); rm(p1, p2); p01


# HVGs - Dim reduction
## Calculate HVGs
s01 <- calculateHVF(s01, expression_values = "normalized")

## Run PCA
s01 <- runPCA(s01, expression_values = "normalized", feats_to_use = "hvf")
s01.pca <- plotPCA(s01); s01.pca

s01.scree <- screePlot(s01, expression_values = "normalized", 
                       feats_to_use = "hvf", ncp = 50); s01.scree

# Run UMAP
s01 <- runUMAP(s01, dimensions_to_use = 1:10, n_components = 2)
s01.umap <- plotUMAP(s01); s01.umap

# Clustering
s01 <- createNearestNetwork(s01, dimensions_to_use = 1:10)
s01 <- doLeidenCluster(s01, name = "leiden_clus")
