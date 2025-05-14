####--LIBRARIES--####

library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(Giotto) #pak::pkg_install("drieslab/Giotto")
# library(spacexr) #devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
# library(nicheDE) #devtools::install_github("kaishumason/NicheDE")
# library(spatialGE) #devtools::install_github("fridleylab/spatialGE")

####--SAMPLE--####

# Normalization
merged.samples <- readRDS("./project/material/filtered_samples/merge.rds") # Load data
merged.samples <- normalizeGiotto(merged.samples); merged.samples

## Add statistics
merged.samples <- addFeatStatistics(merged.samples); merged.samples <- addCellStatistics(merged.samples)


# # Visualization
# p1 <- spatPlot2D(s01, cell_color = "nr_feats", color_as_factor = F)
# p2 <- spatPlot2D(s01, cell_color = "mito_perc", color_as_factor = F)
# p01 <- ggarrange(p1, p2); rm(p1, p2); p01
# 
# ggsave("./project/outcomes/vis/s01.png", plot = p01, scale = 3, width = 1920, height = 1080, units = "px")

# HVGs - Dim reduction
## Calculate HVGs
merged.samples <- calculateHVF(merged.samples, expression_values = "normalized")

## Run dim reduction
### PCA
merged.samples <- runPCA(merged.samples, expression_values = "normalized", feats_to_use = "hvf", ncp = 50)

### UMAP
merged.samples <- runUMAP(merged.samples, dimensions_to_use = 1:7, n_components = 2)


# Clustering
merged.samples <- createNearestNetwork(merged.samples, dimensions_to_use = 1:7)
merged.samples <- doLeidenCluster(merged.samples, name = "leiden_clus", resolution = 0.25)

# Plot dim reduction
## PCA
merged.pca <- plotPCA(merged.samples); merged.pca

### Cumulative variance explained
merge.scree <- screePlot(merged.samples, expression_values = "normalized",
                       feats_to_use = "hvf", ncp = 50); merge.scree

## UMAP
merge.umap <- plotUMAP(merged.samples, cell_color = "leiden_clus", point_size = 2,
                     point_shape = "no_border", label_size = 0,
                     title = "Clusters"); merge.umap
