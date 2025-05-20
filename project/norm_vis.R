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
# Timer start
t0 <- Sys.time()


# Normalization
s01 <- readRDS("./project/material/filtered_samples/s01_filtered.rds") # Sample 01 (example)
s01 <- normalizeGiotto(s01); s01

merged.samples <- readRDS("./project/material/filtered_samples/merge.rds") # Merged data
merged.samples <- normalizeGiotto(merged.samples); merged.samples

## Add statistics
s01 <- addFeatStatistics(s01); s01 <- addCellStatistics(s01)
merged.samples <- addFeatStatistics(merged.samples); merged.samples <- addCellStatistics(merged.samples)


# Visualization
p1 <- spatPlot2D(s01, cell_color = "nr_feats", color_as_factor = F)
p2 <- spatPlot2D(s01, cell_color = "mito_perc", color_as_factor = F)
p01 <- ggarrange(p1, p2); rm(p1, p2); p01

ggsave("./project/outcomes/vis/s01.png", plot = p01, scale = 3, width = 1920, height = 1080, units = "px")


# HVGs - Dim reduction
## Calculate HVGs
merged.samples <- calculateHVF(merged.samples, expression_values = "normalized")

## Run dim reduction
### PCA
merged.samples <- runPCA(merged.samples, expression_values = "normalized", feats_to_use = "hvf", ncp = 50,
                         name = "pca")

### UMAP
merged.samples <- runUMAP(merged.samples, dimensions_to_use = 1:7, n_components = 2, 
                          dim_reduction_name = "pca", feats_to_use = "hvf")


# Clustering
merged.samples <- createNearestNetwork(merged.samples, spat_unit = "cell", feat_type = "rna", 
                                       dimensions_to_use = 1:7, feats_to_use = "hvf", 
                                       dim_reduction_name = "pca", name = "sNN.pca")
merged.samples <- doLeidenCluster(merged.samples, name = "leiden_clus", resolution = 0.25) # Leiden


# Plot dim reduction
## PCA
merged.pca <- plotPCA(merged.samples, dim_reduction_name = "pca"); merged.pca


### Cumulative variance explained
merge.scree <- screePlot(merged.samples, expression_values = "normalized",
                         dim_reduction_name = "pca", feats_to_use = "hvf", ncp = 50); merge.scree


## UMAP
merge.umap <- plotUMAP(merged.samples, cell_color = "km_clus", point_size = 2,
                     point_shape = "no_border", label_size = 0, title = ""); merge.umap


# Save sample
saveRDS(merged.samples, file = "./project/material/filtered_samples/merge_norm.rds")

## Timer stop
t1 <- Sys.time() - t0; t1
