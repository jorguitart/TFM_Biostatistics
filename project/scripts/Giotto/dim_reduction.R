####--LIBRARIES--####
library(tidyverse)
library(ggpubr)
library(Giotto) #pak::pkg_install("drieslab/Giotto")
# library(spacexr) #devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
# library(nicheDE) #devtools::install_github("kaishumason/NicheDE")
# library(spatialGE) #devtools::install_github("fridleylab/spatialGE")

setwd("~/TFM")

####--SAMPLE--####
# Timer start
t0 <- Sys.time()


# Normalization
merged.samples <- loadGiotto("./project/material/filtered_samples/merged_sample")
merged.samples <- normalizeGiotto(merged.samples, norm_methods = "standard"); merged.samples

## Add statistics
merged.samples <- addFeatStatistics(merged.samples); merged.samples <- addCellStatistics(merged.samples)


# HVGs - Dim reduction
## Calculate HVGs
merged.samples <- calculateHVF(merged.samples, expression_values = "normalized",
                               method = "cov_loess", show_plot = T)

## Run dim reduction
### PCA
merged.samples <- runPCA(merged.samples, expression_values = "normalized", feats_to_use = "hvf", 
                         ncp = 50, name = "pca")

### UMAP
merged.samples <- runUMAP(merged.samples, dimensions_to_use = 1:20, dim_reduction_name = "pca", 
                          feats_to_use = "hvf")


# Clustering
merged.samples <- runGiottoHarmony(merged.samples, dim_reduction_to_use = "pca", 
                                   vars_use = "list_ID", dimensions_to_use = 1:10)
merged.samples <- createNearestNetwork(merged.samples, spat_unit = "cell", feat_type = "rna", 
                                       feats_to_use = "hvf", dim_reduction_name = "pca", name = "sNN.pca")
merged.samples <- doLeidenCluster(merged.samples, name = "leiden_clus", resolution = 0.25)


# Markers
markers <- findMarkers_one_vs_all(merged.samples, method = "scran", expression_values = "normalized", 
                                  cluster_column = "leiden_clus")
save(markers, file = "./project/material/marker_genes.rda")

# Plot dim reduction
## PCA
merged.pca <- plotPCA(merged.samples, dim_reduction_name = "pca"); merged.pca


### Cumulative variance explained
merge.scree <- screePlot(merged.samples, expression_values = "normalized",
                         dim_reduction_name = "pca", feats_to_use = "hvf", ncp = 50)


## UMAP
merge.umap <- plotUMAP(merged.samples, cell_color = "type", point_size = 2,
                     point_shape = "no_border", label_size = 0, title = "", point_alpha = 0.5)


# Save sample
saveGiotto(merged.samples, foldername = "normalized_sample", 
           dir = "./project/material/filtered_samples", overwrite = T)

## Timer stop
t1 <- Sys.time() - t0; t1
