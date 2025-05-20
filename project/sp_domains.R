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


# Spatial extraction of genes
sample <- readRDS("./project/material/filtered_samples/merge_norm.rds") # Sample

## Create network
sample <- createSpatialNetwork(sample, minimum_k = 2, name = "spat_network", method = "Delaunay")

## Execute extraction
sample <- binSpect(sample, expression_values = "normalized", bin_method = "rank",
                   spatial_network_name = "spat_network", do_parallel = T, cores = 4, return_gobject = T)

## Initialize HMRF 
sample@dimension_reduction$cells$cell$rna$spatial$spatial_feat <- sample@dimension_reduction$cells$cell$rna$pca$pca
sample.hmrf <- initHMRF_V2(sample, spat_unit = "cell", feat_type = "rna", cl.method = "leiden", 
                           metadata_to_use = "leiden_clus" , spatial_network_name = "spat_network")

HMRF.model <- doHMRF_V2(sample.hmrf, betas = c(0, 5, 20))

## Timer stop
t1 <- Sys.time() - t0; t1

