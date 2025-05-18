####--LIBRARIES--####
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(Giotto) #pak::pkg_install("drieslab/Giotto")
# library(spacexr) #devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
# library(nicheDE) #devtools::install_github("kaishumason/NicheDE")
# library(spatialGE) #devtools::install_github("fridleylab/spatialGE")
source("./project/initHMRF_debugged.R") # Debugged initHMRF_V2() function from Giotto


####--SAMPLE--####
# Timer start
t0 <- Sys.time()


# Spatial extraction of genes
sample <- readRDS("./project/material/filtered_samples/merge_norm.rds") # Sample

## Create network
sample <- createSpatialNetwork(sample, minimum_k = 2, name = "spat_network", method = "Delaunay")

## Execute extraction
sample <- binSpect(sample, expression_values = "normalized", bin_method = "kmeans",
                   spatial_network_name = "spat_network", kmeans_algo = "kmeans", 
                   do_parallel = T, cores = 4, return_gobject = T)

# NEEDS DEBUGGING
# sample.hmrf <- initHMRF_V2(sample, expression_values = "normalized", filter_method = "none",
#                                  spatial_network_name = "spat_network", use_spatial_genes = "binSpect", 
#                                  cl.method = "leiden", use_pca = T, use_pca_dim = 1:7,
#                                  spatial_network_name_for_neighborhood = "sNN.pca")

## Timer stop
t1 <- Sys.time() - t0; t1

