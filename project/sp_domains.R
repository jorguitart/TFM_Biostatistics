####--LIBRARIES--####
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(Giotto) #pak::pkg_install("drieslab/Giotto")
# library(spacexr) #devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
# library(nicheDE) #devtools::install_github("kaishumason/NicheDE")
# library(spatialGE) #devtools::install_github("fridleylab/spatialGE")


####--HMRF--####
if(!file.exists("./project/material/enrichment.RData")) {
  t0 <- Sys.time()
  
  cat("Loading sample...")
  sample <- readRDS("./project/material/filtered_samples/merge_norm.rds")
  
  cat("Creating spatial network")
  sample <- createSpatialNetwork(sample, minimum_k = 2, name = "spat_network", method = "Delaunay") 
  
  cat("Extracting spatial genes...")
  sample <- binSpect(sample, expression_values = "normalized", bin_method = "rank",
                     spatial_network_name = "spat_network", do_parallel = T, cores = 4, return_gobject = T)
  
  cat("Creating HMRF object...")
  sample@dimension_reduction$cells$cell$rna$spatial$spatial_feat <- sample@dimension_reduction$cells$cell$rna$pca$pca
  sample.hmrf <- initHMRF_V2(sample, spat_unit = "cell", feat_type = "rna", cl.method = "leiden",
                             metadata_to_use = "leiden_clus" , spatial_network_name = "spat_network")
  
  t1 <- Sys.time() - t0
  cat("HMRF object obtained. "); t1
  
  cat("Saving environment image...")
  rm(t0, t1); save.image(file = "./project/material/enrichment.RData")
  cat("Done.")
}

# Run HMRF model
load("./project/material/enrichment.RData"); t0 <- Sys.time()
HMRF.model <- doHMRF_V2(sample.hmrf, betas = c(0, 5, 20))
t1 <- Sys.time() - t0; t1

