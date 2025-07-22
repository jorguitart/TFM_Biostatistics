####--LIBRARIES--####
library(tidyverse)
library(ggpubr)
library(Giotto) #pak::pkg_install("drieslab/Giotto")
# library(spacexr) #devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
# library(nicheDE) #devtools::install_github("kaishumason/NicheDE")
# library(spatialGE) #devtools::install_github("fridleylab/spatialGE")

setwd("~/TFM")

####--ANALYSIS--####
# Timer start
t0 <- Sys.time()

# Normalization
message("Loading sample...")
sample <- loadGiotto("./project/material/Giotto/merged_sample",
                     python_path = "/usr/bin/python36")

message("Normalizing sample...")
sample <- normalizeGiotto(sample, norm_methods = "standard")
sample <- addFeatStatistics(sample); sample <- addCellStatistics(sample)

# HVGs - Dim reduction
## Calculate HVGs
message("Calculating HVF...")
sample <- calculateHVF(sample, expression_values = "normalized",
                       method = "cov_loess", show_plot = T)

## Run dim reduction
message("Reducing dimensionality...")
### PCA
sample <- runPCA(sample, expression_values = "normalized", feats_to_use = "hvf", 
                 ncp = 50, name = "pca")

### UMAP
sample <- runUMAP(sample, dimensions_to_use = 1:20, dim_reduction_name = "pca", 
                  feats_to_use = "hvf")


# Clustering
message("Running Harmony integration...")
sample <- runGiottoHarmony(sample, dim_reduction_to_use = "pca", 
                           vars_use = "list_ID", dimensions_to_use = 1:10)

message("Leiden clustering...")
sample <- createNearestNetwork(sample, spat_unit = "cell", feat_type = "rna", 
                               feats_to_use = "hvf", dim_reduction_name = "pca", name = "sNN.pca")
sample <- doLeidenCluster(sample, name = "leiden_clus", resolution = 0.25)

# Spatial domains
## binSpect
message("Creating spatial network...")
sample <- createSpatialDelaunayNetwork(sample, name = "spat_network")

message("Extracting spatial genes...")
sample <- binSpect(sample, expression_values = "normalized", bin_method = "kmeans", 
                   spatial_network_name = "spat_network", return_gobject = T)

## HMRF
message("Creating HMRF object...")
sample.hmrf <- initHMRF_V2(sample, use_spatial_genes = "binSpect", gene_list_from_top = 500, use_pca = F,
                           gene_samples = 500, gene_sampling_rate = 1, hmrf_seed = 100, k = 9,
                           spatial_network_name = "spat_network", cl.method = "km")

### Run HMRF model
message("Running HMRF...")
HMRF.model <- doHMRF_V2(sample.hmrf, betas = c(0, 5, 5))

prob <- HMRF.model$`k=9 b=20.00`$prob
assigned <- apply(prob, 1, which.max)
assigned <- assigned[sample@cell_ID$cell]

sample <- addCellMetadata(sample, new_metadata = as.numeric(assigned), 
                          vector_name = "HMRF_B20")

cell.meta <- pDataDT(sample)
cell.meta$domain <- with(
  cell.meta,
  ifelse(HMRF_B20 %in% c(3, 5, 7, 8) & type == "CTRL", "WM",
         ifelse(HMRF_B20 %in% 9 & type == "CTRL", "GM",
                ifelse(HMRF_B20 %in% 4 & type == "MSCA", "LR",
                       ifelse(HMRF_B20 %in% 6 & type == "MSCA", "GM",
                              ifelse(HMRF_B20 %in% 5 & type == "MSCA", "PPWM",
                                     ifelse(HMRF_B20 %in% c(1, 2) & type == "MSCA", "LC",
                                            ifelse(HMRF_B20 %in% c(8, 9) & type == "MSCA", "VI",
                                                   ifelse(HMRF_B20 %in% c(1, 8) & type == "MSCI", "LC",
                                                          ifelse(HMRF_B20 %in% 6 & type == "MSCI", "GM",
                                                                 ifelse(HMRF_B20 %in% 5 & type == "MSCI", "LR",
                                                                        ifelse(HMRF_B20 %in% 2 & type == "MSCI", "GM",
                                                                               ifelse(HMRF_B20 %in% 7 & type == "MSCI", "VI",
                                                                                      ifelse(HMRF_B20 %in% c(3, 9) & type == "MSCI", "PPWM", NA))))))))))))))


cell.meta <- createCellMetaObj(cell.meta)
sample <- setCellMetadata(sample, x = cell.meta)

message("Done. Time: "); t1 <- Sys.time() - t0; t1

message("Saving gobject...")
save(HMRF.model, file = "./project/material/Giotto/HMRF.RData")
saveGiotto(sample, foldername = "HMRF_sample", 
           dir = "./project/material/Giotto", overwrite = T)

message("Done.")
