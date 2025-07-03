####--LIBRARIES--####
library(tidyverse)
library(ggpubr)
library(Giotto) #pak::pkg_install("drieslab/Giotto")
# library(spacexr) #devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
# library(nicheDE) #devtools::install_github("kaishumason/NicheDE")
# library(spatialGE) #devtools::install_github("fridleylab/spatialGE")

setwd("~/TFM")

####--HMRF--####
t0 <- Sys.time()

if(!file.exists("./project/material/preHMRF.RData")) {
  
  message("Loading sample...")
  sample <- loadGiotto("./project/material/filtered_samples/normalized_sample")
  
  message("Creating spatial network...")
  sample <- createSpatialDelaunayNetwork(sample, name = "spat_network")
  
  message("Extracting spatial genes...")
  sample <- binSpect(sample, expression_values = "normalized", bin_method = "kmeans", 
                     spatial_network_name = "spat_network", return_gobject = T)
  
  message("Creating HMRF object...")
  saveGiotto(sample, foldername = "preinit_sample", dir = "./project/material/filtered_samples", overwrite = T)
  sample.hmrf <- initHMRF_V2(sample, use_spatial_genes = "binSpect", gene_list_from_top = 500, use_pca = F,
                             gene_samples = 500, gene_sampling_rate = 1, hmrf_seed = 100, k = 13)
  
  message("Saving environment image...")
  save(sample.hmrf, file = "./project/material/preHMRF.RData")
  message("Done."); cat("\n")
} else {message("File found. Starting HMRF...")}

# Run HMRF model
sample <- loadGiotto("./project/material/filtered_samples/preinit_sample")
load("./project/material/preHMRF.RData")
HMRF.model <- doHMRF_V2(sample.hmrf, c(0, 5, 6)); save(HMRF.model, file = "./project/material/HMRF.RData")
sample <- addHMRF_V2(sample, HMRFoutput = HMRF.model)
saveGiotto(sample, foldername = "resolved_sample", dir = "./project/material/filtered_samples", overwrite = T)
message("Done."); t1 <- Sys.time() - t0; t1
