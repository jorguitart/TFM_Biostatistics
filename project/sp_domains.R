####--LIBRARIES--####
library(tidyverse)
library(ggpubr)
library(Giotto) #pak::pkg_install("drieslab/Giotto")
# library(spacexr) #devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
# library(nicheDE) #devtools::install_github("kaishumason/NicheDE")
# library(spatialGE) #devtools::install_github("fridleylab/spatialGE")


####--HMRF--####
t0 <- Sys.time()

if(!file.exists("./project/material/preHMRF.RData")) {
  
  message("Loading sample...")
  sample <- loadGiotto(path_to_folder = "./project/material/filtered_samples/normalized_sample", 
                       python_path = "C:/ProgramData/anaconda3/python.exe")
  
  message("Creating spatial network...")
  sample <- createSpatialNetwork(sample, minimum_k = 6, name = "spat_network", method = "Delaunay", maximum_distance_delaunay = 40) 
  
  message("Extracting spatial genes...")
  sample <- binSpect(sample, expression_values = "normalized", bin_method = "kmeans", spatial_network_name = "spat_network", return_gobject = T)
  sp.feats <- sample@feat_metadata$cell$rna$feat_ID[order(sample@feat_metadata$cell$rna$binSpect.pval)][1:500]
  
  message("Calculating spatial correlation genes...")
  sp.cor <- detectSpatialCorFeats(sample, method = "network", spatial_network_name = "spat_network", subset_feats = sp.feats)
  sp.cor <- clusterSpatialCorFeats(sp.cor, name = "spnetwork_cl", k = 13)
  cl.feats <- showSpatialCorFeats(sp.cor, use_clus_name = "spnetwork_cl", show_top_feats = 1)
  
  message("Creating HMRF object...")
  sample@dimension_reduction$cells$cell$rna$spatial$spatial_feat <- sample@dimension_reduction$cells$cell$rna$pca$pca
  saveGiotto(sample, foldername = "preinit_sample", dir = "./project/material/filtered_samples", overwrite = T)
  sample.hmrf <- initHMRF_V2(sample, spat_unit = "cell", feat_type = "rna", cl.method = "leiden", user_gene_list = sp.feats, 
                             spatial_network_name = "spat_network", k = 7, resolution.cl = 0.25)
  
  message("Saving environment image...")
  save(sample.hmrf, file = "./project/material/preHMRF.RData")
  message("Done."); cat("\n")
} else {message("File found. Starting HMRF...")}

# Run HMRF model
sample <- loadGiotto(path_to_folder = "./project/material/filtered_samples/init_sample", 
                     python_path = "/usr/bin/python36")
load("./project/material/preHMRF.RData")
HMRF.model <- doHMRF_V2(sample.hmrf); save(HMRF.model, file = "./project/material/HMRF.RData")
sample <- addHMRF_V2(sample, HMRFoutput = HMRF.model)
saveGiotto(sample, foldername = "resolved_sample", dir = "./project/material/filtered_samples", overwrite = T)
message("Done."); t1 <- Sys.time() - t0; t1
