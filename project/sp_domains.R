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
  sample <- createSpatialNetwork(sample, minimum_k = 2, name = "spat_network", method = "Delaunay") 
  
  message("Extracting spatial genes...")
  sample <- binSpect(sample, expression_values = "normalized", bin_method = "rank",
                     spatial_network_name = "spat_network", do_parallel = T, cores = 4, return_gobject = T)
  feat.meta <- data.frame(fDataDT(sample)); top500 <- feat.meta[order(feat.meta$binSpect.pval), ][1:500, 1]
  
  message("Creating HMRF object...")
  sample@dimension_reduction$cells$cell$rna$spatial$spatial_feat <- sample@dimension_reduction$cells$cell$rna$pca$pca
  sample.hmrf <- initHMRF_V2(sample, spat_unit = "cell", feat_type = "rna", cl.method = "leiden", user_gene_list = top500, 
                             metadata_to_use = "leiden_clus" , spatial_network_name = "spat_network")
  
  message("HMRF object created.")
  
  message("Saving environment image...")
  saveGiotto(sample, foldername = "init_sample", dir = "./project/material/filtered_samples")
  save(sample.hmrf, file = "./project/material/preHMRF.RData")
  message("Done."); cat("\n")
} else {message("File found. Starting HMRF...")}

# Run HMRF model
sample <- loadGiotto(path_to_folder = "./project/material/filtered_samples/init_sample", 
                     python_path = "C:/ProgramData/anaconda3/python.exe")
load("./project/material/preHMRF.RData")
HMRF.model <- doHMRF_V2(sample.hmrf); save(HMRF.model, file = "./project/material/HMRF.RData")
sample <- addHMRF_V2(sample, HMRFoutput = HMRF.model)
saveGiotto(sample, foldername = "resolved_sample", dir = "./project/material/filtered_samples")
message("Done."); t1 <- Sys.time() - t0; t1

