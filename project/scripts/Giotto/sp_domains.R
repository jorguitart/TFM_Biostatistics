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
                             gene_samples = 500, gene_sampling_rate = 1, hmrf_seed = 100, k = 9,
                             spatial_network_name = "spat_network", cl.method = "km")
  
  message("Saving environment image...")
  save(sample.hmrf, file = "./project/material/preHMRF.RData")
  message("Done."); cat("\n")
} else {message("File found. Starting HMRF...")}

# Run HMRF model
sample <- loadGiotto("./project/material/filtered_samples/preinit_sample")
load("./project/material/preHMRF.RData")

HMRF.model <- doHMRF_V2(sample.hmrf, betas = c(0, 5, 4))
save(HMRF.model, file = "./project/material/HMRF.RData")

prob <- HMRF.model$`k=9 b=15.00`$prob
assigned <- apply(prob, 1, which.max)
assigned <- assigned[sample@cell_ID$cell]

sample <- addCellMetadata(sample, new_metadata = as.numeric(assigned), 
                          vector_name = "HMRF_B15")

spatPlot2D(sample, group_by = "list_ID", cell_color = "HMRF_B15", 
           point_shape = "no_border", point_size = 1.2)

cell.meta <- pDataDT(sample)
cell.meta$domain <- with(
  cell.meta,
  ifelse(HMRF_B15 %in% c(2, 4, 9) & type == "CTRL", "WM",
         ifelse(HMRF_B15 %in% c(1, 6) & type == "CTRL", "GM",
                ifelse(HMRF_B15 %in% 3 & type ==  "MSCA", "LC",
                       ifelse(HMRF_B15 %in% 6 & type == "MSCA", "GM",
                              ifelse(HMRF_B15 %in% c(1, 5) & type == "MSCA", "LR",
                                     ifelse(HMRF_B15 %in% 7 & type == "MSCA", "PPWM",
                                            ifelse(HMRF_B15 %in% 4 & type == "MSCA", "VI",
                                                   ifelse(HMRF_B15 %in% c(3, 4) & type == "MSCI", "LC",
                                                          ifelse(HMRF_B15 %in% c(1, 8) & type == "MSCI", "GM",
                                                                 ifelse(HMRF_B15 %in% c(5, 9) & type == "MSCI", "VI",
                                                                        ifelse(HMRF_B15 %in% c(2, 6) & type == "MSCI", "PPWM",
                                                                               ifelse(HMRF_B15 %in% 7 & type == "MSCI", "LR", NA)))))))))))))

cell.meta <- createCellMetaObj(cell.meta)
sample <- setCellMetadata(sample, x = cell.meta)

spatPlot2D(sample, group_by = "list_ID", cell_color = "domain", 
           point_shape = "no_border", point_size = 1.2)

saveGiotto(sample, foldername = "resolved_sample", 
           dir = "./project/material/filtered_samples", overwrite = T)

message("Done."); t1 <- Sys.time() - t0; t1