####--LIBRARIES--####
library(tidyverse)
library(ggpubr)
library(Giotto) #pak::pkg_install("drieslab/Giotto")
# library(spacexr) #devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
# library(nicheDE) #devtools::install_github("kaishumason/NicheDE")
# library(spatialGE) #devtools::install_github("fridleylab/spatialGE")


####--DATA--####
sample <- loadGiotto(path_to_folder = "./project/material/filtered_samples/resolved_sample", 
                     python_path = "C:/ProgramData/anaconda3/python.exe")
sample@spatial_info <- sample@spatial_locs

spatPlot2D(sample, group_by = "list_ID", cell_color = "leiden_clus", show_plot = T)


