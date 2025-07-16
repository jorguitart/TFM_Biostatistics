####--LIBRARIES--####
library(tidyverse)
library(ggpubr)
library(Giotto) #pak::pkg_install("drieslab/Giotto")
# library(spacexr) #devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
# library(nicheDE) #devtools::install_github("kaishumason/NicheDE")
# library(spatialGE) #devtools::install_github("fridleylab/spatialGE")

setwd("~/TFM")

####--DATA--####
sample <- loadGiotto(path_to_folder = "./project/material/filtered_samples/resolved_sample",
                     python_path = "C:/ProgramData/anaconda3/python.exe")

spatPlot2D(sample, group_by = "list_ID", cell_color = "HMRF_B15", point_shape = "no_border")

