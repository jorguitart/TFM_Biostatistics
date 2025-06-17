####--LIBRARIES--####
library(tidyverse)
library(ggpubr)
library(Giotto) #pak::pkg_install("drieslab/Giotto")
# library(spacexr) #devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
# library(nicheDE) #devtools::install_github("kaishumason/NicheDE")
# library(spatialGE) #devtools::install_github("fridleylab/spatialGE")


####--DATA--####
sample <- loadGiotto(path_to_folder = "./project/material/filtered_samples/resolved_sample", 
                     python_path = "/usr/bin/python36")

spatPlot2D(sample, group_by = "list_ID", cell_color = "hmrf k=23 b=50.00", show_plot = T, point_shape = "no_border")


