####--LIBRARIES--####
library(tidyverse)
library(ggpubr)
# library(Giotto) #pak::pkg_install("drieslab/Giotto")
# library(spacexr) #devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
# library(nicheDE) #devtools::install_github("kaishumason/NicheDE")
library(spatialGE) #devtools::install_github("fridleylab/spatialGE")

setwd("~/TFM")

####--DATA--####
load("./project/material/spatialGE/newsGEobject.RData")

sGE.obj <- transform_data(sGE.obj, scale_f = 6000)
sGE.obj <- STclust(sGE.obj, ws = 0.5, topgenes = 500, cores = 10)
save(sGE.obj, file = "./project/material/spatialGE/sGEclust.RData")