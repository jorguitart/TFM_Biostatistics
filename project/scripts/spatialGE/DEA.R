####--LIBRARIES--####
library(tidyverse)
library(ggpubr)
# library(Giotto) #pak::pkg_install("drieslab/Giotto")
# library(spacexr) #devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
# library(nicheDE) #devtools::install_github("kaishumason/NicheDE")
library(spatialGE) #devtools::install_github("fridleylab/spatialGE")

setwd("~/TFM")

####--DATA--####
load("./project/material/spatialGE/typesGEobject.RData")

DEA <- STdiff(typesGE.obj, samples = "MSCA", annot = "domain", clusters = c("LC", "PPWM"),
              topgenes = 500, test_type = "t_test", cores = 12)

save(DEA, file = "./project/material/spatialGE/DEA.RData")
