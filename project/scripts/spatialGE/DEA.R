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

DEA <- STdiff(newsGE.obj, samples = c("MS377N", "MS549H"), annot = "domain", clusters = c("LC", "LR"),
              topgenes = 500, test_type = "mm", cores = 12)

save(DEA, file = "./project/material/spatialGE/DEA.RData")
