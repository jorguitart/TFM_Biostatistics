####--LIBRARIES--####
library(tidyverse)
library(ggpubr)
library(Giotto) #pak::pkg_install("drieslab/Giotto")
# library(spacexr) #devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
# library(nicheDE) #devtools::install_github("kaishumason/NicheDE")
library(spatialGE) #devtools::install_github("fridleylab/spatialGE")

setwd("~/TFM")

dir <- "./project/material/GSE279181"
sam <- c("CO37", "CO40", "CO41", "CO74", "CO85", "CO96", 
         "MS197D", "MS197U", "MS377N", "MS377T", "MS377I", "MS411", 
         "MS497I", "MS497T", "MS549H", "MS549T")
sample.path <- file.path(dir, sam)

sGE.obj <- STlist(rnacounts = sample.path, samples = sam)
sGE.obj <- filter_data(sGE.obj, spot_mingenes = 125, gene_minspots = 28)
sGE.obj <- transform_data(sGE.obj, scale_f = 6000)
save(sGE.obj, file = "./project/material/sGEobject.RData")
