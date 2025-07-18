####--LIBRARIES--####
library(tidyverse)
library(ggpubr)
library(Giotto) #pak::pkg_install("drieslab/Giotto")
# library(spacexr) #devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
# library(nicheDE) #devtools::install_github("kaishumason/NicheDE")
library(spatialGE) #devtools::install_github("fridleylab/spatialGE")

setwd("~/TFM")

####--DATA--####
sample <- loadGiotto("./project/material/filtered_samples/resolved_sample",
                     python_path = "C:/ProgramData/anaconda3/python.exe")

rnacounts <- as.data.frame(as.matrix(sample@expression$cell$rna$normalized@exprMat))
spotcoords <- sample@spatial_locs$cell$raw@coordinates



rna <- list(list(rnacounts)); names(rna) <- "sample"
spot <- list(list(spotcoords)); names(spot) <- "sample"

sGE.obj <- STlist(rnacounts = rna, spotcoords = spot)
