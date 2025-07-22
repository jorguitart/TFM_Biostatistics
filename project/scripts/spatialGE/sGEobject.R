####--LIBRARIES--####
library(tidyverse)
library(ggpubr)
library(Giotto) #pak::pkg_install("drieslab/Giotto")
# library(spacexr) #devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
# library(nicheDE) #devtools::install_github("kaishumason/NicheDE")
library(spatialGE) #devtools::install_github("fridleylab/spatialGE")

setwd("~/TFM")

####--DATA--####
sample <- loadGiotto("./project/material/Giotto/HMRF_sample")

rnacounts <- as.data.frame(as.matrix(sample@expression$cell$rna$raw@exprMat))
counts <- rnacounts
counts <- data.frame(genes = rownames(counts), counts)
names.counts <- gsub("\\.", "-", colnames(counts)[-1])
colnames(counts) <- c("genes", names.counts)

spotcoords <- as.data.frame(sample@spatial_locs$cell$raw@coordinates)
spotcoords <- spotcoords[, c(3, 2, 1)]
spotcoords <- spotcoords[spotcoords$cell_ID == colnames(counts[, -1]), ]

rna <- list(sample = counts)
spot <- list(sample = spotcoords)

sGE.obj <- STlist(rnacounts = rna, spotcoords = spot)
save(sGE.obj, file = "./project/material/spatialGE/sGEobject.RData")
