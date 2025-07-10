####--LIBRARIES--####
library(tidyverse)
library(ggpubr)
library(Giotto) #pak::pkg_install("drieslab/Giotto")
# library(spacexr) #devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
library(nicheDE) #devtools::install_github("kaishumason/NicheDE")
# library(spatialGE) #devtools::install_github("fridleylab/spatialGE")

setwd("~/TFM")

####--SAMPLE--####
sample <- loadGiotto("./project/material/filtered_samples/preinit_sample")
load("./project/material/library.RData")
load("./project/material/DWLS.rds")

counts <- t(sample@expression$cell$rna$raw@exprMat)
counts <- as.matrix(counts)
coord <- sample@spatial_locs$cell$raw@coordinates
coords <- coord[, -3]; rownames(coords) <- coord$cell_ID

int <- intersect(rownames(counts), DWLS$cell_ID)
DWLS <- DWLS[cell_ID %in% int]; deconv.mat <- as.matrix(DWLS[, -1])
rownames(deconv.mat) <- DWLS$cell_ID; deconv.mat <- deconv.mat[rownames(counts), ]


niche.obj <- CreateNicheDEObject(counts, coords, lib.mat, deconv.mat, sigma = c(1, 400, 1000))
