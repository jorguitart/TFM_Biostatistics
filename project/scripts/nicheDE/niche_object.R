####--LIBRARIES--####
library(tidyverse)
library(ggpubr)
library(Giotto) #pak::pkg_install("drieslab/Giotto")
library(nicheDE) #devtools::install_github("kaishumason/NicheDE")

setwd("~/TFM")

####--SAMPLE--####
sample <- loadGiotto("./project/material/Giotto/HMRF_sample")
load("./project/material/nicheDE/library.RData")
load("./project/material/nicheDE/DWLS.rds")

counts <- t(sample@expression$cell$rna$raw@exprMat)
counts <- as.matrix(counts)
coord <- sample@spatial_locs$cell$raw@coordinates
coords <- coord[, -3]; rownames(coords) <- coord$cell_ID

int <- intersect(rownames(counts), DWLS$cell_ID)
DWLS <- DWLS[cell_ID %in% int]; deconv.mat <- as.matrix(DWLS[, -1])
rownames(deconv.mat) <- DWLS$cell_ID; deconv.mat <- deconv.mat[rownames(counts), ]


niche.obj <- CreateNicheDEObject(counts, coords, lib.mat, deconv.mat, sigma = c(1, 100, 400))

niche.obj <- CalculateEffectiveNicheLargeScale(niche.obj)

save(niche.obj, file = "./project/material/nicheDE/nicheObj.RData")

message("Starting niche-DE calculation...")
niche.obj <- niche_DE(niche.obj, num_cores = 64, outfile = "./project/material/nicheDE_track2.out")
save(niche.obj, file = "./project/material/nicheDE/nicheObj.RData")
