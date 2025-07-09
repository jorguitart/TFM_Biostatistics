library(zellkonverter)
library(Seurat)
library(SeuratObject)
library(SummarizedExperiment)
library(Matrix)
library(nicheDE)
library(spacexr)
library(Giotto)
library(dplyr)

setwd("~/TFM")

# Data
sample <- loadGiotto("./project/material/filtered_samples/preinit_sample")
atlas <- readH5AD("./project/material/GSE279181/GSE279180_sn_atlas.h5ad", reader = "R")

# Create reference
expr.mat <- assay(atlas, "X")
cell.meta <- as.data.frame(colData(atlas))
celltypes <- as.factor(cell.meta$celltype)
names(celltypes) <- rownames(cell.meta)
ref <- Reference(expr.mat, celltypes)

# Create puck
counts <- sample@expression$cell$rna$raw@exprMat
coord <- sample@spatial_locs$cell$raw@coordinates
coords <- coord[, -3]; rownames(coords) <- coord$cell_ID; colnames(coords) <- c("imagerow", "imagecol")
puck <- SpatialRNA(coords, counts)

# Create RCTD object
myRCTD <- create.RCTD(puck, ref, max_cores = 10, UMI_min = 0)

# Deconvolution
CT <- unique(celltypes)
lib.siz <- rep(NA, length(CT))
cell.lib <- colSums(expr.mat)

for (j in c(1:length(CT))) {
  ct <- CT[j]
  lib.siz[j] = mean(cell.lib[celltypes == ct])
}

deconv.mat <- matrix(0, nrow(coords), length(CT))
colnames(deconv.mat) <- CT; rownames(deconv.mat) <- rownames(coords)

for (j in c(1:length(myRCTD@results))) {
  fills = match(myRCTD@results[[j]]$cell_type_list, CT)
  deconv.mat[j, fills] = myRCTD@results[[j]]$sub_weights
  deconv.mat[j, ] = deconv.mat[j, ]/sum(deconv.mat[j, ])
  spot.lib = deconv.mat[j, ]*1/lib.siz
  deconv.mat[j, ] = spot.lib/sum(spot.lib)
}

save(deconv.mat, file = "./project/material/deconv.RData")