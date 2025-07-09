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
atlas <- readH5AD("./project/material/GSE279181/GSE279180_sn_atlas.h5ad", reader = "R")

expr.mat <- assay(atlas, "X")
expr.mat <- t(as.matrix(expr.mat))

cell.meta <- as.data.frame(colData(atlas))
celltypes <- data.frame(colnames.refr. = rownames(cell.meta), Idents.refr. = cell.meta$celltype)
rownames(celltypes) <- celltypes$colnames.refr.

lib.mat <- CreateLibraryMatrix(expr.mat, celltypes)
save(lib.mat, file = "./project/material/library.RData")
