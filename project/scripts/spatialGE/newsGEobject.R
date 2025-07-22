####--LIBRARIES--####
library(tidyverse)
library(ggpubr)
library(Giotto) #pak::pkg_install("drieslab/Giotto")
# library(spacexr) #devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
# library(nicheDE) #devtools::install_github("kaishumason/NicheDE")
library(spatialGE) #devtools::install_github("fridleylab/spatialGE")

setwd("~/TFM")

####--DATA--####
# Create STlist
sample <- loadGiotto("./project/material/Giotto/HMRF_sample")

samples <- list(); rnacounts <- list(); names.counts <- list(); spotcoords <- list()

for (i in unique(sample@cell_metadata$cell$rna$list_ID)) {
  samples[[i]] <- subsetGiotto(sample, cell_ids = pDataDT(sample)[list_ID %in% i]$cell_ID)

  rnacounts[[i]] <- as.data.frame(as.matrix(samples[[i]]@expression$cell$rna$raw@exprMat))
  rnacounts[[i]] <- data.frame(genes = rownames(rnacounts[[i]]), rnacounts[[i]])
  names.counts[[i]] <- gsub("\\.", "-", colnames(rnacounts[[i]])[-1])
  colnames(rnacounts[[i]]) <- c("genes", names.counts[[i]])

  spotcoords[[i]] <- as.data.frame(samples[[i]]@spatial_locs$cell$raw@coordinates)
  spotcoords[[i]] <- spotcoords[[i]][, c(3, 2, 1)]
  spotcoords[[i]] <- spotcoords[[i]][spotcoords[[i]]$cell_ID == colnames(rnacounts[[i]][, -1]), ]
}

newsGE.obj <- STlist(rnacounts = rnacounts, spotcoords = spotcoords)

# Add metadata
load("./project/material/spatialGE/sGEclust.RData")

## Normalized counts
tr_counts <- list()
for (i in unique(sample@cell_metadata$cell$rna$list_ID)) {
  tr_counts[[i]] <- sGE.obj@tr_counts$sample[, grep(i, sGE.obj@tr_counts$sample@Dimnames[[2]])]
}
newsGE.obj@tr_counts <- tr_counts

## Domains
sp.meta <- sGE.obj@spatial_meta$sample
names(sp.meta)[names(sp.meta) == "stclust_spw0.025_dsplFalse"] <- "cluster"

clusters <- list()
for (i in unique(sample@cell_metadata$cell$rna$list_ID)) {
  clusters[[i]] <- sp.meta$cluster[grep(i, sp.meta$libname)]
}

for (i in unique(sample@cell_metadata$cell$rna$list_ID)) {
  newsGE.obj@spatial_meta[[i]][["cluster"]] <- clusters[[i]]
}

## Sample type
sample_type <- c(rep("CTRL", 5), rep("MSCA", 4), rep("MSCI", 3))
newsGE.obj@sample_meta$sample_type <- sample_type

save(newsGE.obj, file = "./project/material/spatialGE/newsGEobject.RData")
