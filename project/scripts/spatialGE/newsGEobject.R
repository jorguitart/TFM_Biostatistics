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
newsGE.obj <- transform_data(newsGE.obj, scale_f = 6000)

# Add metadata
load("./project/material/spatialGE/sGEclust.RData")

## Normalized counts
tr_counts <- list()
for (i in unique(sample@cell_metadata$cell$rna$list_ID)) {
  tr_counts[[i]] <- sGE.obj@tr_counts$sample[, grep(i, sGE.obj@tr_counts$sample@Dimnames[[2]])]
}
newsGE.obj@tr_counts <- tr_counts

## Spatial metadata
sp.meta <- sGE.obj@spatial_meta$sample
names(sp.meta)[names(sp.meta) == "stclust_spw0.025_dsplFalse"] <- "cluster"

sp.meta$type <- with(
  sp.meta,
  ifelse(libname %in% sp.meta$libname[grep("CO37|CO40|CO41|CO74|CO85", libname)], "CTRL", 
         ifelse(libname %in% sp.meta$libname[grep("MS497T|MS549H|MS549T", libname)], "MSCI", "MSCA")))

sp.meta$domain <- with(
  sp.meta,
  ifelse(cluster %in% c(1, 2, 3, 4, 5, 6, 7, 9) & type == "CTRL", "WM",
         ifelse(cluster %in% 8 & type == "CTRL", "GM",
                ifelse(cluster %in% c(1, 3, 4) & type == "MSCA", "VI",
                       ifelse(cluster %in% 2 & type == "MSCA", "LC",
                              ifelse(cluster %in% 6 & type == "MSCA", "PPWM",
                                     ifelse(cluster %in% 5 & type == "MSCA", "LC",
                                            ifelse(cluster %in% 7 & type == "MSCA", "LR",
                                                   ifelse(cluster %in% c(8, 9) & type == "MSCA", "GM", 
                                                          ifelse(cluster %in% 9 & type == "MSCI", "GM",
                                                                 ifelse(cluster %in% c(3, 8) & type == "MSCI", "LC",
                                                                        ifelse(cluster %in% 2 & type == "MSCI", "LR",
                                                                               ifelse(cluster %in% 4 & type == "MSCI", "PPWM",
                                                                                      ifelse(cluster %in% c(1, 5, 6 & type == "MSCI", "VI", NA)))))))))))))))

### Clusters
clusters <- list()
for (i in unique(sample@cell_metadata$cell$rna$list_ID)) {
  clusters[[i]] <- sp.meta$cluster[grep(i, sp.meta$libname)]
}

for (i in unique(sample@cell_metadata$cell$rna$list_ID)) {
  newsGE.obj@spatial_meta[[i]][["cluster"]] <- clusters[[i]]
}

### Sample types
types <- list()
for (i in unique(sample@cell_metadata$cell$rna$list_ID)) {
  types[[i]] <- sp.meta$type[grep(i, sp.meta$libname)]
}

for (i in unique(sample@cell_metadata$cell$rna$list_ID)) {
  newsGE.obj@spatial_meta[[i]][["type"]] <- types[[i]]
}

### Spatial domains
domains <- list()
for (i in unique(sample@cell_metadata$cell$rna$list_ID)) {
  domains[[i]] <- sp.meta$domain[grep(i, sp.meta$libname)]
}

for (i in unique(sample@cell_metadata$cell$rna$list_ID)) {
  newsGE.obj@spatial_meta[[i]][["domain"]] <- domains[[i]]
}

## Sample type
sample_type <- c(rep("CTRL", 5), rep("MSCA", 4), rep("MSCI", 3))
newsGE.obj@sample_meta$sample_type <- sample_type

save(newsGE.obj, file = "./project/material/spatialGE/newsGEobject.RData")
