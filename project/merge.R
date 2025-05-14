####--LIBRARIES--####
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(Giotto) #pak::pkg_install("drieslab/Giotto")
# library(spacexr) #devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
# library(nicheDE) #devtools::install_github("kaishumason/NicheDE")
# library(spatialGE) #devtools::install_github("fridleylab/spatialGE")


####--DATA--####
s01 <- readRDS("./project/material/filtered_samples/s01_filtered.rds")
s02 <- readRDS("./project/material/filtered_samples/s02_filtered.rds")
s03 <- readRDS("./project/material/filtered_samples/s03_filtered.rds")
s04 <- readRDS("./project/material/filtered_samples/s04_filtered.rds")
s05 <- readRDS("./project/material/filtered_samples/s05_filtered.rds")
s06 <- readRDS("./project/material/filtered_samples/s06_filtered.rds")
s08 <- readRDS("./project/material/filtered_samples/s08_filtered.rds")
s09 <- readRDS("./project/material/filtered_samples/s09_filtered.rds")
s11 <- readRDS("./project/material/filtered_samples/s11_filtered.rds")
s12 <- readRDS("./project/material/filtered_samples/s12_filtered.rds")
s14 <- readRDS("./project/material/filtered_samples/s14_filtered.rds")
s15 <- readRDS("./project/material/filtered_samples/s15_filtered.rds")
s16 <- readRDS("./project/material/filtered_samples/s16_filtered.rds")
s17 <- readRDS("./project/material/filtered_samples/s17_filtered.rds")
s18 <- readRDS("./project/material/filtered_samples/s18_filtered.rds")


####--MERGED SAMPLE--####
# Obtain combined expression matrix
s.list <- list(s01, s02, s03, s04, s05, s06, s08, s09, s11, s12, s14, s15, s16, s17, s18)

matrix.list <- lapply(s.list, function(s) getExpression(s, values = "raw", output = "matrix"))
com.genes <- Reduce(intersect, lapply(matrix.list, rownames))
matrix.list.com <- lapply(seq_along(matrix.list), function(i) {
  mat <- matrix.list[[i]][com.genes, ]
  colnames(mat) <- paste0("sample", i, "_", colnames(mat))
  return(mat)
})

comb.matrix <- do.call(cbind, matrix.list.com)

# Obtain combined spatial coords
coord.list <- lapply(seq_along(s.list), function(i) {
  coords <- getSpatialLocations(s.list[[i]])
  coords$cell_ID <- paste0("sample", i, "_", coords$cell_ID)
  coords$sample <- paste0("sample", i)
  coords
})

comb.coords <- do.call(rbind, coord.list)

# Obtain combined metadata
meta.list <- lapply(seq_along(s.list), function (i) {
  meta <- pDataDT(s.list[[i]])
  meta$cell_ID <- paste0("sample", i, "_", meta$cell_ID)
  meta$sample <- paste0("sample", i)
  meta
})

comb.meta <- do.call(rbind, meta.list)

# Create and save gobject
merged.samples <- createGiottoObject(
  expression = comb.matrix, spatial_locs = comb.coords, cell_metadata = comb.meta, 
  instructions = createGiottoInstructions(python_path = "C:/ProgramData/anaconda3/python.exe"))

saveRDS(merged.samples, file = "./project/material/filtered_samples/merge.rds")

## Remove unnecessary objetcs
rm(s01, s02, s03, s04, s05, s06, s08, s09, s11, s12, s14, s15, s16, s17, s18, 
   s.list, com.genes, matrix.list, matrix.list.com, coord.list, meta.list, comb.matrix, comb.coords, comb.meta)
