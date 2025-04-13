####--LIBRARIES--####

library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(Giotto) #pak::pkg_install("drieslab/Giotto")
# library(spacexr) #devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
# library(nicheDE) #devtools::install_github("kaishumason/NicheDE")
# library(spatialGE) #devtools::install_github("fridleylab/spatialGE")

####--SAMPLE 01--####

# Normalization
load("./project/material/filtered_samples/s01_filtered.R") # Load data
s01 <- normalizeGiotto(s01.filtered); s01; rm(s01.filtered)

## Add statistics
s01 <- addFeatStatistics(s01); s01 <- addCellStatistics(s01)


# Visualization

## Find mito genes percentatge per spot
mito_genes <- grep("MT-", s01@feat_ID, value = T)
s01 <- addFeatsPerc(s01, feats = mito_genes, vector_name = "mito"); rm(mito_genes)

## Create plot
p1 <- spatPlot2D(s01, cell_color = "nr_feats", color_as_factor = F)
p2 <- spatPlot2D(s01, cell_color = "mito", color_as_factor = F)
p01 <- ggarrange(p1, p2); rm(p1, p2); p01


####--SAMPLE 02--####

# Normalization
load("./project/material/filtered_samples/s02_filtered.R") # Load data
s02 <- normalizeGiotto(s02.filtered); s02; rm(s02.filtered)

## Add statistics
s02 <- addFeatStatistics(s02); s02 <- addCellStatistics(s02)


# Visualization

## Find mito genes percentatge per spot
mito_genes <- grep("MT-", s02@feat_ID, value = T)
s02 <- addFeatsPerc(s02, feats = mito_genes, vector_name = "mito"); rm(mito_genes)

## Create plot
p1 <- spatPlot2D(s02, cell_color = "nr_feats", color_as_factor = F)
p2 <- spatPlot2D(s02, cell_color = "mito", color_as_factor = F)
p02 <- ggarrange(p1, p2); rm(p1, p2); p02


