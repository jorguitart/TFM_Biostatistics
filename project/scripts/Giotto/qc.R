####--LIBRARIES--####
library(tidyverse)
library(ggpubr)
library(Giotto) #pak::pkg_install("drieslab/Giotto")
# library(spacexr) #devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
# library(nicheDE) #devtools::install_github("kaishumason/NicheDE")
# library(spatialGE) #devtools::install_github("fridleylab/spatialGE")

setwd("~/TFM")

####--DIRs & INSTRUCTIONS--####
#Directories
dir <- "./project/material/GSE279181"
sam <- c("CO37", "CO40", "CO41", "CO74", "CO85", "CO96", 
         "MS94", "MS197D", "MS197U", "MS229", "MS377N", "MS377T", "MS377I", "MS411", 
         "MS497I", "MS497T", "MS549H", "MS549T")
sample.path <- file.path(dir, sam); rm(dir, sam)

# Giotto instructions
instr <- createGiottoInstructions(save_plot = F)


####--SAMPLE 01--####
# Obtain data
s01 <- createGiottoVisiumObject(
  visium_dir = sample.path[1], #ct01
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = instr)
s01 <- addCellMetadata(s01, new_metadata = rep("CTRL", ncol(s01)), vector_name = "type")

# Quality control 
## Filter over tissue genes
in.tissue <- pDataDT(s01)[in_tissue == 1]$cell_ID 
if (sum((in.tissue == pDataDT(s01)$cell_ID) == F) == 0) {
  cat("All genes are over tissue")
  rm(in.tissue)
} else {
  s01 <- subsetGiotto(s01, cell_ids = in.tissue)
  rm(in.tissue)
} # All genes are over tissue

## Mito percentage
s01.matrix <- getExpression(s01, values = "raw", output = "matrix") # Extract expression matrix
is_mito <- grep("(^MT-)|(^Mt-)|(^mt-)", rownames(s01.matrix), value = T) # Identify mito genes
mito.sum <- apply(s01.matrix[is_mito, , drop = F], 2, sum)
expr.sum <- apply(s01.matrix, 2, sum)
mito.perc <- mito.sum / expr.sum * 100

s01 <- addCellMetadata(s01, new_metadata = as.numeric(mito.perc), vector_name = "mito_perc")
s01.meta <- pDataDT(s01); rm(is_mito, mito.sum, expr.sum, mito.perc)


## Create cell metadata
reads.depth <- apply(s01.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s01.matrix != 0, 2, sum) # Detected genes per cell
s01.meta$n_reads <- as.numeric(reads.depth); s01.meta$n_genes <- as.numeric(genes.per.cell)
rm(reads.depth, genes.per.cell)

### Libsize vs Detected genes correlation
s01.cor <- as.numeric(round(with(s01.meta, cor.test(n_reads, n_genes))$estimate, 4)) 


## Plots
### Detected genes
s01.detected <- filterDistributions(s01, detection = "cells", nr_bins = 50,
                                    method = "threshold", expression_threshold = 1, 
                                    default_save_name = "detected_genes")

### Genes detection along cells
s01.det.along <- filterDistributions(s01, detection = "feats", nr_bins = 50,
                                     method = "threshold", 
                                     default_save_name = "library_size")

### Libsize vs Detected genes
s01.libvdet <- ggplot(s01.meta, aes(n_reads, n_genes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = paste0("Corr = ", s01.cor), x = 11000, y = 0); rm(s01.cor)

### Threshold evaluation
s01.thresholds <- filterCombinations(s01, expression_thresholds = 1,
                                     feat_det_in_min_cells = rep(round(ncol(s01.matrix) / 20, 0), 4), # 5% 0f cells
                                     min_det_feats_per_cell = c(125, 250, 500, 1000), 
                                     show_plot = F,
                                     default_save_name = "thresholds")

## Filter sample
s01.filtered <- filterGiotto(s01, expression_values = "raw", expression_threshold = 1,
                            feat_det_in_min_cells = round(ncol(s01.matrix) / 20, 0),
                            min_det_feats_per_cell = 1000)

### Visualization plot
s01.spots.plot <- spatPlot2D(s01, cell_color = ("lightgrey"), point_size = 2,
                            select_cells = s01.filtered@cell_ID$cell, # Kept spots
                            other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                            title = "")

s01.metric.plots <- ggarrange(s01.detected, s01.det.along, 
                              s01.libvdet, s01.spots.plot) # Combine all plots
rm(s01.detected, s01.det.along, s01.libvdet, s01.thresholds); s01.metric.plots

ggsave("./project/outcomes/qc/ct01_metrics.jpg", plot = s01.metric.plots, 
       scale = 2.5, width = 1920, height = 1080, units = "px")

####--SAMPLE 02--####
# Obtain data
s02 <- createGiottoVisiumObject(
  visium_dir = sample.path[2], #ct02
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = instr)
s02 <- addCellMetadata(s02, new_metadata = rep("CTRL", ncol(s02)), vector_name = "type")

# Quality control 
## Filter over tissue genes
in.tissue <- pDataDT(s02)[in_tissue == 1]$cell_ID 
if (sum((in.tissue == pDataDT(s02)$cell_ID) == F) == 0) {
  cat("All genes are over tissue")
  rm(in.tissue)
} else {
  s02 <- subsetGiotto(s02, cell_ids = in.tissue)
  rm(in.tissue)
} # All genes are over tissue

## Mito percentage
s02.matrix <- getExpression(s02, values = "raw", output = "matrix") # Extract expression matrix
is_mito <- grep("(^MT-)|(^Mt-)|(^mt-)", rownames(s02.matrix), value = T) # Identify mito genes
mito.sum <- apply(s02.matrix[is_mito, , drop = F], 2, sum)
expr.sum <- apply(s02.matrix, 2, sum)
mito.perc <- mito.sum / expr.sum * 100

s02 <- addCellMetadata(s02, new_metadata = as.numeric(mito.perc), vector_name = "mito_perc")
s02.meta <- pDataDT(s02); rm(is_mito, mito.sum, expr.sum, mito.perc)


## Create cell metadata
reads.depth <- apply(s02.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s02.matrix != 0, 2, sum) # Detected genes per cell
s02.meta$n_reads <- as.numeric(reads.depth); s02.meta$n_genes <- as.numeric(genes.per.cell)
rm(reads.depth, genes.per.cell) 

### Libsize vs Detected genes correlation
s02.cor <- as.numeric(round(with(s02.meta, cor.test(n_reads, n_genes))$estimate, 4)) 


## Plots
### Detected genes
s02.detected <- filterDistributions(s02, detection = "cells", nr_bins = 50,
                                    method = "threshold", expression_threshold = 1, 
                                    default_save_name = "detected_genes")

### Genes detection along cells
s02.det.along <- filterDistributions(s02, detection = "feats", nr_bins = 50,
                                     method = "threshold", 
                                     default_save_name = "library_size")

### Libsize vs Detected genes
s02.libvdet <- ggplot(s02.meta, aes(n_reads, n_genes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = paste0("Corr = ", s02.cor), x = 20000, y = 0); rm(s02.cor)

### Threshold evaluation
s02.thresholds <- filterCombinations(s02, expression_thresholds = 1,
                                     feat_det_in_min_cells = rep(round(ncol(s02.matrix) / 20, 0), 4),
                                     min_det_feats_per_cell = c(125, 250, 500, 1000), 
                                     show_plot = F,
                                     default_save_name = "thresholds")

## Filter sample
s02.filtered <- filterGiotto(s02, expression_values = "raw", expression_threshold = 1,
                             feat_det_in_min_cells = round(ncol(s02.matrix) / 20, 0),
                             min_det_feats_per_cell = 750)

### Visualization plot
s02.spots.plot <- spatPlot2D(s02, cell_color = ("lightgrey"), point_size = 2,
                             select_cells = s02.filtered@cell_ID$cell, # Kept spots
                             other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                             title = "")

s02.metric.plots <- ggarrange(s02.detected, s02.det.along, 
                              s02.libvdet, s02.spots.plot) # Combine all plots
rm(s02.detected, s02.det.along, s02.libvdet, s02.thresholds); s02.metric.plots

ggsave("./project/outcomes/qc/ct02_metrics.jpg", plot = s02.metric.plots, 
       scale = 2.5, width = 1920, height = 1080, units = "px")

####--SAMPLE 03--####
# Obtain data
s03 <- createGiottoVisiumObject(
  visium_dir = sample.path[3], #ct03
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = instr)
s03 <- addCellMetadata(s03, new_metadata = rep("CTRL", ncol(s03)), vector_name = "type")

# Quality control 
## Filter over tissue genes
in.tissue <- pDataDT(s03)[in_tissue == 1]$cell_ID 
if (sum((in.tissue == pDataDT(s03)$cell_ID) == F) == 0) {
  cat("All genes are over tissue")
  rm(in.tissue)
} else {
  s03 <- subsetGiotto(s03, cell_ids = in.tissue)
  rm(in.tissue)
} # All genes are over tissue

## Mito percentage
s03.matrix <- getExpression(s03, values = "raw", output = "matrix") # Extract expression matrix
is_mito <- grep("(^MT-)|(^Mt-)|(^mt-)", rownames(s03.matrix), value = T) # Identify mito genes
mito.sum <- apply(s03.matrix[is_mito, , drop = F], 2, sum)
expr.sum <- apply(s03.matrix, 2, sum)
mito.perc <- mito.sum / expr.sum * 100

s03 <- addCellMetadata(s03, new_metadata = as.numeric(mito.perc), vector_name = "mito_perc")
s03.meta <- pDataDT(s03); rm(is_mito, mito.sum, expr.sum, mito.perc)


## Create cell metadata
reads.depth <- apply(s03.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s03.matrix != 0, 2, sum) # Detected genes per cell
s03.meta$n_reads <- as.numeric(reads.depth); s03.meta$n_genes <- as.numeric(genes.per.cell)
rm(reads.depth, genes.per.cell) 

### Libsize vs Detected genes correlation
s03.cor <- as.numeric(round(with(s03.meta, cor.test(n_reads, n_genes))$estimate, 4)) 


## Plots
### Detected genes
s03.detected <- filterDistributions(s03, detection = "cells", nr_bins = 50,
                                    method = "threshold", expression_threshold = 1, 
                                    default_save_name = "detected_genes")

### Genes detection along cells
s03.det.along <- filterDistributions(s03, detection = "feats", nr_bins = 50,
                                     method = "threshold", 
                                     default_save_name = "library_size")

### Libsize vs Detected genes
s03.libvdet <- ggplot(s03.meta, aes(n_reads, n_genes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = paste0("Corr = ", s03.cor), x = 18000, y = 0); rm(s03.cor)

### Threshold evaluation
s03.thresholds <- filterCombinations(s03, expression_thresholds = 1,
                                     feat_det_in_min_cells = rep(round(ncol(s03.matrix) / 20, 0), 4),
                                     min_det_feats_per_cell = c(125, 250, 500, 1000), 
                                     show_plot = F,
                                     default_save_name = "thresholds")

## Filter sample
s03.filtered <- filterGiotto(s03, expression_values = "raw", expression_threshold = 1,
                             feat_det_in_min_cells = round(ncol(s03.matrix) / 20, 0),
                             min_det_feats_per_cell = 500)

### Visualization plot
s03.spots.plot <- spatPlot2D(s03, cell_color = ("lightgrey"), point_size = 2,
                             select_cells = s03.filtered@cell_ID$cell, # Kept spots
                             other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                             title = "")

s03.metric.plots <- ggarrange(s03.detected, s03.det.along, 
                              s03.libvdet, s03.spots.plot) # Combine all plots
rm(s03.detected, s03.det.along, s03.libvdet, s03.thresholds); s03.metric.plots

ggsave("./project/outcomes/qc/ct03_metrics.jpg", plot = s03.metric.plots, 
       scale = 2.5, width = 1920, height = 1080, units = "px")

####--SAMPLE 04--####
# Obtain data
s04 <- createGiottoVisiumObject(
  visium_dir = sample.path[4], #ct04
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = instr)
s04 <- addCellMetadata(s04, new_metadata = rep("CTRL", ncol(s04)), vector_name = "type")

# Quality control 
## Filter over tissue genes
in.tissue <- pDataDT(s04)[in_tissue == 1]$cell_ID 
if (sum((in.tissue == pDataDT(s04)$cell_ID) == F) == 0) {
  cat("All genes are over tissue")
  rm(in.tissue)
} else {
  s04 <- subsetGiotto(s04, cell_ids = in.tissue)
  rm(in.tissue)
} # All genes are over tissue

## Mito percentage
s04.matrix <- getExpression(s04, values = "raw", output = "matrix") # Extract expression matrix
is_mito <- grep("(^MT-)|(^Mt-)|(^mt-)", rownames(s04.matrix), value = T) # Identify mito genes
mito.sum <- apply(s04.matrix[is_mito, , drop = F], 2, sum)
expr.sum <- apply(s04.matrix, 2, sum)
mito.perc <- mito.sum / expr.sum * 100

s04 <- addCellMetadata(s04, new_metadata = as.numeric(mito.perc), vector_name = "mito_perc")
s04.meta <- pDataDT(s04); rm(is_mito, mito.sum, expr.sum, mito.perc)


## Create cell metadata
reads.depth <- apply(s04.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s04.matrix != 0, 2, sum) # Detected genes per cell
s04.meta$n_reads <- as.numeric(reads.depth); s04.meta$n_genes <- as.numeric(genes.per.cell)
rm(reads.depth, genes.per.cell) 

### Libsize vs Detected genes correlation
s04.cor <- as.numeric(round(with(s04.meta, cor.test(n_reads, n_genes))$estimate, 4)) 


## Plots
### Detected genes
s04.detected <- filterDistributions(s04, detection = "cells", nr_bins = 50,
                                    method = "threshold", expression_threshold = 1, 
                                    default_save_name = "detected_genes")

### Genes detection along cells
s04.det.along <- filterDistributions(s04, detection = "feats", nr_bins = 50,
                                     method = "threshold", 
                                     default_save_name = "library_size")

### Libsize vs Detected genes
s04.libvdet <- ggplot(s04.meta, aes(n_reads, n_genes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = paste0("Corr = ", s04.cor), x = 25000, y = 0); rm(s04.cor)

### Threshold evaluation
s04.thresholds <- filterCombinations(s04, expression_thresholds = 1,
                                     feat_det_in_min_cells = rep(round(ncol(s04.matrix) / 20, 0), 4),
                                     min_det_feats_per_cell = c(125, 250, 500, 1000), 
                                     show_plot = F,
                                     default_save_name = "thresholds")

## Filter sample
s04.filtered <- filterGiotto(s04, expression_values = "raw", expression_threshold = 1,
                             feat_det_in_min_cells = round(ncol(s04.matrix) / 20, 0),
                             min_det_feats_per_cell = 1000)

### Visualization plot
s04.spots.plot <- spatPlot2D(s04, cell_color = ("lightgrey"), point_size = 2,
                             select_cells = s04.filtered@cell_ID$cell, # Kept spots
                             other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                             title = "")

s04.metric.plots <- ggarrange(s04.detected, s04.det.along, 
                              s04.libvdet, s04.spots.plot) # Combine all plots
rm(s04.detected, s04.det.along, s04.libvdet, s04.thresholds); s04.metric.plots

ggsave("./project/outcomes/qc/ct04_metrics.jpg", plot = s04.metric.plots, 
       scale = 2.5, width = 1920, height = 1080, units = "px")

####--SAMPLE 05--####
# Obtain data
s05 <- createGiottoVisiumObject(
  visium_dir = sample.path[5], #ct05
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = instr)
s05 <- addCellMetadata(s05, new_metadata = rep("CTRL", ncol(s05)), vector_name = "type")

# Quality control 
## Filter over tissue genes
in.tissue <- pDataDT(s05)[in_tissue == 1]$cell_ID 
if (sum((in.tissue == pDataDT(s05)$cell_ID) == F) == 0) {
  cat("All genes are over tissue")
  rm(in.tissue)
} else {
  s05 <- subsetGiotto(s05, cell_ids = in.tissue)
  rm(in.tissue)
} # All genes are over tissue

## Mito percentage
s05.matrix <- getExpression(s05, values = "raw", output = "matrix") # Extract expression matrix
is_mito <- grep("(^MT-)|(^Mt-)|(^mt-)", rownames(s05.matrix), value = T) # Identify mito genes
mito.sum <- apply(s05.matrix[is_mito, , drop = F], 2, sum)
expr.sum <- apply(s05.matrix, 2, sum)
mito.perc <- mito.sum / expr.sum * 100

s05 <- addCellMetadata(s05, new_metadata = as.numeric(mito.perc), vector_name = "mito_perc")
s05.meta <- pDataDT(s05); rm(is_mito, mito.sum, expr.sum, mito.perc)


## Create cell metadata
reads.depth <- apply(s05.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s05.matrix != 0, 2, sum) # Detected genes per cell
s05.meta$n_reads <- as.numeric(reads.depth); s05.meta$n_genes <- as.numeric(genes.per.cell)
rm(reads.depth, genes.per.cell) 

### Libsize vs Detected genes correlation
s05.cor <- as.numeric(round(with(s05.meta, cor.test(n_reads, n_genes))$estimate, 4)) 


## Plots
### Detected genes
s05.detected <- filterDistributions(s05, detection = "cells", nr_bins = 50,
                                    method = "threshold", expression_threshold = 1, 
                                    default_save_name = "detected_genes")

### Genes detection along cells
s05.det.along <- filterDistributions(s05, detection = "feats", nr_bins = 50,
                                     method = "threshold", 
                                     default_save_name = "library_size")

### Libsize vs Detected genes
s05.libvdet <- ggplot(s05.meta, aes(n_reads, n_genes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = paste0("Corr = ", s05.cor), x = 13000, y = 0); rm(s05.cor)

### Threshold evaluation
s05.thresholds <- filterCombinations(s05, expression_thresholds = 1,
                                     feat_det_in_min_cells = rep(round(ncol(s05.matrix) / 20, 0), 4),
                                     min_det_feats_per_cell = c(125, 250, 500, 1000), 
                                     show_plot = F,
                                     default_save_name = "thresholds")

## Filter sample
s05.filtered <- filterGiotto(s05, expression_values = "raw", expression_threshold = 1,
                             feat_det_in_min_cells = round(ncol(s05.matrix) / 20, 0),
                             min_det_feats_per_cell = 500)

### Visualization plot
s05.spots.plot <- spatPlot2D(s05, cell_color = ("lightgrey"), point_size = 2,
                             select_cells = s05.filtered@cell_ID$cell, # Kept spots
                             other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                             title = "")

s05.metric.plots <- ggarrange(s05.detected, s05.det.along, 
                              s05.libvdet, s05.spots.plot) # Combine all plots
rm(s05.detected, s05.det.along, s05.libvdet, s05.thresholds); s05.metric.plots

ggsave("./project/outcomes/qc/ct05_metrics.jpg", plot = s05.metric.plots, 
       scale = 2.5, width = 1920, height = 1080, units = "px")

####--SAMPLE 06--####
# Obtain data
s06 <- createGiottoVisiumObject(
  visium_dir = sample.path[6], #ct06
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = instr)
s06 <- addCellMetadata(s06, new_metadata = rep("CTRL", ncol(s06)), vector_name = "type")

# Quality control 
## Filter over tissue genes
in.tissue <- pDataDT(s06)[in_tissue == 1]$cell_ID 
if (sum((in.tissue == pDataDT(s06)$cell_ID) == F) == 0) {
  cat("All genes are over tissue")
  rm(in.tissue)
} else {
  s06 <- subsetGiotto(s06, cell_ids = in.tissue)
  rm(in.tissue)
} # All genes are over tissue

## Mito percentage
s06.matrix <- getExpression(s06, values = "raw", output = "matrix") # Extract expression matrix
is_mito <- grep("(^MT-)|(^Mt-)|(^mt-)", rownames(s06.matrix), value = T) # Identify mito genes
mito.sum <- apply(s06.matrix[is_mito, , drop = F], 2, sum)
expr.sum <- apply(s06.matrix, 2, sum)
mito.perc <- mito.sum / expr.sum * 100

s06 <- addCellMetadata(s06, new_metadata = as.numeric(mito.perc), vector_name = "mito_perc")
s06.meta <- pDataDT(s06); rm(is_mito, mito.sum, expr.sum, mito.perc)


## Create cell metadata
reads.depth <- apply(s06.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s06.matrix != 0, 2, sum) # Detected genes per cell
s06.meta$n_reads <- as.numeric(reads.depth); s06.meta$n_genes <- as.numeric(genes.per.cell)
rm(reads.depth, genes.per.cell) 

### Libsize vs Detected genes correlation
s06.cor <- as.numeric(round(with(s06.meta, cor.test(n_reads, n_genes))$estimate, 4)) 


## Plots
### Detected genes
s06.detected <- filterDistributions(s06, detection = "cells", nr_bins = 50,
                                    method = "threshold", expression_threshold = 1, 
                                    default_save_name = "detected_genes")

### Genes detection along cells
s06.det.along <- filterDistributions(s06, detection = "feats", nr_bins = 50,
                                     method = "threshold", 
                                     default_save_name = "library_size")

### Libsize vs Detected genes
s06.libvdet <- ggplot(s06.meta, aes(n_reads, n_genes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = paste0("Corr = ", s06.cor), x = 7500, y = 0); rm(s06.cor)

### Threshold evaluation
s06.thresholds <- filterCombinations(s06, expression_thresholds = 1,
                                     feat_det_in_min_cells = rep(round(ncol(s06.matrix) / 20, 0), 4),
                                     min_det_feats_per_cell = c(125, 250, 500, 1000), 
                                     show_plot = F,
                                     default_save_name = "thresholds")

## Filter sample
s06.filtered <- filterGiotto(s06, expression_values = "raw", expression_threshold = 1,
                             feat_det_in_min_cells = round(ncol(s06.matrix) / 20, 0),
                             min_det_feats_per_cell = 250)

### Visualization plot
s06.spots.plot <- spatPlot2D(s06, cell_color = ("lightgrey"), point_size = 2,
                             select_cells = s06.filtered@cell_ID$cell, # Kept spots
                             other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                             title = "")

s06.metric.plots <- ggarrange(s06.detected, s06.det.along, 
                              s06.libvdet, s06.spots.plot) # Combine all plots
rm(s06.detected, s06.det.along, s06.libvdet, s06.thresholds); s06.metric.plots

ggsave("./project/outcomes/qc/ct06_metrics.jpg", plot = s06.metric.plots, 
       scale = 2.5, width = 1920, height = 1080, units = "px")

####--SAMPLE 07--####
# Obtain data
s07 <- createGiottoVisiumObject(
  visium_dir = sample.path[7], #ms01
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = instr)
s07 <- addCellMetadata(s07, new_metadata = rep("MSCA", ncol(s07)), vector_name = "type")

# Quality control
## Filter over tissue genes
in.tissue <- pDataDT(s07)[in_tissue == 1]$cell_ID
if (sum((in.tissue == pDataDT(s07)$cell_ID) == F) == 0) {
  cat("All genes are over tissue")
  rm(in.tissue)
} else {
  s07 <- subsetGiotto(s07, cell_ids = in.tissue)
  rm(in.tissue)
} # All genes are over tissue

## Mito percentage
s07.matrix <- getExpression(s07, values = "raw", output = "matrix") # Extract expression matrix
is_mito <- grep("(^MT-)|(^Mt-)|(^mt-)", rownames(s07.matrix), value = T) # Identify mito genes
mito.sum <- apply(s07.matrix[is_mito, , drop = F], 2, sum)
expr.sum <- apply(s07.matrix, 2, sum)
mito.perc <- mito.sum / expr.sum * 100

s07 <- addCellMetadata(s07, new_metadata = as.numeric(mito.perc), vector_name = "mito_perc")
s07.meta <- pDataDT(s07); rm(is_mito, mito.sum, expr.sum, mito.perc)


## Create cell metadata
reads.depth <- apply(s07.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s07.matrix != 0, 2, sum) # Detected genes per cell
s07.meta$n_reads <- as.numeric(reads.depth); s07.meta$n_genes <- as.numeric(genes.per.cell)
rm(reads.depth, genes.per.cell)

### Libsize vs Detected genes correlation
s07.cor <- as.numeric(round(with(s07.meta, cor.test(n_reads, n_genes))$estimate, 4))


## Plots
### Detected genes
s07.detected <- filterDistributions(s07, detection = "cells", nr_bins = 50,
                                    method = "threshold", expression_threshold = 1,
                                    default_save_name = "detected_genes")

### Genes detection along cells
s07.det.along <- filterDistributions(s07, detection = "feats", nr_bins = 50,
                                     method = "threshold",
                                     default_save_name = "library_size")

### Libsize vs Detected genes
s07.libvdet <- ggplot(s07.meta, aes(n_reads, n_genes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") +
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = paste0("Corr = ", s07.cor), x = 9000, y = 0); rm(s07.cor)

### Threshold evaluation
s07.thresholds <- filterCombinations(s07, expression_thresholds = 1,
                                     feat_det_in_min_cells = rep(round(ncol(s07.matrix) / 20, 0), 4),
                                     min_det_feats_per_cell = c(125, 250, 500, 1000),
                                     show_plot = F,
                                     default_save_name = "thresholds")

## Filter sample
s07.filtered <- filterGiotto(s07, expression_values = "raw", expression_threshold = 1,
                             feat_det_in_min_cells = round(ncol(s07.matrix) / 20, 0),
                             min_det_feats_per_cell = 1000)

### Visualization plot
s07.spots.plot <- spatPlot2D(s07, cell_color = ("lightgrey"), point_size = 2,
                             select_cells = s07.filtered@cell_ID$cell, # Kept spots
                             other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                             title = "")

s07.metric.plots <- ggarrange(s07.detected, s07.det.along, 
                              s07.libvdet, s07.spots.plot) # Combine all plots
rm(s07.detected, s07.det.along, s07.libvdet, s07.thresholds); s07.metric.plots

ggsave("./project/outcomes/qc/ms01_metrics.jpg", plot = s07.metric.plots, 
       scale = 2.5, width = 1920, height = 1080, units = "px")

####--SAMPLE 08--####
# Obtain data
s08 <- createGiottoVisiumObject(
  visium_dir = sample.path[8], #ms02
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = instr)
s08 <- addCellMetadata(s08, new_metadata = rep("MSCA", ncol(s08)), vector_name = "type")

# Quality control 
## Filter over tissue genes
in.tissue <- pDataDT(s08)[in_tissue == 1]$cell_ID 
if (sum((in.tissue == pDataDT(s08)$cell_ID) == F) == 0) {
  cat("All genes are over tissue")
  rm(in.tissue)
} else {
  s08 <- subsetGiotto(s08, cell_ids = in.tissue)
  rm(in.tissue)
} # All genes are over tissue

## Mito percentage
s08.matrix <- getExpression(s08, values = "raw", output = "matrix") # Extract expression matrix
is_mito <- grep("(^MT-)|(^Mt-)|(^mt-)", rownames(s08.matrix), value = T) # Identify mito genes
mito.sum <- apply(s08.matrix[is_mito, , drop = F], 2, sum)
expr.sum <- apply(s08.matrix, 2, sum)
mito.perc <- mito.sum / expr.sum * 100

s08 <- addCellMetadata(s08, new_metadata = as.numeric(mito.perc), vector_name = "mito_perc")
s08.meta <- pDataDT(s08); rm(is_mito, mito.sum, expr.sum, mito.perc)


## Create cell metadata
reads.depth <- apply(s08.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s08.matrix != 0, 2, sum) # Detected genes per cell
s08.meta$n_reads <- as.numeric(reads.depth); s08.meta$n_genes <- as.numeric(genes.per.cell)
rm(reads.depth, genes.per.cell) 

### Libsize vs Detected genes correlation
s08.cor <- as.numeric(round(with(s08.meta, cor.test(n_reads, n_genes))$estimate, 4)) 


## Plots
### Detected genes
s08.detected <- filterDistributions(s08, detection = "cells", nr_bins = 50,
                                    method = "threshold", expression_threshold = 1, 
                                    default_save_name = "detected_genes")

### Genes detection along cells
s08.det.along <- filterDistributions(s08, detection = "feats", nr_bins = 50,
                                     method = "threshold", 
                                     default_save_name = "library_size")

### Libsize vs Detected genes
s08.libvdet <- ggplot(s08.meta, aes(n_reads, n_genes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = paste0("Corr = ", s08.cor), x = 22000, y = 0); rm(s08.cor)

### Threshold evaluation
s08.thresholds <- filterCombinations(s08, expression_thresholds = 1,
                                     feat_det_in_min_cells = rep(round(ncol(s08.matrix) / 20, 0),
                                                                 4),
                                     min_det_feats_per_cell = c(125, 250, 500, 1000), 
                                     show_plot = F,
                                     default_save_name = "thresholds")

## Filter sample
s08.filtered <- filterGiotto(s08, expression_values = "raw", expression_threshold = 1,
                             feat_det_in_min_cells = round(ncol(s08.matrix) / 20, 0),
                             min_det_feats_per_cell = 125)

### Visualization plot
s08.spots.plot <- spatPlot2D(s08, cell_color = ("lightgrey"), point_size = 2,
                             select_cells = s08.filtered@cell_ID$cell, # Kept spots
                             other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                             title = "")

s08.metric.plots <- ggarrange(s08.detected, s08.det.along, 
                              s08.libvdet, s08.spots.plot) # Combine all plots
rm(s08.detected, s08.det.along, s08.libvdet, s08.thresholds); s08.metric.plots

ggsave("./project/outcomes/qc/ms02_metrics.jpg", plot = s08.metric.plots, 
       scale = 2.5, width = 1920, height = 1080, units = "px")

####--SAMPLE 09--####
# Obtain data
s09 <- createGiottoVisiumObject(
  visium_dir = sample.path[9], #ms03
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = instr)
s09 <- addCellMetadata(s09, new_metadata = rep("MSCA", ncol(s09)), vector_name = "type")

# Quality control 
## Filter over tissue genes
in.tissue <- pDataDT(s09)[in_tissue == 1]$cell_ID 
if (sum((in.tissue == pDataDT(s09)$cell_ID) == F) == 0) {
  cat("All genes are over tissue")
  rm(in.tissue)
} else {
  s09 <- subsetGiotto(s09, cell_ids = in.tissue)
  rm(in.tissue)
} # All genes are over tissue

## Mito percentage
s09.matrix <- getExpression(s09, values = "raw", output = "matrix") # Extract expression matrix
is_mito <- grep("(^MT-)|(^Mt-)|(^mt-)", rownames(s09.matrix), value = T) # Identify mito genes
mito.sum <- apply(s09.matrix[is_mito, , drop = F], 2, sum)
expr.sum <- apply(s09.matrix, 2, sum)
mito.perc <- mito.sum / expr.sum * 100

s09 <- addCellMetadata(s09, new_metadata = as.numeric(mito.perc), vector_name = "mito_perc")
s09.meta <- pDataDT(s09); rm(is_mito, mito.sum, expr.sum, mito.perc)


## Create cell metadata
reads.depth <- apply(s09.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s09.matrix != 0, 2, sum) # Detected genes per cell
s09.meta$n_reads <- as.numeric(reads.depth); s09.meta$n_genes <- as.numeric(genes.per.cell)
rm(reads.depth, genes.per.cell) 

### Libsize vs Detected genes correlation
s09.cor <- as.numeric(round(with(s09.meta, cor.test(n_reads, n_genes))$estimate, 4)) 


## Plots
### Detected genes
s09.detected <- filterDistributions(s09, detection = "cells", nr_bins = 50,
                                    method = "threshold", expression_threshold = 1, 
                                    default_save_name = "detected_genes")

### Genes detection along cells
s09.det.along <- filterDistributions(s09, detection = "feats", nr_bins = 50,
                                     method = "threshold", 
                                     default_save_name = "library_size")

### Libsize vs Detected genes
s09.libvdet <- ggplot(s09.meta, aes(n_reads, n_genes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = paste0("Corr = ", s09.cor), x = 17000, y = 0); rm(s09.cor)

### Threshold evaluation
s09.thresholds <- filterCombinations(s09, expression_thresholds = 1,
                                     feat_det_in_min_cells = rep(round(ncol(s09.matrix) / 20, 0), 4),
                                     min_det_feats_per_cell = c(125, 250, 500, 1000), 
                                     show_plot = F,
                                     default_save_name = "thresholds")

## Filter sample
s09.filtered <- filterGiotto(s09, expression_values = "raw", expression_threshold = 1,
                             feat_det_in_min_cells = round(ncol(s09.matrix) / 20, 0),
                             min_det_feats_per_cell = 250)

### Visualization plot
s09.spots.plot <- spatPlot2D(s09, cell_color = ("lightgrey"), point_size = 2,
                             select_cells = s09.filtered@cell_ID$cell, # Kept spots
                             other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                             title = "")

s09.metric.plots <- ggarrange(s09.detected, s09.det.along, 
                              s09.libvdet, s09.spots.plot) # Combine all plots
rm(s09.detected, s09.det.along, s09.libvdet, s09.thresholds); s09.metric.plots

ggsave("./project/outcomes/qc/ms03_metrics.jpg", plot = s09.metric.plots, 
       scale = 2.5, width = 1920, height = 1080, units = "px")

####--SAMPLE 10--####
# Obtain data
s10 <- createGiottoVisiumObject(
  visium_dir = sample.path[10], #ms04
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = instr)
s10 <- addCellMetadata(s10, new_metadata = rep("MSCA", ncol(s10)), vector_name = "type")

# Quality control 
## Filter over tissue genes
in.tissue <- pDataDT(s10)[in_tissue == 1]$cell_ID 
if (sum((in.tissue == pDataDT(s10)$cell_ID) == F) == 0) {
  cat("All genes are over tissue")
  rm(in.tissue)
} else {
  s10 <- subsetGiotto(s10, cell_ids = in.tissue)
  rm(in.tissue)
} # All genes are over tissue

## Mito percentage
s10.matrix <- getExpression(s10, values = "raw", output = "matrix") # Extract expression matrix
is_mito <- grep("(^MT-)|(^Mt-)|(^mt-)", rownames(s10.matrix), value = T) # Identify mito genes
mito.sum <- apply(s10.matrix[is_mito, , drop = F], 2, sum)
expr.sum <- apply(s10.matrix, 2, sum)
mito.perc <- mito.sum / expr.sum * 100

s10 <- addCellMetadata(s10, new_metadata = as.numeric(mito.perc), vector_name = "mito_perc")
s10.meta <- pDataDT(s10); rm(is_mito, mito.sum, expr.sum, mito.perc)


## Create cell metadata
reads.depth <- apply(s10.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s10.matrix != 0, 2, sum) # Detected genes per cell
s10.meta$n_reads <- as.numeric(reads.depth); s10.meta$n_genes <- as.numeric(genes.per.cell)
rm(reads.depth, genes.per.cell) 

### Libsize vs Detected genes correlation
s10.cor <- as.numeric(round(with(s10.meta, cor.test(n_reads, n_genes))$estimate, 4)) 


## Plots
### Detected genes
s10.detected <- filterDistributions(s10, detection = "cells", nr_bins = 50,
                                    method = "threshold", expression_threshold = 1, 
                                    default_save_name = "detected_genes")

### Genes detection along cells
s10.det.along <- filterDistributions(s10, detection = "feats", nr_bins = 50,
                                     method = "threshold", 
                                     default_save_name = "library_size")

### Libsize vs Detected genes
s10.libvdet <- ggplot(s10.meta, aes(n_reads, n_genes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = paste0("Corr = ", s10.cor), x = 17000, y = 0); rm(s10.cor)

### Threshold evaluation
s10.thresholds <- filterCombinations(s10, expression_thresholds = 1,
                                     feat_det_in_min_cells = rep(round(ncol(s10.matrix) / 20, 0), 4),
                                     min_det_feats_per_cell = c(125, 250, 500, 1000), 
                                     show_plot = F,
                                     default_save_name = "thresholds")

## Filter sample
s10.filtered <- filterGiotto(s10, expression_values = "raw", expression_threshold = 1,
                             feat_det_in_min_cells = round(ncol(s10.matrix) / 20, 0),
                             min_det_feats_per_cell = 250)

### Visualization plot
s10.spots.plot <- spatPlot2D(s10, cell_color = ("lightgrey"), point_size = 2,
                             select_cells = s10.filtered@cell_ID$cell, # Kept spots
                             other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                             title = "Deleted spots (sample 10)"); s10.spots.plot

s10.metric.plots <- ggarrange(s10.detected, s10.det.along, 
                              s10.libvdet, s10.spots.plot) # Combine all plots
rm(s10.detected, s10.det.along, s10.libvdet, s10.thresholds); s10.metric.plots

ggsave("./project/outcomes/qc/ms04_metrics.jpg", plot = s10.metric.plots, 
       scale = 2.5, width = 1920, height = 1080, units = "px")

####--SAMPLE 11--####
# Obtain data
s11 <- createGiottoVisiumObject(
  visium_dir = sample.path[11], #ms05
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = instr)
s11 <- addCellMetadata(s11, new_metadata = rep("MSCA", ncol(s11)), vector_name = "type")

# Quality control 
## Filter over tissue genes
in.tissue <- pDataDT(s11)[in_tissue == 1]$cell_ID 
if (sum((in.tissue == pDataDT(s11)$cell_ID) == F) == 0) {
  cat("All genes are over tissue")
  rm(in.tissue)
} else {
  s11 <- subsetGiotto(s11, cell_ids = in.tissue)
  rm(in.tissue)
} # All genes are over tissue

## Mito percentage
s11.matrix <- getExpression(s11, values = "raw", output = "matrix") # Extract expression matrix
is_mito <- grep("(^MT-)|(^Mt-)|(^mt-)", rownames(s11.matrix), value = T) # Identify mito genes
mito.sum <- apply(s11.matrix[is_mito, , drop = F], 2, sum)
expr.sum <- apply(s11.matrix, 2, sum)
mito.perc <- mito.sum / expr.sum * 100

s11 <- addCellMetadata(s11, new_metadata = as.numeric(mito.perc), vector_name = "mito_perc")
s11.meta <- pDataDT(s11); rm(is_mito, mito.sum, expr.sum, mito.perc)


## Create cell metadata
reads.depth <- apply(s11.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s11.matrix != 0, 2, sum) # Detected genes per cell
s11.meta$n_reads <- as.numeric(reads.depth); s11.meta$n_genes <- as.numeric(genes.per.cell)
rm(reads.depth, genes.per.cell) 

### Libsize vs Detected genes correlation
s11.cor <- as.numeric(round(with(s11.meta, cor.test(n_reads, n_genes))$estimate, 4)) 


## Plots
### Detected genes
s11.detected <- filterDistributions(s11, detection = "cells", nr_bins = 50,
                                    method = "threshold", expression_threshold = 1, 
                                    default_save_name = "detected_genes")

### Genes detection along cells
s11.det.along <- filterDistributions(s11, detection = "feats", nr_bins = 50,
                                     method = "threshold", 
                                     default_save_name = "library_size")

### Libsize vs Detected genes
s11.libvdet <- ggplot(s11.meta, aes(n_reads, n_genes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = paste0("Corr = ", s11.cor), x = 20000, y = 0); rm(s11.cor)

### Threshold evaluation
s11.thresholds <- filterCombinations(s11, expression_thresholds = 1,
                                     feat_det_in_min_cells = rep(round(ncol(s11.matrix) / 20, 0), 4),
                                     min_det_feats_per_cell = c(125, 250, 500, 1000), 
                                     show_plot = F,
                                     default_save_name = "thresholds")

## Filter sample
s11.filtered <- filterGiotto(s11, expression_values = "raw", expression_threshold = 1,
                             feat_det_in_min_cells = round(ncol(s10.matrix) / 20, 0),
                             min_det_feats_per_cell = 125)

### Visualization plot
s11.spots.plot <- spatPlot2D(s11, cell_color = ("lightgrey"), point_size = 2,
                             select_cells = s11.filtered@cell_ID$cell, # Kept spots
                             other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                             title = "")

s11.metric.plots <- ggarrange(s11.detected, s11.det.along, 
                              s11.libvdet, s11.spots.plot) # Combine all plots
rm(s11.detected, s11.det.along, s11.libvdet, s11.thresholds); s11.metric.plots

ggsave("./project/outcomes/qc/ms05_metrics.jpg", plot = s11.metric.plots, 
       scale = 2.5, width = 1920, height = 1080, units = "px")

####--SAMPLE 12--####
# Obtain data
s12 <- createGiottoVisiumObject(
  visium_dir = sample.path[12], #ms06
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = instr)
s12 <- addCellMetadata(s12, new_metadata = rep("MSCA", ncol(s12)), vector_name = "type")

# Quality control 
## Filter over tissue genes
in.tissue <- pDataDT(s12)[in_tissue == 1]$cell_ID 
if (sum((in.tissue == pDataDT(s12)$cell_ID) == F) == 0) {
  cat("All genes are over tissue")
  rm(in.tissue)
} else {
  s12 <- subsetGiotto(s12, cell_ids = in.tissue)
  rm(in.tissue)
} # All genes are over tissue

## Mito percentage
s12.matrix <- getExpression(s12, values = "raw", output = "matrix") # Extract expression matrix
is_mito <- grep("(^MT-)|(^Mt-)|(^mt-)", rownames(s12.matrix), value = T) # Identify mito genes
mito.sum <- apply(s12.matrix[is_mito, , drop = F], 2, sum)
expr.sum <- apply(s12.matrix, 2, sum)
mito.perc <- mito.sum / expr.sum * 100

s12 <- addCellMetadata(s12, new_metadata = as.numeric(mito.perc), vector_name = "mito_perc")
s12.meta <- pDataDT(s12); rm(is_mito, mito.sum, expr.sum, mito.perc)


## Create cell metadata
reads.depth <- apply(s12.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s12.matrix != 0, 2, sum) # Detected genes per cell
s12.meta$n_reads <- as.numeric(reads.depth); s12.meta$n_genes <- as.numeric(genes.per.cell)
rm(reads.depth, genes.per.cell) 

### Libsize vs Detected genes correlation
s12.cor <- as.numeric(round(with(s12.meta, cor.test(n_reads, n_genes))$estimate, 4)) 


## Plots
### Detected genes
s12.detected <- filterDistributions(s12, detection = "cells", nr_bins = 50,
                                    method = "threshold", expression_threshold = 1, 
                                    default_save_name = "detected_genes")

### Genes detection along cells
s12.det.along <- filterDistributions(s12, detection = "feats", nr_bins = 50,
                                     method = "threshold", 
                                     default_save_name = "library_size")

### Libsize vs Detected genes
s12.libvdet <- ggplot(s12.meta, aes(n_reads, n_genes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = paste0("Corr = ", s12.cor), x = 17000, y = 0); rm(s12.cor)

### Threshold evaluation
s12.thresholds <- filterCombinations(s12, expression_thresholds = 1,
                                     feat_det_in_min_cells = rep(round(ncol(s12.matrix) / 20, 0), 4),
                                     min_det_feats_per_cell = c(125, 250, 500, 1000), 
                                     show_plot = F,
                                     default_save_name = "thresholds")

## Filter sample
s12.filtered <- filterGiotto(s12, expression_values = "raw", expression_threshold = 1,
                             feat_det_in_min_cells = round(ncol(s12.matrix) / 20, 0),
                             min_det_feats_per_cell = 250)

### Visualization plot
s12.spots.plot <- spatPlot2D(s12, cell_color = ("lightgrey"), point_size = 2,
                             select_cells = s12.filtered@cell_ID$cell, # Kept spots
                             other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                             title = "")

s12.metric.plots <- ggarrange(s12.detected, s12.det.along, 
                              s12.libvdet, s12.spots.plot) # Combine all plots
rm(s12.detected, s12.det.along, s12.libvdet, s12.thresholds); s12.metric.plots

ggsave("./project/outcomes/qc/ms06_metrics.jpg", plot = s12.metric.plots, 
       scale = 2.5, width = 1920, height = 1080, units = "px")

####--SAMPLE 13--####
# Obtain data
s13 <- createGiottoVisiumObject(
  visium_dir = sample.path[13], #ms07
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = instr)
s13 <- addCellMetadata(s13, new_metadata = rep("MSCA", ncol(s13)), vector_name = "type")

# Quality control 
## Filter over tissue genes
in.tissue <- pDataDT(s13)[in_tissue == 1]$cell_ID 
if (sum((in.tissue == pDataDT(s13)$cell_ID) == F) == 0) {
  cat("All genes are over tissue")
  rm(in.tissue)
} else {
  s13 <- subsetGiotto(s13, cell_ids = in.tissue)
  rm(in.tissue)
} # All genes are over tissue

## Mito percentage
s13.matrix <- getExpression(s13, values = "raw", output = "matrix") # Extract expression matrix
is_mito <- grep("(^MT-)|(^Mt-)|(^mt-)", rownames(s13.matrix), value = T) # Identify mito genes
mito.sum <- apply(s13.matrix[is_mito, , drop = F], 2, sum)
expr.sum <- apply(s13.matrix, 2, sum)
mito.perc <- mito.sum / expr.sum * 100

s13 <- addCellMetadata(s13, new_metadata = as.numeric(mito.perc), vector_name = "mito_perc")
s13.meta <- pDataDT(s13); rm(is_mito, mito.sum, expr.sum, mito.perc)


## Create cell metadata
reads.depth <- apply(s13.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s13.matrix != 0, 2, sum) # Detected genes per cell
s13.meta$n_reads <- as.numeric(reads.depth); s13.meta$n_genes <- as.numeric(genes.per.cell)
rm(reads.depth, genes.per.cell) 

### Libsize vs Detected genes correlation
s13.cor <- as.numeric(round(with(s13.meta, cor.test(n_reads, n_genes))$estimate, 4)) 


## Plots
### Detected genes
s13.detected <- filterDistributions(s13, detection = "cells", nr_bins = 50,
                                    method = "threshold", expression_threshold = 1, 
                                    default_save_name = "detected_genes")

### Genes detection along cells
s13.det.along <- filterDistributions(s13, detection = "feats", nr_bins = 50,
                                     method = "threshold", 
                                     default_save_name = "library_size")

### Libsize vs Detected genes
s13.libvdet <- ggplot(s13.meta, aes(n_reads, n_genes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = paste0("Corr = ", s13.cor), x = 17000, y = 0); rm(s13.cor)

### Threshold evaluation
s13.thresholds <- filterCombinations(s13, expression_thresholds = 1,
                                     feat_det_in_min_cells = rep(round(ncol(s13.matrix) / 20, 0), 4),
                                     min_det_feats_per_cell = c(125, 250, 500, 1000), 
                                     show_plot = F,
                                     default_save_name = "thresholds")

## Filter sample
s13.filtered <- filterGiotto(s13, expression_values = "raw", expression_threshold = 1,
                             feat_det_in_min_cells = round(ncol(s13.matrix) / 20, 0),
                             min_det_feats_per_cell = 250)

### Visualization plot
s13.spots.plot <- spatPlot2D(s13, cell_color = ("lightgrey"), point_size = 2,
                             select_cells = s13.filtered@cell_ID$cell, # Kept spots
                             other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                             title = "")

s13.metric.plots <- ggarrange(s13.detected, s13.det.along, 
                              s13.libvdet, s13.spots.plot) # Combine all plots
rm(s13.detected, s13.det.along, s13.libvdet, s13.thresholds); s13.metric.plots

ggsave("./project/outcomes/qc/ms07_metrics.jpg", plot = s13.metric.plots, 
       scale = 2.5, width = 1920, height = 1080, units = "px")

####--SAMPLE 14--####
# Obtain data
s14 <- createGiottoVisiumObject(
  visium_dir = sample.path[14], #ms08
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = instr)
s14 <- addCellMetadata(s14, new_metadata = rep("MSCA", ncol(s14)), vector_name = "type")

# Quality control 
## Filter over tissue genes
in.tissue <- pDataDT(s14)[in_tissue == 1]$cell_ID 
if (sum((in.tissue == pDataDT(s14)$cell_ID) == F) == 0) {
  cat("All genes are over tissue")
  rm(in.tissue)
} else {
  s14 <- subsetGiotto(s14, cell_ids = in.tissue)
  rm(in.tissue)
} # All genes are over tissue

## Mito percentage
s14.matrix <- getExpression(s14, values = "raw", output = "matrix") # Extract expression matrix
is_mito <- grep("(^MT-)|(^Mt-)|(^mt-)", rownames(s14.matrix), value = T) # Identify mito genes
mito.sum <- apply(s14.matrix[is_mito, , drop = F], 2, sum)
expr.sum <- apply(s14.matrix, 2, sum)
mito.perc <- mito.sum / expr.sum * 100

s14 <- addCellMetadata(s14, new_metadata = as.numeric(mito.perc), vector_name = "mito_perc")
s14.meta <- pDataDT(s14); rm(is_mito, mito.sum, expr.sum, mito.perc)


## Create cell metadata
reads.depth <- apply(s14.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s14.matrix != 0, 2, sum) # Detected genes per cell
s14.meta$n_reads <- as.numeric(reads.depth); s14.meta$n_genes <- as.numeric(genes.per.cell)
rm(reads.depth, genes.per.cell) 

### Libsize vs Detected genes correlation
s14.cor <- as.numeric(round(with(s14.meta, cor.test(n_reads, n_genes))$estimate, 4)) 


## Plots
### Detected genes
s14.detected <- filterDistributions(s14, detection = "cells", nr_bins = 50,
                                    method = "threshold", expression_threshold = 1, 
                                    default_save_name = "detected_genes")

### Genes detection along cells
s14.det.along <- filterDistributions(s14, detection = "feats", nr_bins = 50,
                                     method = "threshold", 
                                     default_save_name = "library_size")

### Libsize vs Detected genes
s14.libvdet <- ggplot(s14.meta, aes(n_reads, n_genes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = paste0("Corr = ", s14.cor), x = 14000, y = 0); rm(s14.cor)

### Threshold evaluation
s14.thresholds <- filterCombinations(s14, expression_thresholds = 1,
                                     feat_det_in_min_cells = rep(round(ncol(s14.matrix) / 20, 0), 4),
                                     min_det_feats_per_cell = c(125, 250, 500, 1000), 
                                     show_plot = F,
                                     default_save_name = "thresholds")

## Filter sample
s14.filtered <- filterGiotto(s14, expression_values = "raw", expression_threshold = 1,
                             feat_det_in_min_cells = round(ncol(s14.matrix) / 20, 0),
                             min_det_feats_per_cell = 500)

### Visualization plot
s14.spots.plot <- spatPlot2D(s14, cell_color = ("lightgrey"), point_size = 2,
                             select_cells = s14.filtered@cell_ID$cell, # Kept spots
                             other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                             title = "")

s14.metric.plots <- ggarrange(s14.detected, s14.det.along, 
                              s14.libvdet, s14.spots.plot) # Combine all plots
rm(s14.detected, s14.det.along, s14.libvdet, s14.thresholds); s14.metric.plots

ggsave("./project/outcomes/qc/ms08_metrics.jpg", plot = s14.metric.plots, 
       scale = 2.5, width = 1920, height = 1080, units = "px")

####--SAMPLE 15--####
# Obtain data
s15 <- createGiottoVisiumObject(
  visium_dir = sample.path[15], #ms09
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = instr)
s15 <- addCellMetadata(s15, new_metadata = rep("MSCI", ncol(s15)), vector_name = "type")

# Quality control 
## Filter over tissue genes
in.tissue <- pDataDT(s15)[in_tissue == 1]$cell_ID 
if (sum((in.tissue == pDataDT(s15)$cell_ID) == F) == 0) {
  cat("All genes are over tissue")
  rm(in.tissue)
} else {
  s15 <- subsetGiotto(s15, cell_ids = in.tissue)
  rm(in.tissue)
} # All genes are over tissue

## Mito percentage
s15.matrix <- getExpression(s15, values = "raw", output = "matrix") # Extract expression matrix
is_mito <- grep("(^MT-)|(^Mt-)|(^mt-)", rownames(s15.matrix), value = T) # Identify mito genes
mito.sum <- apply(s15.matrix[is_mito, , drop = F], 2, sum)
expr.sum <- apply(s15.matrix, 2, sum)
mito.perc <- mito.sum / expr.sum * 100

s15 <- addCellMetadata(s15, new_metadata = as.numeric(mito.perc), vector_name = "mito_perc")
s15.meta <- pDataDT(s15); rm(is_mito, mito.sum, expr.sum, mito.perc)


## Create cell metadata
reads.depth <- apply(s15.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s15.matrix != 0, 2, sum) # Detected genes per cell
s15.meta$n_reads <- as.numeric(reads.depth); s15.meta$n_genes <- as.numeric(genes.per.cell)
rm(reads.depth, genes.per.cell) 

### Libsize vs Detected genes correlation
s15.cor <- as.numeric(round(with(s15.meta, cor.test(n_reads, n_genes))$estimate, 4)) 


## Plots
### Detected genes
s15.detected <- filterDistributions(s15, detection = "cells", nr_bins = 50,
                                    method = "threshold", expression_threshold = 1, 
                                    default_save_name = "detected_genes")

### Genes detection along cells
s15.det.along <- filterDistributions(s15, detection = "feats", nr_bins = 50,
                                     method = "threshold", 
                                     default_save_name = "library_size")

### Libsize vs Detected genes
s15.libvdet <- ggplot(s15.meta, aes(n_reads, n_genes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = paste0("Corr = ", s15.cor), x = 9000, y = 0); rm(s15.cor)

### Threshold evaluation
s15.thresholds <- filterCombinations(s15, expression_thresholds = 1,
                                     feat_det_in_min_cells = rep(round(ncol(s15.matrix) / 20, 0), 4),
                                     min_det_feats_per_cell = c(125, 250, 500, 1000), 
                                     show_plot = F,
                                     default_save_name = "thresholds")

## Filter sample
s15.filtered <- filterGiotto(s15, expression_values = "raw", expression_threshold = 1,
                             feat_det_in_min_cells = round(ncol(s15.matrix) / 20, 0),
                             min_det_feats_per_cell = 500)

### Visualization plot
s15.spots.plot <- spatPlot2D(s15, cell_color = ("lightgrey"), point_size = 2,
                             select_cells = s15.filtered@cell_ID$cell, # Kept spots
                             other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                             title = "Deleted spots (sample 15)"); s15.spots.plot

s15.metric.plots <- ggarrange(s15.detected, s15.det.along, 
                              s15.libvdet, s15.spots.plot) # Combine all plots
rm(s15.detected, s15.det.along, s15.libvdet, s15.thresholds); s15.metric.plots

ggsave("./project/outcomes/qc/ms09_metrics.jpg", plot = s15.metric.plots, 
       scale = 2.5, width = 1920, height = 1080, units = "px")

####--SAMPLE 16--####
# Obtain data
s16 <- createGiottoVisiumObject(
  visium_dir = sample.path[16], #ms10
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = instr)
s16 <- addCellMetadata(s16, new_metadata = rep("MSCI", ncol(s16)), vector_name = "type")

# Quality control 
## Filter over tissue genes
in.tissue <- pDataDT(s16)[in_tissue == 1]$cell_ID 
if (sum((in.tissue == pDataDT(s16)$cell_ID) == F) == 0) {
  cat("All genes are over tissue")
  rm(in.tissue)
} else {
  s16 <- subsetGiotto(s16, cell_ids = in.tissue)
  rm(in.tissue)
} # All genes are over tissue

## Mito percentage
s16.matrix <- getExpression(s16, values = "raw", output = "matrix") # Extract expression matrix
is_mito <- grep("(^MT-)|(^Mt-)|(^mt-)", rownames(s16.matrix), value = T) # Identify mito genes
mito.sum <- apply(s16.matrix[is_mito, , drop = F], 2, sum)
expr.sum <- apply(s16.matrix, 2, sum)
mito.perc <- mito.sum / expr.sum * 100

s16 <- addCellMetadata(s16, new_metadata = as.numeric(mito.perc), vector_name = "mito_perc")
s16.meta <- pDataDT(s16); rm(is_mito, mito.sum, expr.sum, mito.perc)


## Create cell metadata
reads.depth <- apply(s16.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s16.matrix != 0, 2, sum) # Detected genes per cell
s16.meta$n_reads <- as.numeric(reads.depth); s16.meta$n_genes <- as.numeric(genes.per.cell)
rm(reads.depth, genes.per.cell) 

### Libsize vs Detected genes correlation
s16.cor <- as.numeric(round(with(s16.meta, cor.test(n_reads, n_genes))$estimate, 4)) 


## Plots
### Detected genes
s16.detected <- filterDistributions(s16, detection = "cells", nr_bins = 50,
                                    method = "threshold", expression_threshold = 1, 
                                    default_save_name = "detected_genes")

### Genes detection along cells
s16.det.along <- filterDistributions(s16, detection = "feats", nr_bins = 50,
                                     method = "threshold", 
                                     default_save_name = "library_size")

### Libsize vs Detected genes
s16.libvdet <- ggplot(s16.meta, aes(n_reads, n_genes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = paste0("Corr = ", s16.cor), x = 9000, y = 0); rm(s16.cor)

### Threshold evaluation
s16.thresholds <- filterCombinations(s16, expression_thresholds = 1,
                                     feat_det_in_min_cells = rep(round(ncol(s16.matrix) / 20, 0), 4),
                                     min_det_feats_per_cell = c(125, 250, 500, 1000), 
                                     show_plot = F,
                                     default_save_name = "thresholds")

## Filter sample
s16.filtered <- filterGiotto(s16, expression_values = "raw", expression_threshold = 1,
                             feat_det_in_min_cells = round(ncol(s16.matrix) / 20, 0),
                             min_det_feats_per_cell = 500)

### Visualization plot
s16.spots.plot <- spatPlot2D(s16, cell_color = ("lightgrey"), point_size = 2,
                             select_cells = s16.filtered@cell_ID$cell, # Kept spots
                             other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                             title = "")

s16.metric.plots <- ggarrange(s16.detected, s16.det.along, 
                              s16.libvdet, s16.spots.plot) # Combine all plots
rm(s16.detected, s16.det.along, s16.libvdet, s16.thresholds); s16.metric.plots

ggsave("./project/outcomes/qc/ms10_metrics.jpg", plot = s16.metric.plots, 
       scale = 2.5, width = 1920, height = 1080, units = "px")

####--SAMPLE 17--####
# Obtain data
s17 <- createGiottoVisiumObject(
  visium_dir = sample.path[17], #ms11
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = instr)
s17 <- addCellMetadata(s17, new_metadata = rep("MSCI", ncol(s17)), vector_name = "type")

# Quality control 
## Filter over tissue genes
in.tissue <- pDataDT(s17)[in_tissue == 1]$cell_ID 
if (sum((in.tissue == pDataDT(s17)$cell_ID) == F) == 0) {
  cat("All genes are over tissue")
  rm(in.tissue)
} else {
  s17 <- subsetGiotto(s17, cell_ids = in.tissue)
  rm(in.tissue)
} # All genes are over tissue

## Mito percentage
s17.matrix <- getExpression(s17, values = "raw", output = "matrix") # Extract expression matrix
is_mito <- grep("(^MT-)|(^Mt-)|(^mt-)", rownames(s17.matrix), value = T) # Identify mito genes
mito.sum <- apply(s17.matrix[is_mito, , drop = F], 2, sum)
expr.sum <- apply(s17.matrix, 2, sum)
mito.perc <- mito.sum / expr.sum * 100

s17 <- addCellMetadata(s17, new_metadata = as.numeric(mito.perc), vector_name = "mito_perc")
s17.meta <- pDataDT(s17); rm(is_mito, mito.sum, expr.sum, mito.perc)


## Create cell metadata
reads.depth <- apply(s17.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s17.matrix != 0, 2, sum) # Detected genes per cell
s17.meta$n_reads <- as.numeric(reads.depth); s17.meta$n_genes <- as.numeric(genes.per.cell)
rm(reads.depth, genes.per.cell) 

### Libsize vs Detected genes correlation
s17.cor <- as.numeric(round(with(s17.meta, cor.test(n_reads, n_genes))$estimate, 4)) 


## Plots
### Detected genes
s17.detected <- filterDistributions(s17, detection = "cells", nr_bins = 50,
                                    method = "threshold", expression_threshold = 1, 
                                    default_save_name = "detected_genes")

### Genes detection along cells
s17.det.along <- filterDistributions(s17, detection = "feats", nr_bins = 50,
                                     method = "threshold", 
                                     default_save_name = "library_size")

### Libsize vs Detected genes
s17.libvdet <- ggplot(s17.meta, aes(n_reads, n_genes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = paste0("Corr = ", s17.cor), x = 30000, y = 0); rm(s17.cor)

### Threshold evaluation
s17.thresholds <- filterCombinations(s17, expression_thresholds = 1,
                                     feat_det_in_min_cells = rep(round(ncol(s17.matrix) / 20, 0), 4),
                                     min_det_feats_per_cell = c(125, 250, 500, 1000), 
                                     show_plot = F,
                                     default_save_name = "thresholds")

## Filter sample
s17.filtered <- filterGiotto(s17, expression_values = "raw", expression_threshold = 1,
                             feat_det_in_min_cells = round(ncol(s17.matrix) / 20, 0),
                             min_det_feats_per_cell = 500)

### Visualization plot
s17.spots.plot <- spatPlot2D(s17, cell_color = ("lightgrey"), point_size = 2,
                             select_cells = s17.filtered@cell_ID$cell, # Kept spots
                             other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                             title = "")

s17.metric.plots <- ggarrange(s17.detected, s17.det.along, 
                              s17.libvdet, s17.spots.plot) # Combine all plots
rm(s17.detected, s17.det.along, s17.libvdet, s17.thresholds); s17.metric.plots

ggsave("./project/outcomes/qc/ms11_metrics.jpg", plot = s17.metric.plots, 
       scale = 2.5, width = 1920, height = 1080, units = "px")

####--SAMPLE 18--####
# Obtain data
s18 <- createGiottoVisiumObject(
  visium_dir = sample.path[18], #ms12
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = instr)
s18 <- addCellMetadata(s18, new_metadata = rep("MSCI", ncol(s18)), vector_name = "type")

# Quality control 
## Filter over tissue genes
in.tissue <- pDataDT(s18)[in_tissue == 1]$cell_ID 
if (sum((in.tissue == pDataDT(s18)$cell_ID) == F) == 0) {
  cat("All genes are over tissue")
  rm(in.tissue)
} else {
  s18 <- subsetGiotto(s18, cell_ids = in.tissue)
  rm(in.tissue)
} # All genes are over tissue

## Mito percentage
s18.matrix <- getExpression(s18, values = "raw", output = "matrix") # Extract expression matrix
is_mito <- grep("(^MT-)|(^Mt-)|(^mt-)", rownames(s18.matrix), value = T) # Identify mito genes
mito.sum <- apply(s18.matrix[is_mito, , drop = F], 2, sum)
expr.sum <- apply(s18.matrix, 2, sum)
mito.perc <- mito.sum / expr.sum * 100

s18 <- addCellMetadata(s18, new_metadata = as.numeric(mito.perc), vector_name = "mito_perc")
s18.meta <- pDataDT(s18); rm(is_mito, mito.sum, expr.sum, mito.perc)


## Create cell metadata
reads.depth <- apply(s18.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s18.matrix != 0, 2, sum) # Detected genes per cell
s18.meta$n_reads <- as.numeric(reads.depth); s18.meta$n_genes <- as.numeric(genes.per.cell)
rm(reads.depth, genes.per.cell) 

### Libsize vs Detected genes correlation
s18.cor <- as.numeric(round(with(s18.meta, cor.test(n_reads, n_genes))$estimate, 4)) 


## Plots
### Detected genes
s18.detected <- filterDistributions(s18, detection = "cells", nr_bins = 50,
                                    method = "threshold", expression_threshold = 1, 
                                    default_save_name = "detected_genes")

### Genes detection along cells
s18.det.along <- filterDistributions(s18, detection = "feats", nr_bins = 50,
                                     method = "threshold", 
                                     default_save_name = "library_size")

### Libsize vs Detected genes
s18.libvdet <- ggplot(s18.meta, aes(n_reads, n_genes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = paste0("Corr = ", s18.cor), x = 13000, y = 0); rm(s18.cor)

### Threshold evaluation
s18.thresholds <- filterCombinations(s18, expression_thresholds = 1,
                                     feat_det_in_min_cells = rep(round(ncol(s18.matrix) / 20, 0), 4),
                                     min_det_feats_per_cell = c(125, 250, 500, 1000), 
                                     show_plot = F,
                                     default_save_name = "thresholds")

## Filter sample
s18.filtered <- filterGiotto(s18, expression_values = "raw", expression_threshold = 1,
                             feat_det_in_min_cells = round(ncol(s18.matrix) / 20, 0),
                             min_det_feats_per_cell = 250)

### Visualization plot
s18.spots.plot <- spatPlot2D(s18, cell_color = ("lightgrey"), point_size = 2,
                             select_cells = s18.filtered@cell_ID$cell, # Kept spots
                             other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                             title = "")

s18.metric.plots <- ggarrange(s18.detected, s18.det.along, 
                              s18.libvdet, s18.spots.plot) # Combine all plots
rm(s18.detected, s18.det.along, s18.libvdet, s18.thresholds); s18.metric.plots

ggsave("./project/outcomes/qc/ms12_metrics.jpg", plot = s18.metric.plots, 
       scale = 2.5, width = 1920, height = 1080, units = "px")


####--SUMMARY--####
# Delete if more than 25% of cells are out
(length(s01@cell_ID$cell) - length(s01.filtered@cell_ID$cell)) / length(s01@cell_ID$cell)
(length(s02@cell_ID$cell) - length(s02.filtered@cell_ID$cell)) / length(s02@cell_ID$cell)
(length(s03@cell_ID$cell) - length(s03.filtered@cell_ID$cell)) / length(s03@cell_ID$cell)
(length(s04@cell_ID$cell) - length(s04.filtered@cell_ID$cell)) / length(s04@cell_ID$cell)
(length(s05@cell_ID$cell) - length(s05.filtered@cell_ID$cell)) / length(s05@cell_ID$cell)
(length(s06@cell_ID$cell) - length(s06.filtered@cell_ID$cell)) / length(s06@cell_ID$cell) # DELETE
(length(s07@cell_ID$cell) - length(s07.filtered@cell_ID$cell)) / length(s07@cell_ID$cell) # DELETE
(length(s08@cell_ID$cell) - length(s08.filtered@cell_ID$cell)) / length(s08@cell_ID$cell)
(length(s09@cell_ID$cell) - length(s09.filtered@cell_ID$cell)) / length(s09@cell_ID$cell)
(length(s10@cell_ID$cell) - length(s10.filtered@cell_ID$cell)) / length(s10@cell_ID$cell) # DELETE
(length(s11@cell_ID$cell) - length(s11.filtered@cell_ID$cell)) / length(s11@cell_ID$cell)
(length(s12@cell_ID$cell) - length(s12.filtered@cell_ID$cell)) / length(s12@cell_ID$cell)
(length(s13@cell_ID$cell) - length(s13.filtered@cell_ID$cell)) / length(s13@cell_ID$cell) # DELETE
(length(s14@cell_ID$cell) - length(s14.filtered@cell_ID$cell)) / length(s14@cell_ID$cell) # DELETE
(length(s15@cell_ID$cell) - length(s15.filtered@cell_ID$cell)) / length(s15@cell_ID$cell) # DELETE
(length(s16@cell_ID$cell) - length(s16.filtered@cell_ID$cell)) / length(s16@cell_ID$cell)
(length(s17@cell_ID$cell) - length(s17.filtered@cell_ID$cell)) / length(s17@cell_ID$cell)
(length(s18@cell_ID$cell) - length(s18.filtered@cell_ID$cell)) / length(s18@cell_ID$cell)


####--MERGE--####
s.list <- list(s01.filtered, s02.filtered, s03.filtered, s04.filtered, s05.filtered, 
               s08.filtered, s09.filtered, s11.filtered, s12.filtered,
               s16.filtered, s17.filtered, s18.filtered)
names <- c("CO37", "CO40", "CO41", "CO74", "CO85",
           "MS197D", "MS197U", "MS377N", "MS377T",
           "MS497T", "MS549H", "MS549T")

merged.samples <- joinGiottoObjects(gobject_list = s.list, gobject_names = names)
saveGiotto(merged.samples, foldername = "merged_sample", 
           dir = "./project/material/filtered_samples", overwrite = T)
