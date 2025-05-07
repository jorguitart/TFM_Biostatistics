####--LIBRARIES--####

library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(Giotto) #pak::pkg_install("drieslab/Giotto")
# library(spacexr) #devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
# library(nicheDE) #devtools::install_github("kaishumason/NicheDE")
# library(spatialGE) #devtools::install_github("fridleylab/spatialGE")


####--DIRECTORIES & INSTRUCTIONS--####

# Directories
dir <- "./project/material/GSE279181"
sam <- c("ct01", "ct02", "ct03", "ct04", "ct05", "ct06", 
         "ms01", "ms02", "ms03", "ms04", "ms05", "ms06", "ms07", "ms08", "ms09", "ms10", "ms11", "ms12")
sample.path <- file.path(dir, sam); rm(dir, sam)

# Giotto instructions
instr <- createGiottoInstructions(python_path = "C:/ProgramData/anaconda3/python.exe",
                                  show_plot = F, save_plot = F)


####--SAMPLE 01--####

# Obtain data
s01 <- createGiottoVisiumObject(
  visium_dir = sample.path[1], #ct01
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = instr)

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
s01.meta <- pDataDT(s01); rm(s01.matrix, is_mito, mito.sum, expr.sum, mito.perc)

### Keep low mito percentage cells
low.mito <- s01.meta$mito_perc <= 28
s01.nm <- s01[, low.mito]; rm(low.mito)

## Create cell metadata
s01.nm.matrix <- getExpression(s01.nm, values = "raw", output = "matrix")
reads.depth <- apply(s01.nm.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s01.nm.matrix != 0, 2, sum) # Detected genes per cell
s01.nm.meta <- pDataDT(s01.nm)
s01.nm.meta$n_reads <- as.numeric(reads.depth); s01.nm.meta$n_genes <- as.numeric(genes.per.cell)
rm(s01.nm.matrix, reads.depth, genes.per.cell)

### Libsize vs Detected genes correlation
s01.cor <- as.numeric(round(with(s01.nm.meta, cor.test(n_reads, n_genes))$estimate, 4)) 


## Plots
### Detected genes
s01.nm.detected <- filterDistributions(s01.nm, detection = "cells", nr_bins = 50,
                                       method = "threshold", expression_threshold = 1, 
                                       default_save_name = "detected_genes")

### Genes detection along cells
s01.nm.det.along <- filterDistributions(s01.nm, detection = "feats", nr_bins = 50,
                                        method = "threshold", 
                                        default_save_name = "library_size")

### Libsize vs Detected genes
s01.nm.libvdet <- ggplot(s01.nm.meta, aes(n_reads, n_genes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = paste0("Corr = ", s01.cor), x = 11000, y = 0); rm(s01.cor)

### Threshold evaluation
s01.nm.thresholds <- filterCombinations(s01.nm, expression_thresholds = 1,
                                    feat_det_in_min_cells = c(50, 100, 50, 100),
                                    min_det_feats_per_cell = c(500, 500, 1000, 1000), 
                                    show_plot = F,
                                    default_save_name = "thresholds")
s01.nm.thresholds[["ggplot"]][["theme"]][["legend.position"]] <- "none" # Remove legend

s01.nm.metric.plots <- ggarrange(s01.nm.detected, s01.nm.det.along, 
                                 s01.nm.libvdet, s01.nm.thresholds$ggplot) # Combine all plots
rm(s01.nm.detected, s01.nm.det.along, s01.nm.libvdet, s01.nm.thresholds); s01.nm.metric.plots
ggsave("./project/outcomes/qc/ct01_metrics.png", plot = s01.nm.metric.plots, 
       scale = 2, width = 1920, height = 1080, units = "px")

## Filter sample
s01.filtered <- filterGiotto(s01.nm, expression_values = "raw", expression_threshold = 1,
                            feat_det_in_min_cells = 50,
                            min_det_feats_per_cell = 1000)
saveRDS(s01.filtered, file = "./project/material/filtered_samples/s01_filtered.rds") # readRDS()

## Deleted spots
length(s01@cell_ID$cell) - length(s01.filtered@cell_ID$cell) 

### Visualization plot
s01.spots.plot <- spatPlot2D(s01, cell_color = ("lightgrey"), point_size = 2,
                            select_cells = s01.filtered@cell_ID$cell, # Kept spots
                            other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                            title = "Deleted spots (sample 1)"); s01.spots.plot
ggsave("./project/outcomes/qc/ct01_deleted.png", plot = s01.spots.plot, 
       scale = 2.5, width = 1920, height = 1080, units = "px")

####--SAMPLE 02--####

# Obtain data
s02 <- createGiottoVisiumObject(
  visium_dir = sample.path[2], #ct02
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = instr)

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
s02.meta <- pDataDT(s02); rm(s02.matrix, is_mito, mito.sum, expr.sum, mito.perc)

### Keep low mito percentage cells
low.mito <- s02.meta$mito_perc <= 28
s02.nm <- s02[, low.mito]; rm(low.mito)

## Create cell metadata
s02.nm.matrix <- getExpression(s02.nm, values = "raw", output = "matrix") 
reads.depth <- apply(s02.nm.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s02.nm.matrix != 0, 2, sum) # Detected genes per cell
s02.nm.meta <- pDataDT(s02.nm)
s02.nm.meta$n_reads <- as.numeric(reads.depth); s02.nm.meta$n_genes <- as.numeric(genes.per.cell)
rm(s02.nm.matrix, reads.depth, genes.per.cell) 

### Libsize vs Detected genes correlation
s02.cor <- as.numeric(round(with(s02.nm.meta, cor.test(n_reads, n_genes))$estimate, 4)) 


## Plots
### Detected genes
s02.nm.detected <- filterDistributions(s02.nm, detection = "cells", nr_bins = 50,
                                       method = "threshold", expression_threshold = 1, 
                                       default_save_name = "detected_genes")

### Genes detection along cells
s02.nm.det.along <- filterDistributions(s02.nm, detection = "feats", nr_bins = 50,
                                        method = "threshold", 
                                        default_save_name = "library_size")

### Libsize vs Detected genes
s02.nm.libvdet <- ggplot(s02.nm.meta, aes(n_reads, n_genes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = paste0("Corr = ", s02.cor), x = 11000, y = 0); rm(s02.cor)

### Threshold evaluation
s02.nm.thresholds <- filterCombinations(s02.nm, expression_thresholds = 1,
                                        feat_det_in_min_cells = c(50, 100, 50, 100),
                                        min_det_feats_per_cell = c(500, 500, 1000, 1000), 
                                        show_plot = F,
                                        default_save_name = "thresholds")
s02.nm.thresholds[["ggplot"]][["theme"]][["legend.position"]] <- "none" # Remove legend

s02.nm.metric.plots <- ggarrange(s02.nm.detected, s02.nm.det.along, 
                                 s02.nm.libvdet, s02.nm.thresholds$ggplot) # Combine all plots
rm(s02.nm.detected, s02.nm.det.along, s02.nm.libvdet, s02.nm.thresholds); s02.nm.metric.plots
ggsave("./project/outcomes/qc/ct02_metrics.png", plot = s02.nm.metric.plots, 
       scale = 2, width = 1920, height = 1080, units = "px")

## Filter sample
s02.filtered <- filterGiotto(s02.nm, expression_values = "raw", expression_threshold = 1,
                             feat_det_in_min_cells = 50,
                             min_det_feats_per_cell = 1000)
saveRDS(s02.filtered, file = "./project/material/filtered_samples/s02_filtered.rds") # readRDS()

## Deleted spots
length(s02@cell_ID$cell) - length(s02.filtered@cell_ID$cell) # 525 deleted spots

### Visualization plot
s02.spots.plot <- spatPlot2D(s02, cell_color = ("lightgrey"), point_size = 2,
                             select_cells = s02.filtered@cell_ID$cell, # Kept spots
                             other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                             title = "Deleted spots (sample 2)"); s02.spots.plot
ggsave("./project/outcomes/qc/ct02_deleted.png", plot = s02.spots.plot, 
       scale = 2.5, width = 1920, height = 1080, units = "px")

####--SAMPLE 03--####

# Obtain data
s03 <- createGiottoVisiumObject(
  visium_dir = sample.path[3], #ct03
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = instr)

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
s03.meta <- pDataDT(s03); rm(s03.matrix, is_mito, mito.sum, expr.sum, mito.perc)

### Keep low mito percentage cells
low.mito <- s03.meta$mito_perc <= 28
s03.nm <- s03[, low.mito]; rm(low.mito)

## Create cell metadata
s03.nm.matrix <- getExpression(s03.nm, values = "raw", output = "matrix") 
reads.depth <- apply(s03.nm.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s03.nm.matrix != 0, 2, sum) # Detected genes per cell
s03.nm.meta <- pDataDT(s03.nm)
s03.nm.meta$n_reads <- as.numeric(reads.depth); s03.nm.meta$n_genes <- as.numeric(genes.per.cell)
rm(s03.nm.matrix, reads.depth, genes.per.cell) 

### Libsize vs Detected genes correlation
s03.cor <- as.numeric(round(with(s03.nm.meta, cor.test(n_reads, n_genes))$estimate, 4)) 


## Plots
### Detected genes
s03.nm.detected <- filterDistributions(s03.nm, detection = "cells", nr_bins = 50,
                                       method = "threshold", expression_threshold = 1, 
                                       default_save_name = "detected_genes")

### Genes detection along cells
s03.nm.det.along <- filterDistributions(s03.nm, detection = "feats", nr_bins = 50,
                                        method = "threshold", 
                                        default_save_name = "library_size")

### Libsize vs Detected genes
s03.nm.libvdet <- ggplot(s03.nm.meta, aes(n_reads, n_genes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = paste0("Corr = ", s03.cor), x = 11000, y = 0); rm(s03.cor)

### Threshold evaluation
s03.nm.thresholds <- filterCombinations(s03.nm, expression_thresholds = 1,
                                        feat_det_in_min_cells = c(50, 100, 50, 100),
                                        min_det_feats_per_cell = c(500, 500, 1000, 1000), 
                                        show_plot = F,
                                        default_save_name = "thresholds")
s03.nm.thresholds[["ggplot"]][["theme"]][["legend.position"]] <- "none" # Remove legend

s03.nm.metric.plots <- ggarrange(s03.nm.detected, s03.nm.det.along, 
                                 s03.nm.libvdet, s03.nm.thresholds$ggplot) # Combine all plots
rm(s03.nm.detected, s03.nm.det.along, s03.nm.libvdet, s03.nm.thresholds); s03.nm.metric.plots
ggsave("./project/outcomes/qc/ct03_metrics.png", plot = s03.nm.metric.plots, 
       scale = 2, width = 1920, height = 1080, units = "px")

## Filter sample
s03.filtered <- filterGiotto(s03.nm, expression_values = "raw", expression_threshold = 1,
                             feat_det_in_min_cells = 50,
                             min_det_feats_per_cell = 1000)
saveRDS(s03.filtered, file = "./project/material/filtered_samples/s03_filtered.rds") # readRDS()

## Deleted spots
length(s03@cell_ID$cell) - length(s03.filtered@cell_ID$cell) # 525 deleted spots

### Visualization plot
s03.spots.plot <- spatPlot2D(s03, cell_color = ("lightgrey"), point_size = 2,
                             select_cells = s03.filtered@cell_ID$cell, # Kept spots
                             other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                             title = "Deleted spots (sample 3)"); s03.spots.plot
ggsave("./project/outcomes/qc/ct03_deleted.png", plot = s03.spots.plot, 
       scale = 2.5, width = 1920, height = 1080, units = "px")

####--SAMPLE 04--####

# Obtain data
s04 <- createGiottoVisiumObject(
  visium_dir = sample.path[4], #ct04
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = instr)

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
s04.meta <- pDataDT(s04); rm(s04.matrix, is_mito, mito.sum, expr.sum, mito.perc)

### Keep low mito percentage cells
low.mito <- s04.meta$mito_perc <= 28
s04.nm <- s04[, low.mito]; rm(low.mito)

## Create cell metadata
s04.nm.matrix <- getExpression(s04.nm, values = "raw", output = "matrix") 
reads.depth <- apply(s04.nm.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s04.nm.matrix != 0, 2, sum) # Detected genes per cell
s04.nm.meta <- pDataDT(s04.nm)
s04.nm.meta$n_reads <- as.numeric(reads.depth); s04.nm.meta$n_genes <- as.numeric(genes.per.cell)
rm(s04.nm.matrix, reads.depth, genes.per.cell) 

### Libsize vs Detected genes correlation
s04.cor <- as.numeric(round(with(s04.nm.meta, cor.test(n_reads, n_genes))$estimate, 4)) 


## Plots
### Detected genes
s04.nm.detected <- filterDistributions(s04.nm, detection = "cells", nr_bins = 50,
                                       method = "threshold", expression_threshold = 1, 
                                       default_save_name = "detected_genes")

### Genes detection along cells
s04.nm.det.along <- filterDistributions(s04.nm, detection = "feats", nr_bins = 50,
                                        method = "threshold", 
                                        default_save_name = "library_size")

### Libsize vs Detected genes
s04.nm.libvdet <- ggplot(s04.nm.meta, aes(n_reads, n_genes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = paste0("Corr = ", s04.cor), x = 11000, y = 0); rm(s04.cor)

### Threshold evaluation
s04.nm.thresholds <- filterCombinations(s04.nm, expression_thresholds = 1,
                                        feat_det_in_min_cells = c(50, 100, 50, 100),
                                        min_det_feats_per_cell = c(500, 500, 1000, 1000), 
                                        show_plot = F,
                                        default_save_name = "thresholds")
s04.nm.thresholds[["ggplot"]][["theme"]][["legend.position"]] <- "none" # Remove legend

s04.nm.metric.plots <- ggarrange(s04.nm.detected, s04.nm.det.along, 
                                 s04.nm.libvdet, s04.nm.thresholds$ggplot) # Combine all plots
rm(s04.nm.detected, s04.nm.det.along, s04.nm.libvdet, s04.nm.thresholds); s04.nm.metric.plots
ggsave("./project/outcomes/qc/ct04_metrics.png", plot = s04.nm.metric.plots, 
       scale = 2, width = 1920, height = 1080, units = "px")

## Filter sample
s04.filtered <- filterGiotto(s04.nm, expression_values = "raw", expression_threshold = 1,
                             feat_det_in_min_cells = 50,
                             min_det_feats_per_cell = 1000)
saveRDS(s04.filtered, file = "./project/material/filtered_samples/s04_filtered.rds") # readRDS()

## Deleted spots
length(s04@cell_ID$cell) - length(s04.filtered@cell_ID$cell) # 525 deleted spots

### Visualization plot
s04.spots.plot <- spatPlot2D(s04, cell_color = ("lightgrey"), point_size = 2,
                             select_cells = s04.filtered@cell_ID$cell, # Kept spots
                             other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                             title = "Deleted spots (sample 4)"); s04.spots.plot
ggsave("./project/outcomes/qc/ct04_deleted.png", plot = s04.spots.plot, 
       scale = 2.5, width = 1920, height = 1080, units = "px")

####--SAMPLE 05--####

# Obtain data
s05 <- createGiottoVisiumObject(
  visium_dir = sample.path[5], #ct05
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = instr)

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
s05.meta <- pDataDT(s05); rm(s05.matrix, is_mito, mito.sum, expr.sum, mito.perc)

### Keep low mito percentage cells
low.mito <- s05.meta$mito_perc <= 28
s05.nm <- s05[, low.mito]; rm(low.mito)

## Create cell metadata
s05.nm.matrix <- getExpression(s05.nm, values = "raw", output = "matrix") 
reads.depth <- apply(s05.nm.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s05.nm.matrix != 0, 2, sum) # Detected genes per cell
s05.nm.meta <- pDataDT(s05.nm)
s05.nm.meta$n_reads <- as.numeric(reads.depth); s05.nm.meta$n_genes <- as.numeric(genes.per.cell)
rm(s05.nm.matrix, reads.depth, genes.per.cell) 

### Libsize vs Detected genes correlation
s05.cor <- as.numeric(round(with(s05.nm.meta, cor.test(n_reads, n_genes))$estimate, 4)) 


## Plots
### Detected genes
s05.nm.detected <- filterDistributions(s05.nm, detection = "cells", nr_bins = 50,
                                       method = "threshold", expression_threshold = 1, 
                                       default_save_name = "detected_genes")

### Genes detection along cells
s05.nm.det.along <- filterDistributions(s05.nm, detection = "feats", nr_bins = 50,
                                        method = "threshold", 
                                        default_save_name = "library_size")

### Libsize vs Detected genes
s05.nm.libvdet <- ggplot(s05.nm.meta, aes(n_reads, n_genes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = paste0("Corr = ", s05.cor), x = 11000, y = 0); rm(s05.cor)

### Threshold evaluation
s05.nm.thresholds <- filterCombinations(s05.nm, expression_thresholds = 1,
                                        feat_det_in_min_cells = c(50, 100, 50, 100),
                                        min_det_feats_per_cell = c(500, 500, 1000, 1000), 
                                        show_plot = F,
                                        default_save_name = "thresholds")
s05.nm.thresholds[["ggplot"]][["theme"]][["legend.position"]] <- "none" # Remove legend

s05.nm.metric.plots <- ggarrange(s05.nm.detected, s05.nm.det.along, 
                                 s05.nm.libvdet, s05.nm.thresholds$ggplot) # Combine all plots
rm(s05.nm.detected, s05.nm.det.along, s05.nm.libvdet, s05.nm.thresholds); s05.nm.metric.plots
ggsave("./project/outcomes/qc/ct05_metrics.png", plot = s05.nm.metric.plots, 
       scale = 2, width = 1920, height = 1080, units = "px")

## Filter sample
s05.filtered <- filterGiotto(s05.nm, expression_values = "raw", expression_threshold = 1,
                             feat_det_in_min_cells = 50,
                             min_det_feats_per_cell = 1000)
saveRDS(s05.filtered, file = "./project/material/filtered_samples/s05_filtered.rds") # readRDS()

## Deleted spots
length(s05@cell_ID$cell) - length(s05.filtered@cell_ID$cell) # 525 deleted spots

### Visualization plot
s05.spots.plot <- spatPlot2D(s05, cell_color = ("lightgrey"), point_size = 2,
                             select_cells = s05.filtered@cell_ID$cell, # Kept spots
                             other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                             title = "Deleted spots (sample 5)"); s05.spots.plot
ggsave("./project/outcomes/qc/ct05_deleted.png", plot = s05.spots.plot, 
       scale = 2.5, width = 1920, height = 1080, units = "px")

####--SAMPLE 06--####

# Obtain data
s06 <- createGiottoVisiumObject(
  visium_dir = sample.path[6], #ct06
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = instr)

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
s06.meta <- pDataDT(s06); rm(s06.matrix, is_mito, mito.sum, expr.sum, mito.perc)

### Keep low mito percentage cells
low.mito <- s06.meta$mito_perc <= 28
s06.nm <- s06[, low.mito]; rm(low.mito)

## Create cell metadata
s06.nm.matrix <- getExpression(s06.nm, values = "raw", output = "matrix") 
reads.depth <- apply(s06.nm.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s06.nm.matrix != 0, 2, sum) # Detected genes per cell
s06.nm.meta <- pDataDT(s06.nm)
s06.nm.meta$n_reads <- as.numeric(reads.depth); s06.nm.meta$n_genes <- as.numeric(genes.per.cell)
rm(s06.nm.matrix, reads.depth, genes.per.cell) 

### Libsize vs Detected genes correlation
s06.cor <- as.numeric(round(with(s06.nm.meta, cor.test(n_reads, n_genes))$estimate, 4)) 


## Plots
### Detected genes
s06.nm.detected <- filterDistributions(s06.nm, detection = "cells", nr_bins = 50,
                                       method = "threshold", expression_threshold = 1, 
                                       default_save_name = "detected_genes")

### Genes detection along cells
s06.nm.det.along <- filterDistributions(s06.nm, detection = "feats", nr_bins = 50,
                                        method = "threshold", 
                                        default_save_name = "library_size")

### Libsize vs Detected genes
s06.nm.libvdet <- ggplot(s06.nm.meta, aes(n_reads, n_genes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = paste0("Corr = ", s06.cor), x = 11000, y = 0); rm(s06.cor)

### Threshold evaluation
s06.nm.thresholds <- filterCombinations(s06.nm, expression_thresholds = 1,
                                        feat_det_in_min_cells = c(50, 100, 50, 100),
                                        min_det_feats_per_cell = c(500, 500, 1000, 1000), 
                                        show_plot = F,
                                        default_save_name = "thresholds")
s06.nm.thresholds[["ggplot"]][["theme"]][["legend.position"]] <- "none" # Remove legend

s06.nm.metric.plots <- ggarrange(s06.nm.detected, s06.nm.det.along, 
                                 s06.nm.libvdet, s06.nm.thresholds$ggplot) # Combine all plots
rm(s06.nm.detected, s06.nm.det.along, s06.nm.libvdet, s06.nm.thresholds); s06.nm.metric.plots
ggsave("./project/outcomes/qc/ct06_metrics.png", plot = s06.nm.metric.plots, 
       scale = 2, width = 1920, height = 1080, units = "px")

## Filter sample
s06.filtered <- filterGiotto(s06.nm, expression_values = "raw", expression_threshold = 1,
                             feat_det_in_min_cells = 50,
                             min_det_feats_per_cell = 1000)
saveRDS(s06.filtered, file = "./project/material/filtered_samples/s06_filtered.rds") # readRDS()

## Deleted spots
length(s06@cell_ID$cell) - length(s06.filtered@cell_ID$cell) # 525 deleted spots

### Visualization plot
s06.spots.plot <- spatPlot2D(s06, cell_color = ("lightgrey"), point_size = 2,
                             select_cells = s06.filtered@cell_ID$cell, # Kept spots
                             other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                             title = "Deleted spots (sample 6)"); s06.spots.plot
ggsave("./project/outcomes/qc/ct06_deleted.png", plot = s06.spots.plot, 
       scale = 2.5, width = 1920, height = 1080, units = "px")

####--SAMPLE 07--####

# Obtain data
s07 <- createGiottoVisiumObject(
  visium_dir = sample.path[7], #ms01
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = instr)

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
s07.meta <- pDataDT(s07); rm(s07.matrix, is_mito, mito.sum, expr.sum, mito.perc)

### Keep low mito percentage cells
low.mito <- s07.meta$mito_perc <= 28
s07.nm <- s07[, low.mito]; rm(low.mito)

## Create cell metadata
s07.nm.matrix <- getExpression(s07.nm, values = "raw", output = "matrix") 
reads.depth <- apply(s07.nm.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s07.nm.matrix != 0, 2, sum) # Detected genes per cell
s07.nm.meta <- pDataDT(s07.nm)
s07.nm.meta$n_reads <- as.numeric(reads.depth); s07.nm.meta$n_genes <- as.numeric(genes.per.cell)
rm(s07.nm.matrix, reads.depth, genes.per.cell) 

### Libsize vs Detected genes correlation
s07.cor <- as.numeric(round(with(s07.nm.meta, cor.test(n_reads, n_genes))$estimate, 4)) 


## Plots
### Detected genes
s07.nm.detected <- filterDistributions(s07.nm, detection = "cells", nr_bins = 50,
                                       method = "threshold", expression_threshold = 1, 
                                       default_save_name = "detected_genes")

### Genes detection along cells
s07.nm.det.along <- filterDistributions(s07.nm, detection = "feats", nr_bins = 50,
                                        method = "threshold", 
                                        default_save_name = "library_size")

### Libsize vs Detected genes
s07.nm.libvdet <- ggplot(s07.nm.meta, aes(n_reads, n_genes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = paste0("Corr = ", s07.cor), x = 11000, y = 0); rm(s07.cor)

### Threshold evaluation
s07.nm.thresholds <- filterCombinations(s07.nm, expression_thresholds = 1,
                                        feat_det_in_min_cells = c(50, 100, 50, 100),
                                        min_det_feats_per_cell = c(500, 500, 1000, 1000), 
                                        show_plot = F,
                                        default_save_name = "thresholds")
s07.nm.thresholds[["ggplot"]][["theme"]][["legend.position"]] <- "none" # Remove legend

s07.nm.metric.plots <- ggarrange(s07.nm.detected, s07.nm.det.along, 
                                 s07.nm.libvdet, s07.nm.thresholds$ggplot) # Combine all plots
rm(s07.nm.detected, s07.nm.det.along, s07.nm.libvdet, s07.nm.thresholds); s07.nm.metric.plots
ggsave("./project/outcomes/qc/ms01_metrics.png", plot = s07.nm.metric.plots, 
       scale = 2, width = 1920, height = 1080, units = "px")

## Filter sample
s07.filtered <- filterGiotto(s07.nm, expression_values = "raw", expression_threshold = 1,
                             feat_det_in_min_cells = 50,
                             min_det_feats_per_cell = 1000)
saveRDS(s07.filtered, file = "./project/material/filtered_samples/s07_filtered.rds") # readRDS()

## Deleted spots
length(s07@cell_ID$cell) - length(s07.filtered@cell_ID$cell) # 525 deleted spots

### Visualization plot
s07.spots.plot <- spatPlot2D(s07, cell_color = ("lightgrey"), point_size = 2,
                             select_cells = s07.filtered@cell_ID$cell, # Kept spots
                             other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                             title = "Deleted spots (sample 7)"); s07.spots.plot
ggsave("./project/outcomes/qc/ms01_deleted.png", plot = s07.spots.plot, 
       scale = 2.5, width = 1920, height = 1080, units = "px")

####--SAMPLE 08--####

# Obtain data
s08 <- createGiottoVisiumObject(
  visium_dir = sample.path[8], #ms02
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = instr)

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
s08.meta <- pDataDT(s08); rm(s08.matrix, is_mito, mito.sum, expr.sum, mito.perc)

### Keep low mito percentage cells
low.mito <- s08.meta$mito_perc <= 28
s08.nm <- s08[, low.mito]; rm(low.mito)

## Create cell metadata
s08.nm.matrix <- getExpression(s08.nm, values = "raw", output = "matrix") 
reads.depth <- apply(s08.nm.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s08.nm.matrix != 0, 2, sum) # Detected genes per cell
s08.nm.meta <- pDataDT(s08.nm)
s08.nm.meta$n_reads <- as.numeric(reads.depth); s08.nm.meta$n_genes <- as.numeric(genes.per.cell)
rm(s08.nm.matrix, reads.depth, genes.per.cell) 

### Libsize vs Detected genes correlation
s08.cor <- as.numeric(round(with(s08.nm.meta, cor.test(n_reads, n_genes))$estimate, 4)) 


## Plots
### Detected genes
s08.nm.detected <- filterDistributions(s08.nm, detection = "cells", nr_bins = 50,
                                       method = "threshold", expression_threshold = 1, 
                                       default_save_name = "detected_genes")

### Genes detection along cells
s08.nm.det.along <- filterDistributions(s08.nm, detection = "feats", nr_bins = 50,
                                        method = "threshold", 
                                        default_save_name = "library_size")

### Libsize vs Detected genes
s08.nm.libvdet <- ggplot(s08.nm.meta, aes(n_reads, n_genes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = paste0("Corr = ", s08.cor), x = 11000, y = 0); rm(s08.cor)

### Threshold evaluation
s08.nm.thresholds <- filterCombinations(s08.nm, expression_thresholds = 1,
                                        feat_det_in_min_cells = c(50, 100, 50, 100),
                                        min_det_feats_per_cell = c(500, 500, 1000, 1000), 
                                        show_plot = F,
                                        default_save_name = "thresholds")
s08.nm.thresholds[["ggplot"]][["theme"]][["legend.position"]] <- "none" # Remove legend

s08.nm.metric.plots <- ggarrange(s08.nm.detected, s08.nm.det.along, 
                                 s08.nm.libvdet, s08.nm.thresholds$ggplot) # Combine all plots
rm(s08.nm.detected, s08.nm.det.along, s08.nm.libvdet, s08.nm.thresholds); s08.nm.metric.plots
ggsave("./project/outcomes/qc/ms02_metrics.png", plot = s08.nm.metric.plots, 
       scale = 2, width = 1920, height = 1080, units = "px")

## Filter sample
s08.filtered <- filterGiotto(s08.nm, expression_values = "raw", expression_threshold = 1,
                             feat_det_in_min_cells = 50,
                             min_det_feats_per_cell = 1000)
saveRDS(s08.filtered, file = "./project/material/filtered_samples/s08_filtered.rds") # readRDS()

## Deleted spots
length(s08@cell_ID$cell) - length(s08.filtered@cell_ID$cell) # 525 deleted spots

### Visualization plot
s08.spots.plot <- spatPlot2D(s08, cell_color = ("lightgrey"), point_size = 2,
                             select_cells = s08.filtered@cell_ID$cell, # Kept spots
                             other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                             title = "Deleted spots (sample 8)"); s08.spots.plot
ggsave("./project/outcomes/qc/ms02_deleted.png", plot = s08.spots.plot, 
       scale = 2.5, width = 1920, height = 1080, units = "px")

####--SAMPLE 09--####

# Obtain data
s09 <- createGiottoVisiumObject(
  visium_dir = sample.path[9], #ms03
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = instr)

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
s09.meta <- pDataDT(s09); rm(s09.matrix, is_mito, mito.sum, expr.sum, mito.perc)

### Keep low mito percentage cells
low.mito <- s09.meta$mito_perc <= 28
s09.nm <- s09[, low.mito]; rm(low.mito)

## Create cell metadata
s09.nm.matrix <- getExpression(s09.nm, values = "raw", output = "matrix") 
reads.depth <- apply(s09.nm.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s09.nm.matrix != 0, 2, sum) # Detected genes per cell
s09.nm.meta <- pDataDT(s09.nm)
s09.nm.meta$n_reads <- as.numeric(reads.depth); s09.nm.meta$n_genes <- as.numeric(genes.per.cell)
rm(s09.nm.matrix, reads.depth, genes.per.cell) 

### Libsize vs Detected genes correlation
s09.cor <- as.numeric(round(with(s09.nm.meta, cor.test(n_reads, n_genes))$estimate, 4)) 


## Plots
### Detected genes
s09.nm.detected <- filterDistributions(s09.nm, detection = "cells", nr_bins = 50,
                                       method = "threshold", expression_threshold = 1, 
                                       default_save_name = "detected_genes")

### Genes detection along cells
s09.nm.det.along <- filterDistributions(s09.nm, detection = "feats", nr_bins = 50,
                                        method = "threshold", 
                                        default_save_name = "library_size")

### Libsize vs Detected genes
s09.nm.libvdet <- ggplot(s09.nm.meta, aes(n_reads, n_genes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = paste0("Corr = ", s09.cor), x = 11000, y = 0); rm(s09.cor)

### Threshold evaluation
s09.nm.thresholds <- filterCombinations(s09.nm, expression_thresholds = 1,
                                        feat_det_in_min_cells = c(50, 100, 50, 100),
                                        min_det_feats_per_cell = c(500, 500, 1000, 1000), 
                                        show_plot = F,
                                        default_save_name = "thresholds")
s09.nm.thresholds[["ggplot"]][["theme"]][["legend.position"]] <- "none" # Remove legend

s09.nm.metric.plots <- ggarrange(s09.nm.detected, s09.nm.det.along, 
                                 s09.nm.libvdet, s09.nm.thresholds$ggplot) # Combine all plots
rm(s09.nm.detected, s09.nm.det.along, s09.nm.libvdet, s09.nm.thresholds); s09.nm.metric.plots
ggsave("./project/outcomes/qc/ms03_metrics.png", plot = s09.nm.metric.plots, 
       scale = 2, width = 1920, height = 1080, units = "px")

## Filter sample
s09.filtered <- filterGiotto(s09.nm, expression_values = "raw", expression_threshold = 1,
                             feat_det_in_min_cells = 50,
                             min_det_feats_per_cell = 1000)
saveRDS(s09.filtered, file = "./project/material/filtered_samples/s09_filtered.rds") # readRDS()

## Deleted spots
length(s09@cell_ID$cell) - length(s09.filtered@cell_ID$cell) # 525 deleted spots

### Visualization plot
s09.spots.plot <- spatPlot2D(s09, cell_color = ("lightgrey"), point_size = 2,
                             select_cells = s09.filtered@cell_ID$cell, # Kept spots
                             other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                             title = "Deleted spots (sample 9)"); s09.spots.plot
ggsave("./project/outcomes/qc/ms03_deleted.png", plot = s09.spots.plot, 
       scale = 2.5, width = 1920, height = 1080, units = "px")

####--SAMPLE 10--####

# Obtain data
s10 <- createGiottoVisiumObject(
  visium_dir = sample.path[10], #ms04
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = instr)

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
s10.meta <- pDataDT(s10); rm(s10.matrix, is_mito, mito.sum, expr.sum, mito.perc)

### Keep low mito percentage cells
low.mito <- s10.meta$mito_perc <= 28
s10.nm <- s10[, low.mito]; rm(low.mito)

## Create cell metadata
s10.nm.matrix <- getExpression(s10.nm, values = "raw", output = "matrix") 
reads.depth <- apply(s10.nm.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s10.nm.matrix != 0, 2, sum) # Detected genes per cell
s10.nm.meta <- pDataDT(s10.nm)
s10.nm.meta$n_reads <- as.numeric(reads.depth); s10.nm.meta$n_genes <- as.numeric(genes.per.cell)
rm(s10.nm.matrix, reads.depth, genes.per.cell) 

### Libsize vs Detected genes correlation
s10.cor <- as.numeric(round(with(s10.nm.meta, cor.test(n_reads, n_genes))$estimate, 4)) 


## Plots
### Detected genes
s10.nm.detected <- filterDistributions(s10.nm, detection = "cells", nr_bins = 50,
                                       method = "threshold", expression_threshold = 1, 
                                       default_save_name = "detected_genes")

### Genes detection along cells
s10.nm.det.along <- filterDistributions(s10.nm, detection = "feats", nr_bins = 50,
                                        method = "threshold", 
                                        default_save_name = "library_size")

### Libsize vs Detected genes
s10.nm.libvdet <- ggplot(s10.nm.meta, aes(n_reads, n_genes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = paste0("Corr = ", s10.cor), x = 11000, y = 0); rm(s10.cor)

### Threshold evaluation
s10.nm.thresholds <- filterCombinations(s10.nm, expression_thresholds = 1,
                                        feat_det_in_min_cells = c(50, 100, 50, 100),
                                        min_det_feats_per_cell = c(500, 500, 1000, 1000), 
                                        show_plot = F,
                                        default_save_name = "thresholds")
s10.nm.thresholds[["ggplot"]][["theme"]][["legend.position"]] <- "none" # Remove legend

s10.nm.metric.plots <- ggarrange(s10.nm.detected, s10.nm.det.along, 
                                 s10.nm.libvdet, s10.nm.thresholds$ggplot) # Combine all plots
rm(s10.nm.detected, s10.nm.det.along, s10.nm.libvdet, s10.nm.thresholds); s10.nm.metric.plots
ggsave("./project/outcomes/qc/ms04_metrics.png", plot = s10.nm.metric.plots, 
       scale = 2, width = 1920, height = 1080, units = "px")

## Filter sample
s10.filtered <- filterGiotto(s10.nm, expression_values = "raw", expression_threshold = 1,
                             feat_det_in_min_cells = 50,
                             min_det_feats_per_cell = 1000)
saveRDS(s10.filtered, file = "./project/material/filtered_samples/s10_filtered.rds") # readRDS()

## Deleted spots
length(s10@cell_ID$cell) - length(s10.filtered@cell_ID$cell) # 525 deleted spots

### Visualization plot
s10.spots.plot <- spatPlot2D(s10, cell_color = ("lightgrey"), point_size = 2,
                             select_cells = s10.filtered@cell_ID$cell, # Kept spots
                             other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                             title = "Deleted spots (sample 10)"); s10.spots.plot
ggsave("./project/outcomes/qc/ms04_deleted.png", plot = s10.spots.plot, 
       scale = 2.5, width = 1920, height = 1080, units = "px")

####--SAMPLE 11--####

# Obtain data
s11 <- createGiottoVisiumObject(
  visium_dir = sample.path[11], #ms05
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = instr)

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
s11.meta <- pDataDT(s11); rm(s11.matrix, is_mito, mito.sum, expr.sum, mito.perc)

### Keep low mito percentage cells
low.mito <- s11.meta$mito_perc <= 28
s11.nm <- s11[, low.mito]; rm(low.mito)

## Create cell metadata
s11.nm.matrix <- getExpression(s11.nm, values = "raw", output = "matrix") 
reads.depth <- apply(s11.nm.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s11.nm.matrix != 0, 2, sum) # Detected genes per cell
s11.nm.meta <- pDataDT(s11.nm)
s11.nm.meta$n_reads <- as.numeric(reads.depth); s11.nm.meta$n_genes <- as.numeric(genes.per.cell)
rm(s11.nm.matrix, reads.depth, genes.per.cell) 

### Libsize vs Detected genes correlation
s11.cor <- as.numeric(round(with(s11.nm.meta, cor.test(n_reads, n_genes))$estimate, 4)) 


## Plots
### Detected genes
s11.nm.detected <- filterDistributions(s11.nm, detection = "cells", nr_bins = 50,
                                       method = "threshold", expression_threshold = 1, 
                                       default_save_name = "detected_genes")

### Genes detection along cells
s11.nm.det.along <- filterDistributions(s11.nm, detection = "feats", nr_bins = 50,
                                        method = "threshold", 
                                        default_save_name = "library_size")

### Libsize vs Detected genes
s11.nm.libvdet <- ggplot(s11.nm.meta, aes(n_reads, n_genes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = paste0("Corr = ", s11.cor), x = 11000, y = 0); rm(s11.cor)

### Threshold evaluation
s11.nm.thresholds <- filterCombinations(s11.nm, expression_thresholds = 1,
                                        feat_det_in_min_cells = c(50, 100, 50, 100),
                                        min_det_feats_per_cell = c(500, 500, 1000, 1000), 
                                        show_plot = F,
                                        default_save_name = "thresholds")
s11.nm.thresholds[["ggplot"]][["theme"]][["legend.position"]] <- "none" # Remove legend

s11.nm.metric.plots <- ggarrange(s11.nm.detected, s11.nm.det.along, 
                                 s11.nm.libvdet, s11.nm.thresholds$ggplot) # Combine all plots
rm(s11.nm.detected, s11.nm.det.along, s11.nm.libvdet, s11.nm.thresholds); s11.nm.metric.plots
ggsave("./project/outcomes/qc/ms05_metrics.png", plot = s11.nm.metric.plots, 
       scale = 2, width = 1920, height = 1080, units = "px")

## Filter sample
s11.filtered <- filterGiotto(s11.nm, expression_values = "raw", expression_threshold = 1,
                             feat_det_in_min_cells = 50,
                             min_det_feats_per_cell = 1000)
saveRDS(s11.filtered, file = "./project/material/filtered_samples/s11_filtered.rds") # readRDS()

## Deleted spots
length(s11@cell_ID$cell) - length(s11.filtered@cell_ID$cell) # 525 deleted spots

### Visualization plot
s11.spots.plot <- spatPlot2D(s11, cell_color = ("lightgrey"), point_size = 2,
                             select_cells = s11.filtered@cell_ID$cell, # Kept spots
                             other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                             title = "Deleted spots (sample 11)"); s11.spots.plot
ggsave("./project/outcomes/qc/ms05_deleted.png", plot = s11.spots.plot, 
       scale = 2.5, width = 1920, height = 1080, units = "px")

####--SAMPLE 12--####

# Obtain data
s12 <- createGiottoVisiumObject(
  visium_dir = sample.path[12], #ms06
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = instr)

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
s12.meta <- pDataDT(s12); rm(s12.matrix, is_mito, mito.sum, expr.sum, mito.perc)

### Keep low mito percentage cells
low.mito <- s12.meta$mito_perc <= 28
s12.nm <- s12[, low.mito]; rm(low.mito)

## Create cell metadata
s12.nm.matrix <- getExpression(s12.nm, values = "raw", output = "matrix") 
reads.depth <- apply(s12.nm.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s12.nm.matrix != 0, 2, sum) # Detected genes per cell
s12.nm.meta <- pDataDT(s12.nm)
s12.nm.meta$n_reads <- as.numeric(reads.depth); s12.nm.meta$n_genes <- as.numeric(genes.per.cell)
rm(s12.nm.matrix, reads.depth, genes.per.cell) 

### Libsize vs Detected genes correlation
s12.cor <- as.numeric(round(with(s12.nm.meta, cor.test(n_reads, n_genes))$estimate, 4)) 


## Plots
### Detected genes
s12.nm.detected <- filterDistributions(s12.nm, detection = "cells", nr_bins = 50,
                                       method = "threshold", expression_threshold = 1, 
                                       default_save_name = "detected_genes")

### Genes detection along cells
s12.nm.det.along <- filterDistributions(s12.nm, detection = "feats", nr_bins = 50,
                                        method = "threshold", 
                                        default_save_name = "library_size")

### Libsize vs Detected genes
s12.nm.libvdet <- ggplot(s12.nm.meta, aes(n_reads, n_genes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = paste0("Corr = ", s12.cor), x = 11000, y = 0); rm(s12.cor)

### Threshold evaluation
s12.nm.thresholds <- filterCombinations(s12.nm, expression_thresholds = 1,
                                        feat_det_in_min_cells = c(50, 100, 50, 100),
                                        min_det_feats_per_cell = c(500, 500, 1000, 1000), 
                                        show_plot = F,
                                        default_save_name = "thresholds")
s12.nm.thresholds[["ggplot"]][["theme"]][["legend.position"]] <- "none" # Remove legend

s12.nm.metric.plots <- ggarrange(s12.nm.detected, s12.nm.det.along, 
                                 s12.nm.libvdet, s12.nm.thresholds$ggplot) # Combine all plots
rm(s12.nm.detected, s12.nm.det.along, s12.nm.libvdet, s12.nm.thresholds); s12.nm.metric.plots
ggsave("./project/outcomes/qc/ms06_metrics.png", plot = s12.nm.metric.plots, 
       scale = 2, width = 1920, height = 1080, units = "px")

## Filter sample
s12.filtered <- filterGiotto(s12.nm, expression_values = "raw", expression_threshold = 1,
                             feat_det_in_min_cells = 50,
                             min_det_feats_per_cell = 1000)
saveRDS(s12.filtered, file = "./project/material/filtered_samples/s12_filtered.rds") # readRDS()

## Deleted spots
length(s12@cell_ID$cell) - length(s12.filtered@cell_ID$cell) # 525 deleted spots

### Visualization plot
s12.spots.plot <- spatPlot2D(s12, cell_color = ("lightgrey"), point_size = 2,
                             select_cells = s12.filtered@cell_ID$cell, # Kept spots
                             other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                             title = "Deleted spots (sample 12)"); s12.spots.plot
ggsave("./project/outcomes/qc/ms06_deleted.png", plot = s12.spots.plot, 
       scale = 2.5, width = 1920, height = 1080, units = "px")

####--SAMPLE 13--####

# Obtain data
s13 <- createGiottoVisiumObject(
  visium_dir = sample.path[13], #ms07
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = instr)

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
s13.meta <- pDataDT(s13); rm(s13.matrix, is_mito, mito.sum, expr.sum, mito.perc)

### Keep low mito percentage cells
low.mito <- s13.meta$mito_perc <= 28
s13.nm <- s13[, low.mito]; rm(low.mito)

## Create cell metadata
s13.nm.matrix <- getExpression(s13.nm, values = "raw", output = "matrix") 
reads.depth <- apply(s13.nm.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s13.nm.matrix != 0, 2, sum) # Detected genes per cell
s13.nm.meta <- pDataDT(s13.nm)
s13.nm.meta$n_reads <- as.numeric(reads.depth); s13.nm.meta$n_genes <- as.numeric(genes.per.cell)
rm(s13.nm.matrix, reads.depth, genes.per.cell) 

### Libsize vs Detected genes correlation
s13.cor <- as.numeric(round(with(s13.nm.meta, cor.test(n_reads, n_genes))$estimate, 4)) 


## Plots
### Detected genes
s13.nm.detected <- filterDistributions(s13.nm, detection = "cells", nr_bins = 50,
                                       method = "threshold", expression_threshold = 1, 
                                       default_save_name = "detected_genes")

### Genes detection along cells
s13.nm.det.along <- filterDistributions(s13.nm, detection = "feats", nr_bins = 50,
                                        method = "threshold", 
                                        default_save_name = "library_size")

### Libsize vs Detected genes
s13.nm.libvdet <- ggplot(s13.nm.meta, aes(n_reads, n_genes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = paste0("Corr = ", s13.cor), x = 11000, y = 0); rm(s13.cor)

### Threshold evaluation
s13.nm.thresholds <- filterCombinations(s13.nm, expression_thresholds = 1,
                                        feat_det_in_min_cells = c(50, 100, 50, 100),
                                        min_det_feats_per_cell = c(500, 500, 1000, 1000), 
                                        show_plot = F,
                                        default_save_name = "thresholds")
s13.nm.thresholds[["ggplot"]][["theme"]][["legend.position"]] <- "none" # Remove legend

s13.nm.metric.plots <- ggarrange(s13.nm.detected, s13.nm.det.along, 
                                 s13.nm.libvdet, s13.nm.thresholds$ggplot) # Combine all plots
rm(s13.nm.detected, s13.nm.det.along, s13.nm.libvdet, s13.nm.thresholds); s13.nm.metric.plots
ggsave("./project/outcomes/qc/ms07_metrics.png", plot = s13.nm.metric.plots, 
       scale = 2, width = 1920, height = 1080, units = "px")

## Filter sample
s13.filtered <- filterGiotto(s13.nm, expression_values = "raw", expression_threshold = 1,
                             feat_det_in_min_cells = 50,
                             min_det_feats_per_cell = 1000)
saveRDS(s13.filtered, file = "./project/material/filtered_samples/s13_filtered.rds") # readRDS()

## Deleted spots
length(s13@cell_ID$cell) - length(s13.filtered@cell_ID$cell) # 525 deleted spots

### Visualization plot
s13.spots.plot <- spatPlot2D(s13, cell_color = ("lightgrey"), point_size = 2,
                             select_cells = s13.filtered@cell_ID$cell, # Kept spots
                             other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                             title = "Deleted spots (sample 13)"); s13.spots.plot
ggsave("./project/outcomes/qc/ms07_deleted.png", plot = s13.spots.plot, 
       scale = 2.5, width = 1920, height = 1080, units = "px")

####--SAMPLE 14--####

# Obtain data
s14 <- createGiottoVisiumObject(
  visium_dir = sample.path[14], #ms08
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = instr)

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
s14.meta <- pDataDT(s14); rm(s14.matrix, is_mito, mito.sum, expr.sum, mito.perc)

### Keep low mito percentage cells
low.mito <- s14.meta$mito_perc <= 28
s14.nm <- s14[, low.mito]; rm(low.mito)

## Create cell metadata
s14.nm.matrix <- getExpression(s14.nm, values = "raw", output = "matrix") 
reads.depth <- apply(s14.nm.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s14.nm.matrix != 0, 2, sum) # Detected genes per cell
s14.nm.meta <- pDataDT(s14.nm)
s14.nm.meta$n_reads <- as.numeric(reads.depth); s14.nm.meta$n_genes <- as.numeric(genes.per.cell)
rm(s14.nm.matrix, reads.depth, genes.per.cell) 

### Libsize vs Detected genes correlation
s14.cor <- as.numeric(round(with(s14.nm.meta, cor.test(n_reads, n_genes))$estimate, 4)) 


## Plots
### Detected genes
s14.nm.detected <- filterDistributions(s14.nm, detection = "cells", nr_bins = 50,
                                       method = "threshold", expression_threshold = 1, 
                                       default_save_name = "detected_genes")

### Genes detection along cells
s14.nm.det.along <- filterDistributions(s14.nm, detection = "feats", nr_bins = 50,
                                        method = "threshold", 
                                        default_save_name = "library_size")

### Libsize vs Detected genes
s14.nm.libvdet <- ggplot(s14.nm.meta, aes(n_reads, n_genes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = paste0("Corr = ", s14.cor), x = 11000, y = 0); rm(s14.cor)

### Threshold evaluation
s14.nm.thresholds <- filterCombinations(s14.nm, expression_thresholds = 1,
                                        feat_det_in_min_cells = c(50, 100, 50, 100),
                                        min_det_feats_per_cell = c(500, 500, 1000, 1000), 
                                        show_plot = F,
                                        default_save_name = "thresholds")
s14.nm.thresholds[["ggplot"]][["theme"]][["legend.position"]] <- "none" # Remove legend

s14.nm.metric.plots <- ggarrange(s14.nm.detected, s14.nm.det.along, 
                                 s14.nm.libvdet, s14.nm.thresholds$ggplot) # Combine all plots
rm(s14.nm.detected, s14.nm.det.along, s14.nm.libvdet, s14.nm.thresholds); s14.nm.metric.plots
ggsave("./project/outcomes/qc/ms08_metrics.png", plot = s14.nm.metric.plots, 
       scale = 2, width = 1920, height = 1080, units = "px")

## Filter sample
s14.filtered <- filterGiotto(s14.nm, expression_values = "raw", expression_threshold = 1,
                             feat_det_in_min_cells = 50,
                             min_det_feats_per_cell = 1000)
saveRDS(s14.filtered, file = "./project/material/filtered_samples/s14_filtered.rds") # readRDS()

## Deleted spots
length(s14@cell_ID$cell) - length(s14.filtered@cell_ID$cell) # 525 deleted spots

### Visualization plot
s14.spots.plot <- spatPlot2D(s14, cell_color = ("lightgrey"), point_size = 2,
                             select_cells = s14.filtered@cell_ID$cell, # Kept spots
                             other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                             title = "Deleted spots (sample 14)"); s14.spots.plot
ggsave("./project/outcomes/qc/ms08_deleted.png", plot = s14.spots.plot, 
       scale = 2.5, width = 1920, height = 1080, units = "px")

####--SAMPLE 15--####

# Obtain data
s15 <- createGiottoVisiumObject(
  visium_dir = sample.path[15], #ms09
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = instr)

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
s15.meta <- pDataDT(s15); rm(s15.matrix, is_mito, mito.sum, expr.sum, mito.perc)

### Keep low mito percentage cells
low.mito <- s15.meta$mito_perc <= 28
s15.nm <- s15[, low.mito]; rm(low.mito)

## Create cell metadata
s15.nm.matrix <- getExpression(s15.nm, values = "raw", output = "matrix") 
reads.depth <- apply(s15.nm.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s15.nm.matrix != 0, 2, sum) # Detected genes per cell
s15.nm.meta <- pDataDT(s15.nm)
s15.nm.meta$n_reads <- as.numeric(reads.depth); s15.nm.meta$n_genes <- as.numeric(genes.per.cell)
rm(s15.nm.matrix, reads.depth, genes.per.cell) 

### Libsize vs Detected genes correlation
s15.cor <- as.numeric(round(with(s15.nm.meta, cor.test(n_reads, n_genes))$estimate, 4)) 


## Plots
### Detected genes
s15.nm.detected <- filterDistributions(s15.nm, detection = "cells", nr_bins = 50,
                                       method = "threshold", expression_threshold = 1, 
                                       default_save_name = "detected_genes")

### Genes detection along cells
s15.nm.det.along <- filterDistributions(s15.nm, detection = "feats", nr_bins = 50,
                                        method = "threshold", 
                                        default_save_name = "library_size")

### Libsize vs Detected genes
s15.nm.libvdet <- ggplot(s15.nm.meta, aes(n_reads, n_genes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = paste0("Corr = ", s15.cor), x = 11000, y = 0); rm(s15.cor)

### Threshold evaluation
s15.nm.thresholds <- filterCombinations(s15.nm, expression_thresholds = 1,
                                        feat_det_in_min_cells = c(50, 100, 50, 100),
                                        min_det_feats_per_cell = c(500, 500, 1000, 1000), 
                                        show_plot = F,
                                        default_save_name = "thresholds")
s15.nm.thresholds[["ggplot"]][["theme"]][["legend.position"]] <- "none" # Remove legend

s15.nm.metric.plots <- ggarrange(s15.nm.detected, s15.nm.det.along, 
                                 s15.nm.libvdet, s15.nm.thresholds$ggplot) # Combine all plots
rm(s15.nm.detected, s15.nm.det.along, s15.nm.libvdet, s15.nm.thresholds); s15.nm.metric.plots
ggsave("./project/outcomes/qc/ms09_metrics.png", plot = s15.nm.metric.plots, 
       scale = 2, width = 1920, height = 1080, units = "px")

## Filter sample
s15.filtered <- filterGiotto(s15.nm, expression_values = "raw", expression_threshold = 1,
                             feat_det_in_min_cells = 50,
                             min_det_feats_per_cell = 1000)
saveRDS(s15.filtered, file = "./project/material/filtered_samples/s15_filtered.rds") # readRDS()

## Deleted spots
length(s15@cell_ID$cell) - length(s15.filtered@cell_ID$cell) # 525 deleted spots

### Visualization plot
s15.spots.plot <- spatPlot2D(s15, cell_color = ("lightgrey"), point_size = 2,
                             select_cells = s15.filtered@cell_ID$cell, # Kept spots
                             other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                             title = "Deleted spots (sample 15)"); s15.spots.plot
ggsave("./project/outcomes/qc/ms09_deleted.png", plot = s15.spots.plot, 
       scale = 2.5, width = 1920, height = 1080, units = "px")

####--SAMPLE 16--####

# Obtain data
s16 <- createGiottoVisiumObject(
  visium_dir = sample.path[16], #ms10
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = instr)

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
s16.meta <- pDataDT(s16); rm(s16.matrix, is_mito, mito.sum, expr.sum, mito.perc)

### Keep low mito percentage cells
low.mito <- s16.meta$mito_perc <= 28
s16.nm <- s16[, low.mito]; rm(low.mito)

## Create cell metadata
s16.nm.matrix <- getExpression(s16.nm, values = "raw", output = "matrix") 
reads.depth <- apply(s16.nm.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s16.nm.matrix != 0, 2, sum) # Detected genes per cell
s16.nm.meta <- pDataDT(s16.nm)
s16.nm.meta$n_reads <- as.numeric(reads.depth); s16.nm.meta$n_genes <- as.numeric(genes.per.cell)
rm(s16.nm.matrix, reads.depth, genes.per.cell) 

### Libsize vs Detected genes correlation
s16.cor <- as.numeric(round(with(s16.nm.meta, cor.test(n_reads, n_genes))$estimate, 4)) 


## Plots
### Detected genes
s16.nm.detected <- filterDistributions(s16.nm, detection = "cells", nr_bins = 50,
                                       method = "threshold", expression_threshold = 1, 
                                       default_save_name = "detected_genes")

### Genes detection along cells
s16.nm.det.along <- filterDistributions(s16.nm, detection = "feats", nr_bins = 50,
                                        method = "threshold", 
                                        default_save_name = "library_size")

### Libsize vs Detected genes
s16.nm.libvdet <- ggplot(s16.nm.meta, aes(n_reads, n_genes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = paste0("Corr = ", s16.cor), x = 11000, y = 0); rm(s16.cor)

### Threshold evaluation
s16.nm.thresholds <- filterCombinations(s16.nm, expression_thresholds = 1,
                                        feat_det_in_min_cells = c(50, 100, 50, 100),
                                        min_det_feats_per_cell = c(500, 500, 1000, 1000), 
                                        show_plot = F,
                                        default_save_name = "thresholds")
s16.nm.thresholds[["ggplot"]][["theme"]][["legend.position"]] <- "none" # Remove legend

s16.nm.metric.plots <- ggarrange(s16.nm.detected, s16.nm.det.along, 
                                 s16.nm.libvdet, s16.nm.thresholds$ggplot) # Combine all plots
rm(s16.nm.detected, s16.nm.det.along, s16.nm.libvdet, s16.nm.thresholds); s16.nm.metric.plots
ggsave("./project/outcomes/qc/ms10_metrics.png", plot = s16.nm.metric.plots, 
       scale = 2, width = 1920, height = 1080, units = "px")

## Filter sample
s16.filtered <- filterGiotto(s16.nm, expression_values = "raw", expression_threshold = 1,
                             feat_det_in_min_cells = 50,
                             min_det_feats_per_cell = 1000)
saveRDS(s16.filtered, file = "./project/material/filtered_samples/s16_filtered.rds") # readRDS()

## Deleted spots
length(s16@cell_ID$cell) - length(s16.filtered@cell_ID$cell) # 525 deleted spots

### Visualization plot
s16.spots.plot <- spatPlot2D(s16, cell_color = ("lightgrey"), point_size = 2,
                             select_cells = s16.filtered@cell_ID$cell, # Kept spots
                             other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                             title = "Deleted spots (sample 16)"); s16.spots.plot
ggsave("./project/outcomes/qc/ms10_deleted.png", plot = s16.spots.plot, 
       scale = 2.5, width = 1920, height = 1080, units = "px")

####--SAMPLE 17--####

# Obtain data
s17 <- createGiottoVisiumObject(
  visium_dir = sample.path[17], #ms11
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = instr)

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
s17.meta <- pDataDT(s17); rm(s17.matrix, is_mito, mito.sum, expr.sum, mito.perc)

### Keep low mito percentage cells
low.mito <- s17.meta$mito_perc <= 28
s17.nm <- s17[, low.mito]; rm(low.mito)

## Create cell metadata
s17.nm.matrix <- getExpression(s17.nm, values = "raw", output = "matrix") 
reads.depth <- apply(s17.nm.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s17.nm.matrix != 0, 2, sum) # Detected genes per cell
s17.nm.meta <- pDataDT(s17.nm)
s17.nm.meta$n_reads <- as.numeric(reads.depth); s17.nm.meta$n_genes <- as.numeric(genes.per.cell)
rm(s17.nm.matrix, reads.depth, genes.per.cell) 

### Libsize vs Detected genes correlation
s17.cor <- as.numeric(round(with(s17.nm.meta, cor.test(n_reads, n_genes))$estimate, 4)) 


## Plots
### Detected genes
s17.nm.detected <- filterDistributions(s17.nm, detection = "cells", nr_bins = 50,
                                       method = "threshold", expression_threshold = 1, 
                                       default_save_name = "detected_genes")

### Genes detection along cells
s17.nm.det.along <- filterDistributions(s17.nm, detection = "feats", nr_bins = 50,
                                        method = "threshold", 
                                        default_save_name = "library_size")

### Libsize vs Detected genes
s17.nm.libvdet <- ggplot(s17.nm.meta, aes(n_reads, n_genes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = paste0("Corr = ", s17.cor), x = 11000, y = 0); rm(s17.cor)

### Threshold evaluation
s17.nm.thresholds <- filterCombinations(s17.nm, expression_thresholds = 1,
                                        feat_det_in_min_cells = c(50, 100, 50, 100),
                                        min_det_feats_per_cell = c(500, 500, 1000, 1000), 
                                        show_plot = F,
                                        default_save_name = "thresholds")
s17.nm.thresholds[["ggplot"]][["theme"]][["legend.position"]] <- "none" # Remove legend

s17.nm.metric.plots <- ggarrange(s17.nm.detected, s17.nm.det.along, 
                                 s17.nm.libvdet, s17.nm.thresholds$ggplot) # Combine all plots
rm(s17.nm.detected, s17.nm.det.along, s17.nm.libvdet, s17.nm.thresholds); s17.nm.metric.plots
ggsave("./project/outcomes/qc/ms11_metrics.png", plot = s17.nm.metric.plots, 
       scale = 2, width = 1920, height = 1080, units = "px")

## Filter sample
s17.filtered <- filterGiotto(s17.nm, expression_values = "raw", expression_threshold = 1,
                             feat_det_in_min_cells = 50,
                             min_det_feats_per_cell = 1000)
saveRDS(s17.filtered, file = "./project/material/filtered_samples/s17_filtered.rds") # readRDS()

## Deleted spots
length(s17@cell_ID$cell) - length(s17.filtered@cell_ID$cell) # 525 deleted spots

### Visualization plot
s17.spots.plot <- spatPlot2D(s17, cell_color = ("lightgrey"), point_size = 2,
                             select_cells = s17.filtered@cell_ID$cell, # Kept spots
                             other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                             title = "Deleted spots (sample 17)"); s17.spots.plot
ggsave("./project/outcomes/qc/ms11_deleted.png", plot = s17.spots.plot, 
       scale = 2.5, width = 1920, height = 1080, units = "px")

####--SAMPLE 18--####

# Obtain data
s18 <- createGiottoVisiumObject(
  visium_dir = sample.path[18], #ms12
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = instr)

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
s18.meta <- pDataDT(s18); rm(s18.matrix, is_mito, mito.sum, expr.sum, mito.perc)

### Keep low mito percentage cells
low.mito <- s18.meta$mito_perc <= 28
s18.nm <- s18[, low.mito]; rm(low.mito)

## Create cell metadata
s18.nm.matrix <- getExpression(s18.nm, values = "raw", output = "matrix") 
reads.depth <- apply(s18.nm.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s18.nm.matrix != 0, 2, sum) # Detected genes per cell
s18.nm.meta <- pDataDT(s18.nm)
s18.nm.meta$n_reads <- as.numeric(reads.depth); s18.nm.meta$n_genes <- as.numeric(genes.per.cell)
rm(s18.nm.matrix, reads.depth, genes.per.cell) 

### Libsize vs Detected genes correlation
s18.cor <- as.numeric(round(with(s18.nm.meta, cor.test(n_reads, n_genes))$estimate, 4)) 


## Plots
### Detected genes
s18.nm.detected <- filterDistributions(s18.nm, detection = "cells", nr_bins = 50,
                                       method = "threshold", expression_threshold = 1, 
                                       default_save_name = "detected_genes")

### Genes detection along cells
s18.nm.det.along <- filterDistributions(s18.nm, detection = "feats", nr_bins = 50,
                                        method = "threshold", 
                                        default_save_name = "library_size")

### Libsize vs Detected genes
s18.nm.libvdet <- ggplot(s18.nm.meta, aes(n_reads, n_genes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = paste0("Corr = ", s18.cor), x = 11000, y = 0); rm(s18.cor)

### Threshold evaluation
s18.nm.thresholds <- filterCombinations(s18.nm, expression_thresholds = 1,
                                        feat_det_in_min_cells = c(50, 100, 50, 100),
                                        min_det_feats_per_cell = c(500, 500, 1000, 1000), 
                                        show_plot = F,
                                        default_save_name = "thresholds")
s18.nm.thresholds[["ggplot"]][["theme"]][["legend.position"]] <- "none" # Remove legend

s18.nm.metric.plots <- ggarrange(s18.nm.detected, s18.nm.det.along, 
                                 s18.nm.libvdet, s18.nm.thresholds$ggplot) # Combine all plots
rm(s18.nm.detected, s18.nm.det.along, s18.nm.libvdet, s18.nm.thresholds); s18.nm.metric.plots
ggsave("./project/outcomes/qc/ms12_metrics.png", plot = s18.nm.metric.plots, 
       scale = 2, width = 1920, height = 1080, units = "px")

## Filter sample
s18.filtered <- filterGiotto(s18.nm, expression_values = "raw", expression_threshold = 1,
                             feat_det_in_min_cells = 50,
                             min_det_feats_per_cell = 1000)
saveRDS(s18.filtered, file = "./project/material/filtered_samples/s18_filtered.rds") # readRDS()

## Deleted spots
length(s18@cell_ID$cell) - length(s18.filtered@cell_ID$cell) # 525 deleted spots

### Visualization plot
s18.spots.plot <- spatPlot2D(s18, cell_color = ("lightgrey"), point_size = 2,
                             select_cells = s18.filtered@cell_ID$cell, # Kept spots
                             other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                             title = "Deleted spots (sample 18)"); s18.spots.plot
ggsave("./project/outcomes/qc/ms12_deleted.png", plot = s18.spots.plot, 
       scale = 2.5, width = 1920, height = 1080, units = "px")
