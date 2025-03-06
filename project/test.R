####--LIBRARIES--####

library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(Giotto) #pak::pkg_install("drieslab/Giotto")
# library(spacexr) #devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
# library(nicheDE) #devtools::install_github("kaishumason/NicheDE")
# library(spatialGE) #devtools::install_github("fridleylab/spatialGE")


####--DIRECTORIES--####

dir <- "./project/material/GSE279181"
sample.ids <- c("ct01", "ct02", "ct03", "ct04", "ct05", "ct06", 
                "ms01", "ms02", "ms03", "ms04", "ms05", "ms06", "ms07", "ms08", "ms09", "ms10", "ms11", "ms12")
sample.path <- file.path(dir, sample.ids); rm(dir, sample.ids)


####--SAMPLE 01--####

# Obtain data
s01 <- createGiottoVisiumObject(
  visium_dir = sample.path[1], #ct01
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = createGiottoInstructions(save_plot = F, show_plot = F))

## Filter over tissue genes
in.tissue <- pDataDT(s01)[in_tissue == 1]$cell_ID 
if (sum((in.tissue == pDataDT(s01)$cell_ID) == F) == 0) {
  cat("All genes are over tissue")
  rm(in.tissue)
} else {
  s01 <- subsetGiotto(s01, cell_ids = in.tissue)
  rm(in.tissue)
} # All genes are over tissue

s01.matrix <- s01@expression[["cell"]][["rna"]][["raw"]]@exprMat # Extract expression matrix from Giotto object


# Quality control 
## Create cell metadata
reads.depth <- apply(s01.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s01.matrix != 0, 2, sum) # Detected genes per cell

is_mito <- grepl("(^MT-)|(^Mt-)|(^mt-)", s01@feat_ID$rna) # Identify mitochondrial genes
if (sum(is_mito == T) != 0) {
  s01@feat_ID$is_mito <- is_mito
  s01.mito <- subsetGiotto(s01, feat_ids = s01@feat_ID$rna[is_mito == T])
  cat(sum(is_mito == T), "mitochondrial genes found:", s01.mito@feat_ID$rna)
  rm(is_mito)
} else {
  cat('No mitochondrial genes found')
  rm(is_mito)
} # 13 mitochondrial genes found

expr.sum <- apply(s01.matrix, 2, sum) # Sum of gene expression values per cell
mito.sum <- apply(s01.matrix[s01.mito@feat_ID$rna, ], 2, sum) # Sum of mito gene expression values per cell
mito.prop <- round(mito.sum/expr.sum, 3) # Mito genes proportion

s01.cell.meta <- data.frame(nReads = reads.depth, nGenes = genes.per.cell, mitoProp = mito.prop) # Write data
rownames(s01.cell.meta) <- s01@cell_ID$cell # Add cell names
rm(reads.depth, genes.per.cell, expr.sum, mito.sum, mito.prop, s01.mito) # Rm variables


## Calcs
### Mito proportion threshold
which(s01.cell.meta$mitoProp > 0.28) # No cells exceed 0.28 mito genes proportion

## Libsize vs Detected genes correlation
with(s01.cell.meta, cor.test(nReads, nGenes)) # Check correlation (0.9850)


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
s01.libvdet <- ggplot(s01.cell.meta, aes(nReads, nGenes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = "Corr = 0.9584", x = 11000, y = -200)

### Threshold evaluation
s01.thresholds <- filterCombinations(s01, expression_thresholds = 1,
                                    feat_det_in_min_cells = c(50, 100, 50, 100),
                                    min_det_feats_per_cell = c(500, 500, 1000, 1000), 
                                    show_plot = F,
                                    default_save_name = "thresholds")
s01.thresholds[["ggplot"]][["theme"]][["legend.position"]] <- "none" # Remove legend

s01.metric.plots <- ggarrange(s01.detected, s01.det.along, s01.libvdet, s01.thresholds$ggplot) # Combine all plots
rm(s01.detected, s01.det.along, s01.libvdet, s01.thresholds); s01.metric.plots
ggsave("./project/outcomes/ct01_metrics.png", plot = s01.metric.plots, scale = 2, width = 1920, height = 1080, units = "px")

## Filter sample
s01.filtered <- filterGiotto(s01, expression_values = "raw", expression_threshold = 1,
                            feat_det_in_min_cells = 50,
                            min_det_feats_per_cell = 1000)
save(s01.filtered, file = "./project/material/filtered_samples/s01_filtered.R") # load("./project/filtered_samples/s01_filtered.R")

## Deleted spots
length(s01@cell_ID$cell) - length(s01.filtered@cell_ID$cell) # 240 deleted spots

### Visualization plot
s01.spots.plot <- spatPlot2D(s01, cell_color = ("lightgrey"), point_size = 2,
                            select_cells = s01.filtered@cell_ID$cell, # Kept spots
                            other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                            title = "Deleted spots (sample 1)"); s01.spots.plot
ggsave("./project/outcomes/ct01_deleted.png", plot = s01.spots.plot, scale = 2.5, width = 1920, height = 1080, units = "px")


####--SAMPLE 02--####

# Obtain data
s02 <- createGiottoVisiumObject(
  visium_dir = sample.path[2], #ct02
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = createGiottoInstructions(save_plot = F, show_plot = F))

## Filter over tissue genes
in.tissue <- pDataDT(s02)[in_tissue == 1]$cell_ID 
if (sum((in.tissue == pDataDT(s02)$cell_ID) == F) == 0) {
  cat("All genes are over tissue")
  rm(in.tissue)
} else {
  s02 <- subsetGiotto(s02, cell_ids = in.tissue)
  rm(in.tissue)
} # All genes are over tissue

s02.matrix <- s02@expression[["cell"]][["rna"]][["raw"]]@exprMat # Extract expression matrix from Giotto object


# Quality control 
## Create cell metadata
reads.depth <- apply(s02.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s02.matrix != 0, 2, sum) # Detected genes per cell

is_mito <- grepl("(^MT-)|(^Mt-)|(^mt-)", s02@feat_ID$rna) # Identify mitochondrial genes
if (sum(is_mito == T) != 0) {
  s02@feat_ID$is_mito <- is_mito
  s02.mito <- subsetGiotto(s02, feat_ids = s02@feat_ID$rna[is_mito == T])
  cat(sum(is_mito == T), "mitochondrial genes found:", s02.mito@feat_ID$rna)
  rm(is_mito)
} else {
  cat('No mitochondrial genes found')
  rm(is_mito)
} # 13 mitochondrial genes found

expr.sum <- apply(s02.matrix, 2, sum) # Sum of gene expression values per cell
mito.sum <- apply(s02.matrix[s02.mito@feat_ID$rna, ], 2, sum) # Sum of mito gene expression values per cell
mito.prop <- round(mito.sum/expr.sum, 3) # Mito genes proportion

s02.cell.meta <- data.frame(nReads = reads.depth, nGenes = genes.per.cell, mitoProp = mito.prop) # Write data
rownames(s02.cell.meta) <- s02@cell_ID$cell # Add cell names
rm(reads.depth, genes.per.cell, expr.sum, mito.sum, mito.prop, s02.mito) # Rm variables


## Calcs
### Mito proportion threshold
which(s02.cell.meta$mitoProp > 0.28) # 124 cells exceed 0.28 mito genes proportion

## Libsize vs Detected genes correlation
with(s02.cell.meta, cor.test(nReads, nGenes)) # Check correlation (0.9761)


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
s02.libvdet <- ggplot(s02.cell.meta, aes(nReads, nGenes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = "Corr = 0.9761", x = 19000, y = -200)

### Threshold evaluation
s02.thresholds <- filterCombinations(s02, expression_thresholds = 1,
                                    feat_det_in_min_cells = c(50, 100, 50, 100),
                                    min_det_feats_per_cell = c(500, 500, 1000, 1000), 
                                    show_plot = F,
                                    default_save_name = "thresholds")
s02.thresholds[["ggplot"]][["theme"]][["legend.position"]] <- "none" # Remove legend

s02.metric.plots <- ggarrange(s02.detected, s02.det.along, s02.libvdet, s02.thresholds$ggplot) # Combine all plots
rm(s02.detected, s02.det.along, s02.libvdet, s02.thresholds); s02.metric.plots
ggsave("./project/outcomes/ct02_metrics.png", plot = s02.metric.plots, scale = 2, width = 1920, height = 1080, units = "px")

## Filter sample
s02.filtered <- filterGiotto(s02, expression_values = "raw", expression_threshold = 1,
                            feat_det_in_min_cells = 50,
                            min_det_feats_per_cell = 750)
save(s02.filtered, file = "./project/material/filtered_samples/s02_filtered.R") # load("./project/filtered_samples/s02_filtered.R")

## Deleted spots
length(s02@cell_ID$cell) - length(s02.filtered@cell_ID$cell) # 277 deleted spots

### Visualization plot
s02.spots.plot <- spatPlot2D(s02, cell_color = ("lightgrey"), point_size = 2,
                            select_cells = s02.filtered@cell_ID$cell, # Kept spots
                            other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                            title = "Deleted spots (sample 2)"); s02.spots.plot
ggsave("./project/outcomes/ct02_deleted.png", plot = s02.spots.plot, scale = 2.5, width = 1920, height = 1080, units = "px")


####--SAMPLE 03--####

# Obtain data
s03 <- createGiottoVisiumObject(
  visium_dir = sample.path[3], #ct03
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = createGiottoInstructions(save_plot = F, show_plot = F))

## Filter over tissue genes
in.tissue <- pDataDT(s03)[in_tissue == 1]$cell_ID 
if (sum((in.tissue == pDataDT(s03)$cell_ID) == F) == 0) {
  cat("All genes are over tissue")
  rm(in.tissue)
} else {
  s03 <- subsetGiotto(s03, cell_ids = in.tissue)
  rm(in.tissue)
} # All genes are over tissue

s03.matrix <- s03@expression[["cell"]][["rna"]][["raw"]]@exprMat # Extract expression matrix from Giotto object


# Quality control 
## Create cell metadata
reads.depth <- apply(s03.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s03.matrix != 0, 2, sum) # Detected genes per cell

is_mito <- grepl("(^MT-)|(^Mt-)|(^mt-)", s03@feat_ID$rna) # Identify mitochondrial genes
if (sum(is_mito == T) != 0) {
  s03@feat_ID$is_mito <- is_mito
  s03.mito <- subsetGiotto(s03, feat_ids = s03@feat_ID$rna[is_mito == T])
  cat(sum(is_mito == T), "mitochondrial genes found:", s03.mito@feat_ID$rna)
  rm(is_mito)
} else {
  cat('No mitochondrial genes found')
  rm(is_mito)
} # 13 mitochondrial genes found

expr.sum <- apply(s03.matrix, 2, sum) # Sum of gene expression values per cell
mito.sum <- apply(s03.matrix[s03.mito@feat_ID$rna, ], 2, sum) # Sum of mito gene expression values per cell
mito.prop <- round(mito.sum/expr.sum, 3) # Mito genes proportion

s03.cell.meta <- data.frame(nReads = reads.depth, nGenes = genes.per.cell, mitoProp = mito.prop) # Write data
rownames(s03.cell.meta) <- s03@cell_ID$cell # Add cell names
rm(reads.depth, genes.per.cell, expr.sum, mito.sum, mito.prop, s03.mito) # Rm variables


## Calcs
### Mito proportion threshold
which(s03.cell.meta$mitoProp > 0.28) # 84 cells exceed 0.28 mito genes proportion

## Libsize vs Detected genes correlation
with(s03.cell.meta, cor.test(nReads, nGenes)) # Check correlation (0.9796)


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
s03.libvdet <- ggplot(s03.cell.meta, aes(nReads, nGenes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = "Corr = 0.9796", x = 18000, y = -200)

### Threshold evaluation
s03.thresholds <- filterCombinations(s03, expression_thresholds = 1,
                                    feat_det_in_min_cells = c(50, 100, 50, 100),
                                    min_det_feats_per_cell = c(500, 500, 1000, 1000), 
                                    show_plot = F,
                                    default_save_name = "thresholds")
s03.thresholds[["ggplot"]][["theme"]][["legend.position"]] <- "none" # Remove legend

s03.metric.plots <- ggarrange(s03.detected, s03.det.along, s03.libvdet, s03.thresholds$ggplot) # Combine all plots
rm(s03.detected, s03.det.along, s03.libvdet, s03.thresholds); s03.metric.plots
ggsave("./project/outcomes/ct03_metrics.png", plot = s03.metric.plots, scale = 2, width = 1920, height = 1080, units = "px")

## Filter sample
s03.filtered <- filterGiotto(s03, expression_values = "raw", expression_threshold = 1,
                            feat_det_in_min_cells = 50,
                            min_det_feats_per_cell = 500)
save(s03.filtered, file = "./project/material/filtered_samples/s03_filtered.R") # load("./project/filtered_samples/s03_filtered.R")

## Deleted spots
length(s03@cell_ID$cell) - length(s03.filtered@cell_ID$cell) # 340 deleted spots

### Visualization plot
s03.spots.plot <- spatPlot2D(s03, cell_color = ("lightgrey"), point_size = 2,
                            select_cells = s03.filtered@cell_ID$cell, # Kept spots
                            other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                            title = "Deleted spots (sample 3)"); s03.spots.plot
ggsave("./project/outcomes/ct03_deleted.png", plot = s03.spots.plot, scale = 2.5, width = 1920, height = 1080, units = "px")


####--SAMPLE 04--####

# Obtain data
s04 <- createGiottoVisiumObject(
  visium_dir = sample.path[4], #ct04
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = createGiottoInstructions(save_plot = F, show_plot = F))

## Filter over tissue genes
in.tissue <- pDataDT(s04)[in_tissue == 1]$cell_ID 
if (sum((in.tissue == pDataDT(s04)$cell_ID) == F) == 0) {
  cat("All genes are over tissue")
  rm(in.tissue)
} else {
  s04 <- subsetGiotto(s04, cell_ids = in.tissue)
  rm(in.tissue)
} # All genes are over tissue

s04.matrix <- s04@expression[["cell"]][["rna"]][["raw"]]@exprMat # Extract expression matrix from Giotto object


# Quality control 
## Create cell metadata
reads.depth <- apply(s04.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s04.matrix != 0, 2, sum) # Detected genes per cell

is_mito <- grepl("(^MT-)|(^Mt-)|(^mt-)", s04@feat_ID$rna) # Identify mitochondrial genes
if (sum(is_mito == T) != 0) {
  s04@feat_ID$is_mito <- is_mito
  s04.mito <- subsetGiotto(s04, feat_ids = s04@feat_ID$rna[is_mito == T])
  cat(sum(is_mito == T), "mitochondrial genes found:", s04.mito@feat_ID$rna)
  rm(is_mito)
} else {
  cat('No mitochondrial genes found')
  rm(is_mito)
} # 13 mitochondrial genes found

expr.sum <- apply(s04.matrix, 2, sum) # Sum of gene expression values per cell
mito.sum <- apply(s04.matrix[s04.mito@feat_ID$rna, ], 2, sum) # Sum of mito gene expression values per cell
mito.prop <- round(mito.sum/expr.sum, 3) # Mito genes proportion

s04.cell.meta <- data.frame(nReads = reads.depth, nGenes = genes.per.cell, mitoProp = mito.prop) # Write data
rownames(s04.cell.meta) <- s04@cell_ID$cell # Add cell names
rm(reads.depth, genes.per.cell, expr.sum, mito.sum, mito.prop, s04.mito) # Rm variables


## Calcs
### Mito proportion threshold
which(s04.cell.meta$mitoProp > 0.28) # 30 cells exceed 0.28 mito genes proportion

## Libsize vs Detected genes correlation
with(s04.cell.meta, cor.test(nReads, nGenes)) # Check correlation (0.9760)


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
s04.libvdet <- ggplot(s04.cell.meta, aes(nReads, nGenes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = "Corr = 0.9760", x = 26000, y = -200)

### Threshold evaluation
s04.thresholds <- filterCombinations(s04, expression_thresholds = 1,
                                    feat_det_in_min_cells = c(50, 100, 50, 100),
                                    min_det_feats_per_cell = c(500, 500, 1000, 1000), 
                                    show_plot = F,
                                    default_save_name = "thresholds")
s04.thresholds[["ggplot"]][["theme"]][["legend.position"]] <- "none" # Remove legend

s04.metric.plots <- ggarrange(s04.detected, s04.det.along, s04.libvdet, s04.thresholds$ggplot) # Combine all plots
rm(s04.detected, s04.det.along, s04.libvdet, s04.thresholds); s04.metric.plots
ggsave("./project/outcomes/ct04_metrics.png", plot = s04.metric.plots, scale = 2, width = 1920, height = 1080, units = "px")

## Filter sample
s04.filtered <- filterGiotto(s04, expression_values = "raw", expression_threshold = 1,
                            feat_det_in_min_cells = 50,
                            min_det_feats_per_cell = 1000)
save(s04.filtered, file = "./project/material/filtered_samples/s04_filtered.R") # load("./project/filtered_samples/s04_filtered.R")

## Deleted spots
length(s04@cell_ID$cell) - length(s04.filtered@cell_ID$cell) # 240 deleted spots

### Visualization plot
s04.spots.plot <- spatPlot2D(s04, cell_color = ("lightgrey"), point_size = 2,
                            select_cells = s04.filtered@cell_ID$cell, # Kept spots
                            other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                            title = "Deleted spots (sample 4)"); s04.spots.plot
ggsave("./project/outcomes/ct04_deleted.png", plot = s04.spots.plot, scale = 2.5, width = 1920, height = 1080, units = "px")


####--SAMPLE 05--####

# Obtain data
s05 <- createGiottoVisiumObject(
  visium_dir = sample.path[5], #ct05
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = createGiottoInstructions(save_plot = F, show_plot = F))

## Filter over tissue genes
in.tissue <- pDataDT(s05)[in_tissue == 1]$cell_ID 
if (sum((in.tissue == pDataDT(s05)$cell_ID) == F) == 0) {
  cat("All genes are over tissue")
  rm(in.tissue)
} else {
  s05 <- subsetGiotto(s05, cell_ids = in.tissue)
  rm(in.tissue)
} # All genes are over tissue

s05.matrix <- s05@expression[["cell"]][["rna"]][["raw"]]@exprMat # Extract expression matrix from Giotto object


# Quality control 
## Create cell metadata
reads.depth <- apply(s05.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s05.matrix != 0, 2, sum) # Detected genes per cell

is_mito <- grepl("(^MT-)|(^Mt-)|(^mt-)", s05@feat_ID$rna) # Identify mitochondrial genes
if (sum(is_mito == T) != 0) {
  s05@feat_ID$is_mito <- is_mito
  s05.mito <- subsetGiotto(s05, feat_ids = s05@feat_ID$rna[is_mito == T])
  cat(sum(is_mito == T), "mitochondrial genes found:", s05.mito@feat_ID$rna)
  rm(is_mito)
} else {
  cat('No mitochondrial genes found')
  rm(is_mito)
} # 13 mitochondrial genes found

expr.sum <- apply(s05.matrix, 2, sum) # Sum of gene expression values per cell
mito.sum <- apply(s05.matrix[s05.mito@feat_ID$rna, ], 2, sum) # Sum of mito gene expression values per cell
mito.prop <- round(mito.sum/expr.sum, 3) # Mito genes proportion

s05.cell.meta <- data.frame(nReads = reads.depth, nGenes = genes.per.cell, mitoProp = mito.prop) # Write data
rownames(s05.cell.meta) <- s05@cell_ID$cell # Add cell names
rm(reads.depth, genes.per.cell, expr.sum, mito.sum, mito.prop, s05.mito) # Rm variables


## Calcs
### Mito proportion threshold
which(s05.cell.meta$mitoProp > 0.28) # No cells exceed 0.28 mito genes proportion

## Libsize vs Detected genes correlation
with(s05.cell.meta, cor.test(nReads, nGenes)) # Check correlation (0.9864)


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
s05.libvdet <- ggplot(s05.cell.meta, aes(nReads, nGenes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = "Corr = 0.9864", x = 13000, y = -200)

### Threshold evaluation
s05.thresholds <- filterCombinations(s05, expression_thresholds = 1,
                                    feat_det_in_min_cells = c(50, 100, 50, 100),
                                    min_det_feats_per_cell = c(500, 500, 1000, 1000), 
                                    show_plot = F,
                                    default_save_name = "thresholds")
s05.thresholds[["ggplot"]][["theme"]][["legend.position"]] <- "none" # Remove legend

s05.metric.plots <- ggarrange(s05.detected, s05.det.along, s05.libvdet, s05.thresholds$ggplot) # Combine all plots
rm(s05.detected, s05.det.along, s05.libvdet, s05.thresholds); s05.metric.plots
ggsave("./project/outcomes/ct05_metrics.png", plot = s05.metric.plots, scale = 2, width = 1920, height = 1080, units = "px")

## Filter sample
s05.filtered <- filterGiotto(s05, expression_values = "raw", expression_threshold = 1,
                            feat_det_in_min_cells = 50,
                            min_det_feats_per_cell = 500)
save(s05.filtered, file = "./project/material/filtered_samples/s05_filtered.R") # load("./project/filtered_samples/s05_filtered.R")

## Deleted spots
length(s05@cell_ID$cell) - length(s05.filtered@cell_ID$cell) # 236 deleted spots

### Visualization plot
s05.spots.plot <- spatPlot2D(s05, cell_color = ("lightgrey"), point_size = 2,
                            select_cells = s05.filtered@cell_ID$cell, # Kept spots
                            other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                            title = "Deleted spots (sample 5)"); s05.spots.plot
ggsave("./project/outcomes/ct05_deleted.png", plot = s05.spots.plot, scale = 2.5, width = 1920, height = 1080, units = "px")


####--SAMPLE 06--####

# Obtain data
s06 <- createGiottoVisiumObject(
  visium_dir = sample.path[6], #ct06
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = createGiottoInstructions(save_plot = F, show_plot = F))

## Filter over tissue genes
in.tissue <- pDataDT(s06)[in_tissue == 1]$cell_ID 
if (sum((in.tissue == pDataDT(s06)$cell_ID) == F) == 0) {
  cat("All genes are over tissue")
  rm(in.tissue)
} else {
  s06 <- subsetGiotto(s06, cell_ids = in.tissue)
  rm(in.tissue)
} # All genes are over tissue

s06.matrix <- s06@expression[["cell"]][["rna"]][["raw"]]@exprMat # Extract expression matrix from Giotto object


# Quality control 
## Create cell metadata
reads.depth <- apply(s06.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s06.matrix != 0, 2, sum) # Detected genes per cell

is_mito <- grepl("(^MT-)|(^Mt-)|(^mt-)", s06@feat_ID$rna) # Identify mitochondrial genes
if (sum(is_mito == T) != 0) {
  s06@feat_ID$is_mito <- is_mito
  s06.mito <- subsetGiotto(s06, feat_ids = s06@feat_ID$rna[is_mito == T])
  cat(sum(is_mito == T), "mitochondrial genes found:", s06.mito@feat_ID$rna)
  rm(is_mito)
} else {
  cat('No mitochondrial genes found')
  rm(is_mito)
} # 13 mitochondrial genes found

expr.sum <- apply(s06.matrix, 2, sum) # Sum of gene expression values per cell
mito.sum <- apply(s06.matrix[s06.mito@feat_ID$rna, ], 2, sum) # Sum of mito gene expression values per cell
mito.prop <- round(mito.sum/expr.sum, 3) # Mito genes proportion

s06.cell.meta <- data.frame(nReads = reads.depth, nGenes = genes.per.cell, mitoProp = mito.prop) # Write data
rownames(s06.cell.meta) <- s06@cell_ID$cell # Add cell names
rm(reads.depth, genes.per.cell, expr.sum, mito.sum, mito.prop, s06.mito) # Rm variables


## Calcs
### Mito proportion threshold
which(s06.cell.meta$mitoProp > 0.28) # 19 cells exceed 0.28 mito genes proportion

## Libsize vs Detected genes correlation
with(s06.cell.meta, cor.test(nReads, nGenes)) # Check correlation (0.9912)


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
s06.libvdet <- ggplot(s06.cell.meta, aes(nReads, nGenes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = "Corr = 0.9912", x = 7000, y = -200)

### Threshold evaluation
s06.thresholds <- filterCombinations(s06, expression_thresholds = 1,
                                    feat_det_in_min_cells = c(50, 100, 50, 100),
                                    min_det_feats_per_cell = c(125, 125, 250, 250), 
                                    show_plot = F,
                                    default_save_name = "thresholds")
s06.thresholds[["ggplot"]][["theme"]][["legend.position"]] <- "none" # Remove legend

s06.metric.plots <- ggarrange(s06.detected, s06.det.along, s06.libvdet, s06.thresholds$ggplot) # Combine all plots
rm(s06.detected, s06.det.along, s06.libvdet, s06.thresholds); s06.metric.plots
ggsave("./project/outcomes/ct06_metrics.png", plot = s06.metric.plots, scale = 2, width = 1920, height = 1080, units = "px")

## Filter sample
s06.filtered <- filterGiotto(s06, expression_values = "raw", expression_threshold = 1,
                            feat_det_in_min_cells = 50,
                            min_det_feats_per_cell = 125)
save(s06.filtered, file = "./project/material/filtered_samples/s06_filtered.R") # load("./project/filtered_samples/s06_filtered.R")

## Deleted spots
length(s06@cell_ID$cell) - length(s06.filtered@cell_ID$cell) # 74 deleted spots

### Visualization plot
s06.spots.plot <- spatPlot2D(s06, cell_color = ("lightgrey"), point_size = 2,
                            select_cells = s06.filtered@cell_ID$cell, # Kept spots
                            other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                            title = "Deleted spots (sample 6)"); s06.spots.plot
ggsave("./project/outcomes/ct06_deleted.png", plot = s06.spots.plot, scale = 2.5, width = 1920, height = 1080, units = "px")


####--SAMPLE 07 (REVISAR!!!!)--####

# Obtain data
s07 <- createGiottoVisiumObject(
  visium_dir = sample.path[7], #ms01
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = createGiottoInstructions(save_plot = F, show_plot = F))

## Filter over tissue genes
in.tissue <- pDataDT(s07)[in_tissue == 1]$cell_ID 
if (sum((in.tissue == pDataDT(s07)$cell_ID) == F) == 0) {
  cat("All genes are over tissue")
  rm(in.tissue)
} else {
  s07 <- subsetGiotto(s07, cell_ids = in.tissue)
  rm(in.tissue)
} # All genes are over tissue

s07.matrix <- s07@expression[["cell"]][["rna"]][["raw"]]@exprMat # Extract expression matrix from Giotto object


# Quality control 
## Create cell metadata
reads.depth <- apply(s07.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s07.matrix != 0, 2, sum) # Detected genes per cell

is_mito <- grepl("(^MT-)|(^Mt-)|(^mt-)", s07@feat_ID$rna) # Identify mitochondrial genes
if (sum(is_mito == T) != 0) {
  s07@feat_ID$is_mito <- is_mito
  s07.mito <- subsetGiotto(s07, feat_ids = s07@feat_ID$rna[is_mito == T])
  cat(sum(is_mito == T), "mitochondrial genes found:", s07.mito@feat_ID$rna)
  rm(is_mito)
} else {
  cat('No mitochondrial genes found')
  rm(is_mito)
} # 13 mitochondrial genes found

expr.sum <- apply(s07.matrix, 2, sum) # Sum of gene expression values per cell
mito.sum <- apply(s07.matrix[s07.mito@feat_ID$rna, ], 2, sum) # Sum of mito gene expression values per cell
mito.prop <- round(mito.sum/expr.sum, 3) # Mito genes proportion

s07.cell.meta <- data.frame(nReads = reads.depth, nGenes = genes.per.cell, mitoProp = mito.prop) # Write data
rownames(s07.cell.meta) <- s07@cell_ID$cell # Add cell names
rm(reads.depth, genes.per.cell, expr.sum, mito.sum, mito.prop, s07.mito) # Rm variables


## Calcs
### Mito proportion threshold
which(s07.cell.meta$mitoProp > 0.28) # 3824 cells exceed 0.28 mito genes proportion

## Libsize vs Detected genes correlation
with(s07.cell.meta, cor.test(nReads, nGenes)) # Check correlation (0.8293)


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
s07.libvdet <- ggplot(s07.cell.meta, aes(nReads, nGenes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = "Corr = 0.8293", x = 12000, y = -200)

### Threshold evaluation
s07.thresholds <- filterCombinations(s07, expression_thresholds = 1,
                                    feat_det_in_min_cells = c(50, 100, 50, 100),
                                    min_det_feats_per_cell = c(250, 250, 500, 500), 
                                    show_plot = F,
                                    default_save_name = "thresholds")
s07.thresholds[["ggplot"]][["theme"]][["legend.position"]] <- "none" # Remove legend

s07.metric.plots <- ggarrange(s07.detected, s07.det.along, s07.libvdet, s07.thresholds$ggplot) # Combine all plots
rm(s07.detected, s07.det.along, s07.libvdet, s07.thresholds); s07.metric.plots
ggsave("./project/outcomes/ms01_metrics.png", plot = s07.metric.plots, scale = 2, width = 1920, height = 1080, units = "px")

## Filter sample
s07.filtered <- filterGiotto(s07, expression_values = "raw", expression_threshold = 1,
                            feat_det_in_min_cells = 50,
                            min_det_feats_per_cell = 250)
save(s07.filtered, file = "./project/material/filtered_samples/s07_filtered.R") # load("./project/filtered_samples/s07_filtered.R")

## Deleted spots
length(s07@cell_ID$cell) - length(s07.filtered@cell_ID$cell) # 431 deleted spots

### Visualization plot
s07.spots.plot <- spatPlot2D(s07, cell_color = ("lightgrey"), point_size = 2,
                            select_cells = s07.filtered@cell_ID$cell, # Kept spots
                            other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                            title = "Deleted spots (sample 7)"); s07.spots.plot
ggsave("./project/outcomes/ms01_deleted.png", plot = s07.spots.plot, scale = 2.5, width = 1920, height = 1080, units = "px")


####--SAMPLE 08--####

# Obtain data
s08 <- createGiottoVisiumObject(
  visium_dir = sample.path[8], #ms02
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = createGiottoInstructions(save_plot = F, show_plot = F))

## Filter over tissue genes
in.tissue <- pDataDT(s08)[in_tissue == 1]$cell_ID 
if (sum((in.tissue == pDataDT(s08)$cell_ID) == F) == 0) {
  cat("All genes are over tissue")
  rm(in.tissue)
} else {
  s08 <- subsetGiotto(s08, cell_ids = in.tissue)
  rm(in.tissue)
} # All genes are over tissue

s08.matrix <- s08@expression[["cell"]][["rna"]][["raw"]]@exprMat # Extract expression matrix from Giotto object


# Quality control 
## Create cell metadata
reads.depth <- apply(s08.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s08.matrix != 0, 2, sum) # Detected genes per cell

is_mito <- grepl("(^MT-)|(^Mt-)|(^mt-)", s08@feat_ID$rna) # Identify mitochondrial genes
if (sum(is_mito == T) != 0) {
  s08@feat_ID$is_mito <- is_mito
  s08.mito <- subsetGiotto(s08, feat_ids = s08@feat_ID$rna[is_mito == T])
  cat(sum(is_mito == T), "mitochondrial genes found:", s08.mito@feat_ID$rna)
  rm(is_mito)
} else {
  cat('No mitochondrial genes found')
  rm(is_mito)
} # 13 mitochondrial genes found

expr.sum <- apply(s08.matrix, 2, sum) # Sum of gene expression values per cell
mito.sum <- apply(s08.matrix[s08.mito@feat_ID$rna, ], 2, sum) # Sum of mito gene expression values per cell
mito.prop <- round(mito.sum/expr.sum, 3) # Mito genes proportion

s08.cell.meta <- data.frame(nReads = reads.depth, nGenes = genes.per.cell, mitoProp = mito.prop) # Write data
rownames(s08.cell.meta) <- s08@cell_ID$cell # Add cell names
rm(reads.depth, genes.per.cell, expr.sum, mito.sum, mito.prop, s08.mito) # Rm variables


## Calcs
### Mito proportion threshold
which(s08.cell.meta$mitoProp > 0.28) # 4 cells exceed 0.28 mito genes proportion

## Libsize vs Detected genes correlation
with(s08.cell.meta, cor.test(nReads, nGenes)) # Check correlation (0.9827)


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
s08.libvdet <- ggplot(s08.cell.meta, aes(nReads, nGenes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = "Corr = 0.9827", x = 22000, y = -200)

### Threshold evaluation
s08.thresholds <- filterCombinations(s08, expression_thresholds = 1,
                                    feat_det_in_min_cells = c(50, 100, 50, 100),
                                    min_det_feats_per_cell = c(125, 125, 250, 250), 
                                    show_plot = F,
                                    default_save_name = "thresholds")
s08.thresholds[["ggplot"]][["theme"]][["legend.position"]] <- "none" # Remove legend

s08.metric.plots <- ggarrange(s08.detected, s08.det.along, s08.libvdet, s08.thresholds$ggplot) # Combine all plots
rm(s08.detected, s08.det.along, s08.libvdet, s08.thresholds); s08.metric.plots
ggsave("./project/outcomes/ms02_metrics.png", plot = s08.metric.plots, scale = 2, width = 1920, height = 1080, units = "px")

## Filter sample
s08.filtered <- filterGiotto(s08, expression_values = "raw", expression_threshold = 1,
                            feat_det_in_min_cells = 50,
                            min_det_feats_per_cell = 150)
save(s08.filtered, file = "./project/material/filtered_samples/s08_filtered.R") # load("./project/filtered_samples/s08_filtered.R")

## Deleted spots
length(s08@cell_ID$cell) - length(s08.filtered@cell_ID$cell) # 86 deleted spots

### Visualization plot
s08.spots.plot <- spatPlot2D(s08, cell_color = ("lightgrey"), point_size = 2,
                            select_cells = s08.filtered@cell_ID$cell, # Kept spots
                            other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                            title = "Deleted spots (sample 8)"); s08.spots.plot
ggsave("./project/outcomes/ms02_deleted.png", plot = s08.spots.plot, scale = 2.5, width = 1920, height = 1080, units = "px")


####--SAMPLE 09--####

# Obtain data
s09 <- createGiottoVisiumObject(
  visium_dir = sample.path[9], #ms03
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = createGiottoInstructions(save_plot = F, show_plot = F))

## Filter over tissue genes
in.tissue <- pDataDT(s09)[in_tissue == 1]$cell_ID 
if (sum((in.tissue == pDataDT(s09)$cell_ID) == F) == 0) {
  cat("All genes are over tissue")
  rm(in.tissue)
} else {
  s09 <- subsetGiotto(s09, cell_ids = in.tissue)
  rm(in.tissue)
} # All genes are over tissue

s09.matrix <- s09@expression[["cell"]][["rna"]][["raw"]]@exprMat # Extract expression matrix from Giotto object


# Quality control 
## Create cell metadata
reads.depth <- apply(s09.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s09.matrix != 0, 2, sum) # Detected genes per cell

is_mito <- grepl("(^MT-)|(^Mt-)|(^mt-)", s09@feat_ID$rna) # Identify mitochondrial genes
if (sum(is_mito == T) != 0) {
  s09@feat_ID$is_mito <- is_mito
  s09.mito <- subsetGiotto(s09, feat_ids = s09@feat_ID$rna[is_mito == T])
  cat(sum(is_mito == T), "mitochondrial genes found:", s09.mito@feat_ID$rna)
  rm(is_mito)
} else {
  cat('No mitochondrial genes found')
  rm(is_mito)
} # 13 mitochondrial genes found

expr.sum <- apply(s09.matrix, 2, sum) # Sum of gene expression values per cell
mito.sum <- apply(s09.matrix[s09.mito@feat_ID$rna, ], 2, sum) # Sum of mito gene expression values per cell
mito.prop <- round(mito.sum/expr.sum, 3) # Mito genes proportion

s09.cell.meta <- data.frame(nReads = reads.depth, nGenes = genes.per.cell, mitoProp = mito.prop) # Write data
rownames(s09.cell.meta) <- s09@cell_ID$cell # Add cell names
rm(reads.depth, genes.per.cell, expr.sum, mito.sum, mito.prop, s09.mito) # Rm variables


## Calcs
### Mito proportion threshold
which(s09.cell.meta$mitoProp > 0.28) # 3 cells exceed 0.28 mito genes proportion

## Libsize vs Detected genes correlation
with(s09.cell.meta, cor.test(nReads, nGenes)) # Check correlation (0.9859)


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
s09.libvdet <- ggplot(s09.cell.meta, aes(nReads, nGenes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = "Corr = 0.9859", x = 18000, y = -200)

### Threshold evaluation
s09.thresholds <- filterCombinations(s09, expression_thresholds = 1,
                                    feat_det_in_min_cells = c(50, 100, 50, 100),
                                    min_det_feats_per_cell = c(250, 250, 500, 500), 
                                    show_plot = F,
                                    default_save_name = "thresholds")
s09.thresholds[["ggplot"]][["theme"]][["legend.position"]] <- "none" # Remove legend

s09.metric.plots <- ggarrange(s09.detected, s09.det.along, s09.libvdet, s09.thresholds$ggplot) # Combine all plots
rm(s09.detected, s09.det.along, s09.libvdet, s09.thresholds); s09.metric.plots
ggsave("./project/outcomes/ms03_metrics.png", plot = s09.metric.plots, scale = 2, width = 1920, height = 1080, units = "px")

## Filter sample
s09.filtered <- filterGiotto(s09, expression_values = "raw", expression_threshold = 1,
                            feat_det_in_min_cells = 50,
                            min_det_feats_per_cell = 250)
save(s09.filtered, file = "./project/material/filtered_samples/s09_filtered.R") # load("./project/filtered_samples/s09_filtered.R")

## Deleted spots
length(s09@cell_ID$cell) - length(s09.filtered@cell_ID$cell) # 88 deleted spots

### Visualization plot
s09.spots.plot <- spatPlot2D(s09, cell_color = ("lightgrey"), point_size = 2,
                            select_cells = s09.filtered@cell_ID$cell, # Kept spots
                            other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                            title = "Deleted spots (sample 9)"); s09.spots.plot
ggsave("./project/outcomes/ms03_deleted.png", plot = s09.spots.plot, scale = 2.5, width = 1920, height = 1080, units = "px")


####--SAMPLE 10 (REVISAR)--####

# Obtain data
s10 <- createGiottoVisiumObject(
  visium_dir = sample.path[10], #ms04
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = createGiottoInstructions(save_plot = F, show_plot = F))

## Filter over tissue genes
in.tissue <- pDataDT(s10)[in_tissue == 1]$cell_ID 
if (sum((in.tissue == pDataDT(s10)$cell_ID) == F) == 0) {
  cat("All genes are over tissue")
  rm(in.tissue)
} else {
  s10 <- subsetGiotto(s10, cell_ids = in.tissue)
  rm(in.tissue)
} # All genes are over tissue

s10.matrix <- s10@expression[["cell"]][["rna"]][["raw"]]@exprMat # Extract expression matrix from Giotto object


# Quality control 
## Create cell metadata
reads.depth <- apply(s10.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s10.matrix != 0, 2, sum) # Detected genes per cell

is_mito <- grepl("(^MT-)|(^Mt-)|(^mt-)", s10@feat_ID$rna) # Identify mitochondrial genes
if (sum(is_mito == T) != 0) {
  s10@feat_ID$is_mito <- is_mito
  s10.mito <- subsetGiotto(s10, feat_ids = s10@feat_ID$rna[is_mito == T])
  cat(sum(is_mito == T), "mitochondrial genes found:", s10.mito@feat_ID$rna)
  rm(is_mito)
} else {
  cat('No mitochondrial genes found')
  rm(is_mito)
} # 13 mitochondrial genes found

expr.sum <- apply(s10.matrix, 2, sum) # Sum of gene expression values per cell
mito.sum <- apply(s10.matrix[s10.mito@feat_ID$rna, ], 2, sum) # Sum of mito gene expression values per cell
mito.prop <- round(mito.sum/expr.sum, 3) # Mito genes proportion

s10.cell.meta <- data.frame(nReads = reads.depth, nGenes = genes.per.cell, mitoProp = mito.prop) # Write data
rownames(s10.cell.meta) <- s10@cell_ID$cell # Add cell names
rm(reads.depth, genes.per.cell, expr.sum, mito.sum, mito.prop, s10.mito) # Rm variables


## Calcs
### Mito proportion threshold
which(s10.cell.meta$mitoProp > 0.28) # 1550 cells exceed 0.28 mito genes proportion

## Libsize vs Detected genes correlation
with(s10.cell.meta, cor.test(nReads, nGenes)) # Check correlation (0.9721)


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
s10.libvdet <- ggplot(s10.cell.meta, aes(nReads, nGenes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = "Corr = 0.9721", x = 14000, y = -200)

### Threshold evaluation
s10.thresholds <- filterCombinations(s10, expression_thresholds = 1,
                                    feat_det_in_min_cells = c(50, 100, 50, 100),
                                    min_det_feats_per_cell = c(75, 75, 250, 250), 
                                    show_plot = F,
                                    default_save_name = "thresholds")
s10.thresholds[["ggplot"]][["theme"]][["legend.position"]] <- "none" # Remove legend

s10.metric.plots <- ggarrange(s10.detected, s10.det.along, s10.libvdet, s10.thresholds$ggplot) # Combine all plots
rm(s10.detected, s10.det.along, s10.libvdet, s10.thresholds); s10.metric.plots
ggsave("./project/outcomes/ms04_metrics.png", plot = s10.metric.plots, scale = 2, width = 1920, height = 1080, units = "px")

## Filter sample
s10.filtered <- filterGiotto(s10, expression_values = "raw", expression_threshold = 1,
                            feat_det_in_min_cells = 50,
                            min_det_feats_per_cell = 75)
save(s10.filtered, file = "./project/material/filtered_samples/s10_filtered.R") # load("./project/filtered_samples/s10_filtered.R")

## Deleted spots
length(s10@cell_ID$cell) - length(s10.filtered@cell_ID$cell) # 187 deleted spots

### Visualization plot
s10.spots.plot <- spatPlot2D(s10, cell_color = ("lightgrey"), point_size = 2,
                            select_cells = s10.filtered@cell_ID$cell, # Kept spots
                            other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                            title = "Deleted spots (sample 10)"); s10.spots.plot
ggsave("./project/outcomes/ms04_deleted.png", plot = s10.spots.plot, scale = 2.5, width = 1920, height = 1080, units = "px")


####--SAMPLE 11--####

# Obtain data
s11 <- createGiottoVisiumObject(
  visium_dir = sample.path[11], #ms05
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = createGiottoInstructions(save_plot = F, show_plot = F))

## Filter over tissue genes
in.tissue <- pDataDT(s11)[in_tissue == 1]$cell_ID 
if (sum((in.tissue == pDataDT(s11)$cell_ID) == F) == 0) {
  cat("All genes are over tissue")
  rm(in.tissue)
} else {
  s11 <- subsetGiotto(s11, cell_ids = in.tissue)
  rm(in.tissue)
} # All genes are over tissue

s11.matrix <- s11@expression[["cell"]][["rna"]][["raw"]]@exprMat # Extract expression matrix from Giotto object


# Quality control 
## Create cell metadata
reads.depth <- apply(s11.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s11.matrix != 0, 2, sum) # Detected genes per cell

is_mito <- grepl("(^MT-)|(^Mt-)|(^mt-)", s11@feat_ID$rna) # Identify mitochondrial genes
if (sum(is_mito == T) != 0) {
  s11@feat_ID$is_mito <- is_mito
  s11.mito <- subsetGiotto(s11, feat_ids = s11@feat_ID$rna[is_mito == T])
  cat(sum(is_mito == T), "mitochondrial genes found:", s11.mito@feat_ID$rna)
  rm(is_mito)
} else {
  cat('No mitochondrial genes found')
  rm(is_mito)
} # 13 mitochondrial genes found

expr.sum <- apply(s11.matrix, 2, sum) # Sum of gene expression values per cell
mito.sum <- apply(s11.matrix[s11.mito@feat_ID$rna, ], 2, sum) # Sum of mito gene expression values per cell
mito.prop <- round(mito.sum/expr.sum, 3) # Mito genes proportion

s11.cell.meta <- data.frame(nReads = reads.depth, nGenes = genes.per.cell, mitoProp = mito.prop) # Write data
rownames(s11.cell.meta) <- s11@cell_ID$cell # Add cell names
rm(reads.depth, genes.per.cell, expr.sum, mito.sum, mito.prop, s11.mito) # Rm variables


## Calcs
### Mito proportion threshold
which(s11.cell.meta$mitoProp > 0.28) # 28 cells exceed 0.28 mito genes proportion

## Libsize vs Detected genes correlation
with(s11.cell.meta, cor.test(nReads, nGenes)) # Check correlation (0.9861)


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
s11.libvdet <- ggplot(s11.cell.meta, aes(nReads, nGenes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = "Corr = 0.9861", x = 19000, y = -200)

### Threshold evaluation
s11.thresholds <- filterCombinations(s11, expression_thresholds = 1,
                                    feat_det_in_min_cells = c(50, 100, 50, 100),
                                    min_det_feats_per_cell = c(125, 125, 250, 250), 
                                    show_plot = F,
                                    default_save_name = "thresholds")
s11.thresholds[["ggplot"]][["theme"]][["legend.position"]] <- "none" # Remove legend

s11.metric.plots <- ggarrange(s11.detected, s11.det.along, s11.libvdet, s11.thresholds$ggplot) # Combine all plots
rm(s11.detected, s11.det.along, s11.libvdet, s11.thresholds); s11.metric.plots
ggsave("./project/outcomes/ms05_metrics.png", plot = s11.metric.plots, scale = 2, width = 1920, height = 1080, units = "px")

## Filter sample
s11.filtered <- filterGiotto(s11, expression_values = "raw", expression_threshold = 1,
                            feat_det_in_min_cells = 50,
                            min_det_feats_per_cell = 150)
save(s11.filtered, file = "./project/material/filtered_samples/s11_filtered.R") # load("./project/filtered_samples/s11_filtered.R")

## Deleted spots
length(s11@cell_ID$cell) - length(s11.filtered@cell_ID$cell) # 59 deleted spots

### Visualization plot
s11.spots.plot <- spatPlot2D(s11, cell_color = ("lightgrey"), point_size = 2,
                            select_cells = s11.filtered@cell_ID$cell, # Kept spots
                            other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                            title = "Deleted spots (sample 11)"); s11.spots.plot
ggsave("./project/outcomes/ms05_deleted.png", plot = s11.spots.plot, scale = 2.5, width = 1920, height = 1080, units = "px")


####--SAMPLE 12--####

# Obtain data
s12 <- createGiottoVisiumObject(
  visium_dir = sample.path[12], #ms06
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = createGiottoInstructions(save_plot = F, show_plot = F))

## Filter over tissue genes
in.tissue <- pDataDT(s12)[in_tissue == 1]$cell_ID 
if (sum((in.tissue == pDataDT(s12)$cell_ID) == F) == 0) {
  cat("All genes are over tissue")
  rm(in.tissue)
} else {
  s12 <- subsetGiotto(s12, cell_ids = in.tissue)
  rm(in.tissue)
} # All genes are over tissue

s12.matrix <- s12@expression[["cell"]][["rna"]][["raw"]]@exprMat # Extract expression matrix from Giotto object


# Quality control 
## Create cell metadata
reads.depth <- apply(s12.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s12.matrix != 0, 2, sum) # Detected genes per cell

is_mito <- grepl("(^MT-)|(^Mt-)|(^mt-)", s12@feat_ID$rna) # Identify mitochondrial genes
if (sum(is_mito == T) != 0) {
  s12@feat_ID$is_mito <- is_mito
  s12.mito <- subsetGiotto(s12, feat_ids = s12@feat_ID$rna[is_mito == T])
  cat(sum(is_mito == T), "mitochondrial genes found:", s12.mito@feat_ID$rna)
  rm(is_mito)
} else {
  cat('No mitochondrial genes found')
  rm(is_mito)
} # 13 mitochondrial genes found

expr.sum <- apply(s12.matrix, 2, sum) # Sum of gene expression values per cell
mito.sum <- apply(s12.matrix[s12.mito@feat_ID$rna, ], 2, sum) # Sum of mito gene expression values per cell
mito.prop <- round(mito.sum/expr.sum, 3) # Mito genes proportion

s12.cell.meta <- data.frame(nReads = reads.depth, nGenes = genes.per.cell, mitoProp = mito.prop) # Write data
rownames(s12.cell.meta) <- s12@cell_ID$cell # Add cell names
rm(reads.depth, genes.per.cell, expr.sum, mito.sum, mito.prop, s12.mito) # Rm variables


## Calcs
### Mito proportion threshold
which(s12.cell.meta$mitoProp > 0.28) # 335 cells exceed 0.28 mito genes proportion

## Libsize vs Detected genes correlation
with(s12.cell.meta, cor.test(nReads, nGenes)) # Check correlation (0.9683)


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
s12.libvdet <- ggplot(s12.cell.meta, aes(nReads, nGenes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = "Corr = 0.9683", x = 18000, y = -200)

### Threshold evaluation
s12.thresholds <- filterCombinations(s12, expression_thresholds = 1,
                                    feat_det_in_min_cells = c(50, 100, 50, 100),
                                    min_det_feats_per_cell = c(250, 250, 500, 500), 
                                    show_plot = F,
                                    default_save_name = "thresholds")
s12.thresholds[["ggplot"]][["theme"]][["legend.position"]] <- "none" # Remove legend

s12.metric.plots <- ggarrange(s12.detected, s12.det.along, s12.libvdet, s12.thresholds$ggplot) # Combine all plots
rm(s12.detected, s12.det.along, s12.libvdet, s12.thresholds); s12.metric.plots
ggsave("./project/outcomes/ms06_metrics.png", plot = s12.metric.plots, scale = 2, width = 1920, height = 1080, units = "px")

## Filter sample
s12.filtered <- filterGiotto(s12, expression_values = "raw", expression_threshold = 1,
                            feat_det_in_min_cells = 50,
                            min_det_feats_per_cell = 250)
save(s12.filtered, file = "./project/material/filtered_samples/s12_filtered.R") # load("./project/filtered_samples/s12_filtered.R")

## Deleted spots
length(s12@cell_ID$cell) - length(s12.filtered@cell_ID$cell) # 131 deleted spots

### Visualization plot
s12.spots.plot <- spatPlot2D(s12, cell_color = ("lightgrey"), point_size = 2,
                            select_cells = s12.filtered@cell_ID$cell, # Kept spots
                            other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                            title = "Deleted spots (sample 12)"); s12.spots.plot
ggsave("./project/outcomes/ms06_deleted.png", plot = s12.spots.plot, scale = 2.5, width = 1920, height = 1080, units = "px")


####--SAMPLE 13 (REVISAR)--####

# Obtain data
s13 <- createGiottoVisiumObject(
  visium_dir = sample.path[13], #ms07
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = createGiottoInstructions(save_plot = F, show_plot = F))

## Filter over tissue genes
in.tissue <- pDataDT(s13)[in_tissue == 1]$cell_ID 
if (sum((in.tissue == pDataDT(s13)$cell_ID) == F) == 0) {
  cat("All genes are over tissue")
  rm(in.tissue)
} else {
  s13 <- subsetGiotto(s13, cell_ids = in.tissue)
  rm(in.tissue)
} # All genes are over tissue

s13.matrix <- s13@expression[["cell"]][["rna"]][["raw"]]@exprMat # Extract expression matrix from Giotto object


# Quality control 
## Create cell metadata
reads.depth <- apply(s13.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s13.matrix != 0, 2, sum) # Detected genes per cell

is_mito <- grepl("(^MT-)|(^Mt-)|(^mt-)", s13@feat_ID$rna) # Identify mitochondrial genes
if (sum(is_mito == T) != 0) {
  s13@feat_ID$is_mito <- is_mito
  s13.mito <- subsetGiotto(s13, feat_ids = s13@feat_ID$rna[is_mito == T])
  cat(sum(is_mito == T), "mitochondrial genes found:", s13.mito@feat_ID$rna)
  rm(is_mito)
} else {
  cat('No mitochondrial genes found')
  rm(is_mito)
} # 13 mitochondrial genes found

expr.sum <- apply(s13.matrix, 2, sum) # Sum of gene expression values per cell
mito.sum <- apply(s13.matrix[s13.mito@feat_ID$rna, ], 2, sum) # Sum of mito gene expression values per cell
mito.prop <- round(mito.sum/expr.sum, 3) # Mito genes proportion

s13.cell.meta <- data.frame(nReads = reads.depth, nGenes = genes.per.cell, mitoProp = mito.prop) # Write data
rownames(s13.cell.meta) <- s13@cell_ID$cell # Add cell names
rm(reads.depth, genes.per.cell, expr.sum, mito.sum, mito.prop, s13.mito) # Rm variables


## Calcs
### Mito proportion threshold
which(s13.cell.meta$mitoProp > 0.28) # 2010 cells exceed 0.28 mito genes proportion

## Libsize vs Detected genes correlation
with(s13.cell.meta, cor.test(nReads, nGenes)) # Check correlation (0.9849)


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
s13.libvdet <- ggplot(s13.cell.meta, aes(nReads, nGenes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = "Corr = 0.9849", x = 18000, y = -200)

### Threshold evaluation
s13.thresholds <- filterCombinations(s13, expression_thresholds = 1,
                                    feat_det_in_min_cells = c(50, 100, 50, 100),
                                    min_det_feats_per_cell = c(125, 125, 250, 250), 
                                    show_plot = F,
                                    default_save_name = "thresholds")
s13.thresholds[["ggplot"]][["theme"]][["legend.position"]] <- "none" # Remove legend

s13.metric.plots <- ggarrange(s13.detected, s13.det.along, s13.libvdet, s13.thresholds$ggplot) # Combine all plots
rm(s13.detected, s13.det.along, s13.libvdet, s13.thresholds); s13.metric.plots
ggsave("./project/outcomes/ms07_metrics.png", plot = s13.metric.plots, scale = 2, width = 1920, height = 1080, units = "px")

## Filter sample
s13.filtered <- filterGiotto(s13, expression_values = "raw", expression_threshold = 1,
                            feat_det_in_min_cells = 50,
                            min_det_feats_per_cell = 120)
save(s13.filtered, file = "./project/material/filtered_samples/s13_filtered.R") # load("./project/filtered_samples/s13_filtered.R")

## Deleted spots
length(s13@cell_ID$cell) - length(s13.filtered@cell_ID$cell) # 98 deleted spots

### Visualization plot
s13.spots.plot <- spatPlot2D(s13, cell_color = ("lightgrey"), point_size = 2,
                            select_cells = s13.filtered@cell_ID$cell, # Kept spots
                            other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                            title = "Deleted spots (sample 13)"); s13.spots.plot
ggsave("./project/outcomes/ms07_deleted.png", plot = s13.spots.plot, scale = 2.5, width = 1920, height = 1080, units = "px")


####--SAMPLE 14--####

# Obtain data
s14 <- createGiottoVisiumObject(
  visium_dir = sample.path[14], #ms08
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = createGiottoInstructions(save_plot = F, show_plot = F))

## Filter over tissue genes
in.tissue <- pDataDT(s14)[in_tissue == 1]$cell_ID 
if (sum((in.tissue == pDataDT(s14)$cell_ID) == F) == 0) {
  cat("All genes are over tissue")
  rm(in.tissue)
} else {
  s14 <- subsetGiotto(s14, cell_ids = in.tissue)
  rm(in.tissue)
} # All genes are over tissue

s14.matrix <- s14@expression[["cell"]][["rna"]][["raw"]]@exprMat # Extract expression matrix from Giotto object


# Quality control 
## Create cell metadata
reads.depth <- apply(s14.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s14.matrix != 0, 2, sum) # Detected genes per cell

is_mito <- grepl("(^MT-)|(^Mt-)|(^mt-)", s14@feat_ID$rna) # Identify mitochondrial genes
if (sum(is_mito == T) != 0) {
  s14@feat_ID$is_mito <- is_mito
  s14.mito <- subsetGiotto(s14, feat_ids = s14@feat_ID$rna[is_mito == T])
  cat(sum(is_mito == T), "mitochondrial genes found:", s14.mito@feat_ID$rna)
  rm(is_mito)
} else {
  cat('No mitochondrial genes found')
  rm(is_mito)
} # 13 mitochondrial genes found

expr.sum <- apply(s14.matrix, 2, sum) # Sum of gene expression values per cell
mito.sum <- apply(s14.matrix[s14.mito@feat_ID$rna, ], 2, sum) # Sum of mito gene expression values per cell
mito.prop <- round(mito.sum/expr.sum, 3) # Mito genes proportion

s14.cell.meta <- data.frame(nReads = reads.depth, nGenes = genes.per.cell, mitoProp = mito.prop) # Write data
rownames(s14.cell.meta) <- s14@cell_ID$cell # Add cell names
rm(reads.depth, genes.per.cell, expr.sum, mito.sum, mito.prop, s14.mito) # Rm variables


## Calcs
### Mito proportion threshold
which(s14.cell.meta$mitoProp > 0.28) # 48 cells exceed 0.28 mito genes proportion

## Libsize vs Detected genes correlation
with(s14.cell.meta, cor.test(nReads, nGenes)) # Check correlation (0.9844)


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
s14.libvdet <- ggplot(s14.cell.meta, aes(nReads, nGenes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = "Corr = 0.9844", x = 13000, y = -200)

### Threshold evaluation
s14.thresholds <- filterCombinations(s14, expression_thresholds = 1,
                                    feat_det_in_min_cells = c(50, 100, 50, 100),
                                    min_det_feats_per_cell = c(250, 250, 500, 500), 
                                    show_plot = F,
                                    default_save_name = "thresholds")
s14.thresholds[["ggplot"]][["theme"]][["legend.position"]] <- "none" # Remove legend

s14.metric.plots <- ggarrange(s14.detected, s14.det.along, s14.libvdet, s14.thresholds$ggplot) # Combine all plots
rm(s14.detected, s14.det.along, s14.libvdet, s14.thresholds); s14.metric.plots
ggsave("./project/outcomes/ms08_metrics.png", plot = s14.metric.plots, scale = 2, width = 1920, height = 1080, units = "px")

## Filter sample
s14.filtered <- filterGiotto(s14, expression_values = "raw", expression_threshold = 1,
                            feat_det_in_min_cells = 50,
                            min_det_feats_per_cell = 500)
save(s14.filtered, file = "./project/material/filtered_samples/s14_filtered.R") # load("./project/filtered_samples/s14_filtered.R")

## Deleted spots
length(s14@cell_ID$cell) - length(s14.filtered@cell_ID$cell) # 284 deleted spots

### Visualization plot
s14.spots.plot <- spatPlot2D(s14, cell_color = ("lightgrey"), point_size = 2,
                            select_cells = s14.filtered@cell_ID$cell, # Kept spots
                            other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                            title = "Deleted spots (sample 14)"); s14.spots.plot
ggsave("./project/outcomes/ms08_deleted.png", plot = s14.spots.plot, scale = 2.5, width = 1920, height = 1080, units = "px")


####--SAMPLE 15--####

# Obtain data
s15 <- createGiottoVisiumObject(
  visium_dir = sample.path[15], #ms09
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = createGiottoInstructions(save_plot = F, show_plot = F))

## Filter over tissue genes
in.tissue <- pDataDT(s15)[in_tissue == 1]$cell_ID 
if (sum((in.tissue == pDataDT(s15)$cell_ID) == F) == 0) {
  cat("All genes are over tissue")
  rm(in.tissue)
} else {
  s15 <- subsetGiotto(s15, cell_ids = in.tissue)
  rm(in.tissue)
} # All genes are over tissue

s15.matrix <- s15@expression[["cell"]][["rna"]][["raw"]]@exprMat # Extract expression matrix from Giotto object


# Quality control 
## Create cell metadata
reads.depth <- apply(s15.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s15.matrix != 0, 2, sum) # Detected genes per cell

is_mito <- grepl("(^MT-)|(^Mt-)|(^mt-)", s15@feat_ID$rna) # Identify mitochondrial genes
if (sum(is_mito == T) != 0) {
  s15@feat_ID$is_mito <- is_mito
  s15.mito <- subsetGiotto(s15, feat_ids = s15@feat_ID$rna[is_mito == T])
  cat(sum(is_mito == T), "mitochondrial genes found:", s15.mito@feat_ID$rna)
  rm(is_mito)
} else {
  cat('No mitochondrial genes found')
  rm(is_mito)
} # 13 mitochondrial genes found

expr.sum <- apply(s15.matrix, 2, sum) # Sum of gene expression values per cell
mito.sum <- apply(s15.matrix[s15.mito@feat_ID$rna, ], 2, sum) # Sum of mito gene expression values per cell
mito.prop <- round(mito.sum/expr.sum, 3) # Mito genes proportion

s15.cell.meta <- data.frame(nReads = reads.depth, nGenes = genes.per.cell, mitoProp = mito.prop) # Write data
rownames(s15.cell.meta) <- s15@cell_ID$cell # Add cell names
rm(reads.depth, genes.per.cell, expr.sum, mito.sum, mito.prop, s15.mito) # Rm variables


## Calcs
### Mito proportion threshold
which(s15.cell.meta$mitoProp > 0.28) # No cells exceed 0.28 mito genes proportion

## Libsize vs Detected genes correlation
with(s15.cell.meta, cor.test(nReads, nGenes)) # Check correlation (0.9889)


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
s15.libvdet <- ggplot(s15.cell.meta, aes(nReads, nGenes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = "Corr = 0.9889", x = 6000, y = -200)

### Threshold evaluation
s15.thresholds <- filterCombinations(s15, expression_thresholds = 1,
                                    feat_det_in_min_cells = c(50, 100, 50, 100),
                                    min_det_feats_per_cell = c(250, 250, 500, 500), 
                                    show_plot = F,
                                    default_save_name = "thresholds")
s15.thresholds[["ggplot"]][["theme"]][["legend.position"]] <- "none" # Remove legend

s15.metric.plots <- ggarrange(s15.detected, s15.det.along, s15.libvdet, s15.thresholds$ggplot) # Combine all plots
rm(s15.detected, s15.det.along, s15.libvdet, s15.thresholds); s15.metric.plots
ggsave("./project/outcomes/ms09_metrics.png", plot = s15.metric.plots, scale = 2, width = 1920, height = 1080, units = "px")

## Filter sample
s15.filtered <- filterGiotto(s15, expression_values = "raw", expression_threshold = 1,
                            feat_det_in_min_cells = 50,
                            min_det_feats_per_cell = 250)
save(s15.filtered, file = "./project/material/filtered_samples/s15_filtered.R") # load("./project/filtered_samples/s15_filtered.R")

## Deleted spots
length(s15@cell_ID$cell) - length(s15.filtered@cell_ID$cell) # 42 deleted spots

### Visualization plot
s15.spots.plot <- spatPlot2D(s15, cell_color = ("lightgrey"), point_size = 2,
                            select_cells = s15.filtered@cell_ID$cell, # Kept spots
                            other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                            title = "Deleted spots (sample 15)"); s15.spots.plot
ggsave("./project/outcomes/ms09_deleted.png", plot = s15.spots.plot, scale = 2.5, width = 1920, height = 1080, units = "px")


####--SAMPLE 16--####

# Obtain data
s16 <- createGiottoVisiumObject(
  visium_dir = sample.path[16], #ms10
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = createGiottoInstructions(save_plot = F, show_plot = F))

## Filter over tissue genes
in.tissue <- pDataDT(s16)[in_tissue == 1]$cell_ID 
if (sum((in.tissue == pDataDT(s16)$cell_ID) == F) == 0) {
  cat("All genes are over tissue")
  rm(in.tissue)
} else {
  s16 <- subsetGiotto(s16, cell_ids = in.tissue)
  rm(in.tissue)
} # All genes are over tissue

s16.matrix <- s16@expression[["cell"]][["rna"]][["raw"]]@exprMat # Extract expression matrix from Giotto object


# Quality control 
## Create cell metadata
reads.depth <- apply(s16.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s16.matrix != 0, 2, sum) # Detected genes per cell

is_mito <- grepl("(^MT-)|(^Mt-)|(^mt-)", s16@feat_ID$rna) # Identify mitochondrial genes
if (sum(is_mito == T) != 0) {
  s16@feat_ID$is_mito <- is_mito
  s16.mito <- subsetGiotto(s16, feat_ids = s16@feat_ID$rna[is_mito == T])
  cat(sum(is_mito == T), "mitochondrial genes found:", s16.mito@feat_ID$rna)
  rm(is_mito)
} else {
  cat('No mitochondrial genes found')
  rm(is_mito)
} # 13 mitochondrial genes found

expr.sum <- apply(s16.matrix, 2, sum) # Sum of gene expression values per cell
mito.sum <- apply(s16.matrix[s16.mito@feat_ID$rna, ], 2, sum) # Sum of mito gene expression values per cell
mito.prop <- round(mito.sum/expr.sum, 3) # Mito genes proportion

s16.cell.meta <- data.frame(nReads = reads.depth, nGenes = genes.per.cell, mitoProp = mito.prop) # Write data
rownames(s16.cell.meta) <- s16@cell_ID$cell # Add cell names
rm(reads.depth, genes.per.cell, expr.sum, mito.sum, mito.prop, s16.mito) # Rm variables


## Calcs
### Mito proportion threshold
which(s16.cell.meta$mitoProp > 0.28) # 117 cells exceed 0.28 mito genes proportion

## Libsize vs Detected genes correlation
with(s16.cell.meta, cor.test(nReads, nGenes)) # Check correlation (0.9855)


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
s16.libvdet <- ggplot(s16.cell.meta, aes(nReads, nGenes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = "Corr = 0.9855", x = 7500, y = -200)

### Threshold evaluation
s16.thresholds <- filterCombinations(s16, expression_thresholds = 1,
                                    feat_det_in_min_cells = c(50, 100, 50, 100),
                                    min_det_feats_per_cell = c(250, 250, 500, 500), 
                                    show_plot = F,
                                    default_save_name = "thresholds")
s16.thresholds[["ggplot"]][["theme"]][["legend.position"]] <- "none" # Remove legend

s16.metric.plots <- ggarrange(s16.detected, s16.det.along, s16.libvdet, s16.thresholds$ggplot) # Combine all plots
rm(s16.detected, s16.det.along, s16.libvdet, s16.thresholds); s16.metric.plots
ggsave("./project/outcomes/ms10_metrics.png", plot = s16.metric.plots, scale = 2, width = 1920, height = 1080, units = "px")

## Filter sample
s16.filtered <- filterGiotto(s16, expression_values = "raw", expression_threshold = 1,
                            feat_det_in_min_cells = 50,
                            min_det_feats_per_cell = 300)
save(s16.filtered, file = "./project/material/filtered_samples/s16_filtered.R") # load("./project/filtered_samples/s16_filtered.R")

## Deleted spots
length(s16@cell_ID$cell) - length(s16.filtered@cell_ID$cell) # 121 deleted spots

### Visualization plot
s16.spots.plot <- spatPlot2D(s16, cell_color = ("lightgrey"), point_size = 2,
                            select_cells = s16.filtered@cell_ID$cell, # Kept spots
                            other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                            title = "Deleted spots (sample 16)"); s16.spots.plot
ggsave("./project/outcomes/ms10_deleted.png", plot = s16.spots.plot, scale = 2.5, width = 1920, height = 1080, units = "px")


####--SAMPLE 17 (REVISAR)--####

# Obtain data
s17 <- createGiottoVisiumObject(
  visium_dir = sample.path[17], #ms11
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = createGiottoInstructions(save_plot = F, show_plot = F))

## Filter over tissue genes
in.tissue <- pDataDT(s17)[in_tissue == 1]$cell_ID 
if (sum((in.tissue == pDataDT(s17)$cell_ID) == F) == 0) {
  cat("All genes are over tissue")
  rm(in.tissue)
} else {
  s17 <- subsetGiotto(s17, cell_ids = in.tissue)
  rm(in.tissue)
} # All genes are over tissue

s17.matrix <- s17@expression[["cell"]][["rna"]][["raw"]]@exprMat # Extract expression matrix from Giotto object


# Quality control 
## Create cell metadata
reads.depth <- apply(s17.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s17.matrix != 0, 2, sum) # Detected genes per cell

is_mito <- grepl("(^MT-)|(^Mt-)|(^mt-)", s17@feat_ID$rna) # Identify mitochondrial genes
if (sum(is_mito == T) != 0) {
  s17@feat_ID$is_mito <- is_mito
  s17.mito <- subsetGiotto(s17, feat_ids = s17@feat_ID$rna[is_mito == T])
  cat(sum(is_mito == T), "mitochondrial genes found:", s17.mito@feat_ID$rna)
  rm(is_mito)
} else {
  cat('No mitochondrial genes found')
  rm(is_mito)
} # 13 mitochondrial genes found

expr.sum <- apply(s17.matrix, 2, sum) # Sum of gene expression values per cell
mito.sum <- apply(s17.matrix[s17.mito@feat_ID$rna, ], 2, sum) # Sum of mito gene expression values per cell
mito.prop <- round(mito.sum/expr.sum, 3) # Mito genes proportion

s17.cell.meta <- data.frame(nReads = reads.depth, nGenes = genes.per.cell, mitoProp = mito.prop) # Write data
rownames(s17.cell.meta) <- s17@cell_ID$cell # Add cell names
rm(reads.depth, genes.per.cell, expr.sum, mito.sum, mito.prop, s17.mito) # Rm variables


## Calcs
### Mito proportion threshold
which(s17.cell.meta$mitoProp > 0.28) # 988 cells exceed 0.28 mito genes proportion

## Libsize vs Detected genes correlation
with(s17.cell.meta, cor.test(nReads, nGenes)) # Check correlation (0.9756)


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
s17.libvdet <- ggplot(s17.cell.meta, aes(nReads, nGenes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = "Corr = 0.9756", x = 28000, y = -200)

### Threshold evaluation
s17.thresholds <- filterCombinations(s17, expression_thresholds = 1,
                                    feat_det_in_min_cells = c(50, 100, 50, 100),
                                    min_det_feats_per_cell = c(500, 500, 1000, 1000), 
                                    show_plot = F,
                                    default_save_name = "thresholds")
s17.thresholds[["ggplot"]][["theme"]][["legend.position"]] <- "none" # Remove legend

s17.metric.plots <- ggarrange(s17.detected, s17.det.along, s17.libvdet, s17.thresholds$ggplot) # Combine all plots
rm(s17.detected, s17.det.along, s17.libvdet, s17.thresholds); s17.metric.plots
ggsave("./project/outcomes/ms11_metrics.png", plot = s17.metric.plots, scale = 2, width = 1920, height = 1080, units = "px")

## Filter sample
s17.filtered <- filterGiotto(s17, expression_values = "raw", expression_threshold = 1,
                            feat_det_in_min_cells = 50,
                            min_det_feats_per_cell = 1000)
save(s17.filtered, file = "./project/material/filtered_samples/s17_filtered.R") # load("./project/filtered_samples/s17_filtered.R")

## Deleted spots
length(s17@cell_ID$cell) - length(s17.filtered@cell_ID$cell) # 240 deleted spots

### Visualization plot
s17.spots.plot <- spatPlot2D(s17, cell_color = ("lightgrey"), point_size = 2,
                            select_cells = s17.filtered@cell_ID$cell, # Kept spots
                            other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                            title = "Deleted spots (sample 17)"); s17.spots.plot
ggsave("./project/outcomes/ms11_deleted.png", plot = s17.spots.plot, scale = 2.5, width = 1920, height = 1080, units = "px")


####--SAMPLE 18 (REVISAR)--####

# Obtain data
s18 <- createGiottoVisiumObject(
  visium_dir = sample.path[18], #ms12
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = createGiottoInstructions(save_plot = F, show_plot = F))

## Filter over tissue genes
in.tissue <- pDataDT(s18)[in_tissue == 1]$cell_ID 
if (sum((in.tissue == pDataDT(s18)$cell_ID) == F) == 0) {
  cat("All genes are over tissue")
  rm(in.tissue)
} else {
  s18 <- subsetGiotto(s18, cell_ids = in.tissue)
  rm(in.tissue)
} # All genes are over tissue

s18.matrix <- s18@expression[["cell"]][["rna"]][["raw"]]@exprMat # Extract expression matrix from Giotto object


# Quality control 
## Create cell metadata
reads.depth <- apply(s18.matrix, 2, sum) # Reads per cell
genes.per.cell <- apply(s18.matrix != 0, 2, sum) # Detected genes per cell

is_mito <- grepl("(^MT-)|(^Mt-)|(^mt-)", s18@feat_ID$rna) # Identify mitochondrial genes
if (sum(is_mito == T) != 0) {
  s18@feat_ID$is_mito <- is_mito
  s18.mito <- subsetGiotto(s18, feat_ids = s18@feat_ID$rna[is_mito == T])
  cat(sum(is_mito == T), "mitochondrial genes found:", s18.mito@feat_ID$rna)
  rm(is_mito)
} else {
  cat('No mitochondrial genes found')
  rm(is_mito)
} # 13 mitochondrial genes found

expr.sum <- apply(s18.matrix, 2, sum) # Sum of gene expression values per cell
mito.sum <- apply(s18.matrix[s18.mito@feat_ID$rna, ], 2, sum) # Sum of mito gene expression values per cell
mito.prop <- round(mito.sum/expr.sum, 3) # Mito genes proportion

s18.cell.meta <- data.frame(nReads = reads.depth, nGenes = genes.per.cell, mitoProp = mito.prop) # Write data
rownames(s18.cell.meta) <- s18@cell_ID$cell # Add cell names
rm(reads.depth, genes.per.cell, expr.sum, mito.sum, mito.prop, s18.mito) # Rm variables


## Calcs
### Mito proportion threshold
which(s18.cell.meta$mitoProp > 0.28) # 353 cells exceed 0.28 mito genes proportion

## Libsize vs Detected genes correlation
with(s18.cell.meta, cor.test(nReads, nGenes)) # Check correlation (0.9901)


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
s18.libvdet <- ggplot(s18.cell.meta, aes(nReads, nGenes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = "Corr = 0.9901", x = 12000, y = -200)

### Threshold evaluation
s18.thresholds <- filterCombinations(s18, expression_thresholds = 1,
                                    feat_det_in_min_cells = c(50, 100, 50, 100),
                                    min_det_feats_per_cell = c(250, 250, 500, 500), 
                                    show_plot = F,
                                    default_save_name = "thresholds")
s18.thresholds[["ggplot"]][["theme"]][["legend.position"]] <- "none" # Remove legend

s18.metric.plots <- ggarrange(s18.detected, s18.det.along, s18.libvdet, s18.thresholds$ggplot) # Combine all plots
rm(s18.detected, s18.det.along, s18.libvdet, s18.thresholds); s18.metric.plots
ggsave("./project/outcomes/ms12_metrics.png", plot = s18.metric.plots, scale = 2, width = 1920, height = 1080, units = "px")

## Filter sample
s18.filtered <- filterGiotto(s18, expression_values = "raw", expression_threshold = 1,
                            feat_det_in_min_cells = 50,
                            min_det_feats_per_cell = 300)
save(s18.filtered, file = "./project/material/filtered_samples/s18_filtered.R") # load("./project/filtered_samples/s18_filtered.R")

## Deleted spots
length(s18@cell_ID$cell) - length(s18.filtered@cell_ID$cell) # 257 deleted spots

### Visualization plot
s18.spots.plot <- spatPlot2D(s18, cell_color = ("lightgrey"), point_size = 2,
                            select_cells = s18.filtered@cell_ID$cell, # Kept spots
                            other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                            title = "Deleted spots (sample 18)"); s18.spots.plot
ggsave("./project/outcomes/ms12_deleted.png", plot = s18.spots.plot, scale = 2.5, width = 1920, height = 1080, units = "px")

