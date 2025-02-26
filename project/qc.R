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

dir <- "./project/material/GSE165098"
sample.ids <- c("LPS1", "LPS2", "LPS3", "saline1", "saline2", "saline3")
sample.path <- file.path(dir, sample.ids); rm(dir, sample.ids)


####--SAMPLE 1--####

# Obtain data
s1 <- createGiottoVisiumObject(
  visium_dir = sample.path[1], #LPS1
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = createGiottoInstructions(save_plot = F, show_plot = F))

## Filter over tissue genes
in.tissue <- pDataDT(s1)[in_tissue == 1]$cell_ID 
if (sum((in.tissue == pDataDT(s1)$cell_ID) == F) == 0) {
  cat("All genes are over tissue")
  rm(in.tissue)
} else {
  s1 <- subsetGiotto(s1, cell_ids = in.tissue)
  rm(in.tissue)
} # All genes are over tissue

s1.matrix <- s1@expression[["cell"]][["rna"]][["raw"]]@exprMat # Extract expression matrix from Giotto object


# Quality control 
## Create cell metadata
reads.depth <- apply(s1.matrix, 2, sum) # Reads per cell
genes.per.cell <- 21215 - apply(s1.matrix == 0, 2, sum) # Detected genes per cell
length(genes.per.cell); length(reads.depth) # Check objects length

is_mito <- grepl("(^MT-)|(^Mt-)|(^mt-)", s1@feat_ID$rna) # Identify mitochondrial genes
if (sum(is_mito == T) != 0) {
  s1@feat_ID$is_mito <- is_mito
  s1.mito <- subsetGiotto(s1, feat_ids = s1@feat_ID$rna[is_mito == T])
  cat(sum(is_mito == T), "mitochondrial genes found:", s1.mito@feat_ID$rna)
  rm(is_mito)
} else {
  cat('No mitochondrial genes found')
  rm(is_mito)
} # 13 mitochondrial genes found

expr.sum <- apply(s1.matrix, 2, sum) # Sum of gene expression values per cell
mito.sum <- apply(s1.matrix[s1.mito@feat_ID$rna, ], 2, sum) # Sum of mito gene expression values per cell
mito.prop <- round(mito.sum/expr.sum, 3) # Mito genes proportion

s1.cell.meta <- data.frame(nReads = reads.depth, nGenes = genes.per.cell, mitoProp = mito.prop) # Write data
rownames(s1.cell.meta) <- s1@cell_ID$cell # Add cell names
rm(reads.depth, genes.per.cell, expr.sum, mito.sum, mito.prop, s1.mito) # Rm variables


## Calcs
### Mito proportion threshold
which(s1.cell.meta$mitoProp > 0.28) # No cells exceed 0.28 mito genes proportion

## Libsize vs Detected genes correlation
with(s1.cell.meta, cor.test(nReads, nGenes)) # Check correlation (0.9355)


## Plots
### Detected genes
s1.detected <- filterDistributions(s1, detection = "cells", nr_bins = 50,
                                   method = "threshold", expression_threshold = 1, 
                                   default_save_name = "detected_genes")

### Genes detection along cells
s1.det.along <- filterDistributions(s1, detection = "feats", nr_bins = 50,
                                  method = "threshold", 
                                  default_save_name = "library_size")

### Libsize vs Detected genes
s1.libvdet <- ggplot(s1.cell.meta, aes(nReads, nGenes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = "Corr = 0.9355", x = 60000, y = 500)

### Threshold evaluation
s1.thresholds <- filterCombinations(s1, expression_thresholds = 1,
                                    feat_det_in_min_cells = c(50, 100, 50, 100),
                                    min_det_feats_per_cell = c(2500, 2500, 5000, 5000), 
                                    show_plot = F,
                                    default_save_name = "thresholds")
s1.thresholds[["ggplot"]][["theme"]][["legend.position"]] <- "none" # Remove legend

s1.metric.plots <- ggarrange(s1.detected, s1.det.along, s1.libvdet, s1.thresholds$ggplot) # Combine all plots
rm(s1.detected, s1.det.along, s1.libvdet, s1.thresholds); s1.metric.plots
ggsave("./project/outcomes/LPS1_metrics.png", plot = s1.metric.plots, scale = 2, width = 1920, height = 1080, units = "px")

## Filter sample
s1.filtered <- filterGiotto(s1, expression_values = "raw", expression_threshold = 1,
                            feat_det_in_min_cells = 50,
                            min_det_feats_per_cell = 2500)
save(s1.filtered, file = "./project/material/filtered_samples/s1_filtered.R") # load("./project/filtered_samples/s1_filtered.R")

## Deleted spots
length(s1@cell_ID$cell) - length(s1.filtered@cell_ID$cell) # 243 deleted spots

### Visualization plot
s1.spots.plot <- spatPlot2D(s1, cell_color = ("lightgrey"), point_size = 2,
                         select_cells = s1.filtered@cell_ID$cell, # Kept spots
                         other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                         title = "Deleted spots (sample 1)"); s1.spots.plot
ggsave("./project/outcomes/LPS1_deleted.png", plot = s1.spots.plot, scale = 2.5, width = 1920, height = 1080, units = "px")


####--SAMPLE 2--####

# Obtain data
s2 <- createGiottoVisiumObject(
  visium_dir = sample.path[2], #LPS2
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = createGiottoInstructions(save_plot = F, show_plot = F))

## Filter over tissue genes
in.tissue <- pDataDT(s2)[in_tissue == 1]$cell_ID 
if (sum((in.tissue == pDataDT(s2)$cell_ID) == F) == 0) {
  cat("All genes are over tissue")
  rm(in.tissue)
} else {
  s1 <- subsetGiotto(s2, cell_ids = in.tissue)
  rm(in.tissue)
} # All genes are over tissue

s2.matrix <- s2@expression[["cell"]][["rna"]][["raw"]]@exprMat # Extract expression matrix from Giotto object


# Quality control 
## Create cell metadata
reads.depth <- apply(s2.matrix, 2, sum) # Reads per cell
genes.per.cell <- 21675 - apply(s2.matrix == 0, 2, sum) # Detected genes per cell
length(genes.per.cell); length(reads.depth) # Check objects length

is_mito <- grepl("(^MT-)|(^Mt-)|(^mt-)", s2@feat_ID$rna) # Identify mitochondrial genes
if (sum(is_mito == T) != 0) {
  s2@feat_ID$is_mito <- is_mito
  s2.mito <- subsetGiotto(s2, feat_ids = s2@feat_ID$rna[is_mito == T])
  cat(sum(is_mito == T), "mitochondrial genes found:", s2.mito@feat_ID$rna)
  rm(is_mito)
} else {
  cat('No mitochondrial genes found')
  rm(is_mito)
} # 13 mitochondrial genes found

expr.sum <- apply(s2.matrix, 2, sum) # Sum of gene expression values per cell
mito.sum <- apply(s2.matrix[s2.mito@feat_ID$rna, ], 2, sum) # Sum of mito gene expression values per cell
mito.prop <- round(mito.sum/expr.sum, 3) # Mito genes proportion

s2.cell.meta <- data.frame(nReads = reads.depth, nGenes = genes.per.cell, mitoProp = mito.prop) # Write data
rownames(s2.cell.meta) <- s2@cell_ID$cell # Add cell names
rm(reads.depth, genes.per.cell, expr.sum, mito.sum, mito.prop, s2.mito) # Rm variables


## Calcs
### Mito proportion threshold
which(s2.cell.meta$mitoProp > 0.28) # No cells exceed 0.28 mito genes proportion

## Libsize vs Detected genes correlation
with(s2.cell.meta, cor.test(nReads, nGenes)) # Check correlation (0.9314)


## Plots
### Detected genes
s2.detected <- filterDistributions(s2, detection = "cells", nr_bins = 50,
                                   method = "threshold", expression_threshold = 1, 
                                   default_save_name = "detected_genes")

### Genes detection along cells
s2.det.along <- filterDistributions(s2, detection = "feats", nr_bins = 50,
                                    method = "threshold", 
                                    default_save_name = "library_size")

### Libsize vs Detected genes
s2.libvdet <- ggplot(s2.cell.meta, aes(nReads, nGenes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = "Corr = 0.9314", x = 70000, y = 500)

### Threshold evaluation
s2.thresholds <- filterCombinations(s2, expression_thresholds = 1,
                                    feat_det_in_min_cells = c(50, 100, 50, 100),
                                    min_det_feats_per_cell = c(3000, 3000, 4000, 4000), 
                                    show_plot = F,
                                    default_save_name = "thresholds")
s2.thresholds[["ggplot"]][["theme"]][["legend.position"]] <- "none" # Remove legend

s2.metric.plots <- ggarrange(s2.detected, s2.det.along, s2.libvdet, s2.thresholds$ggplot) # Combine all plots
rm(s2.detected, s2.det.along, s2.libvdet, s2.thresholds); s2.metric.plots
ggsave("./project/outcomes/LPS2_metrics.png", plot = s2.metric.plots, scale = 2, width = 1920, height = 1080, units = "px")

## Filter sample
s2.filtered <- filterGiotto(s2, expression_values = "raw", expression_threshold = 1,
                            feat_det_in_min_cells = 50,
                            min_det_feats_per_cell = 3000)
save(s2.filtered, file = "./project/material/filtered_samples/s2_filtered.R") # load("./project/filtered_samples/s2_filtered.R")

## Deleted spots
length(s2@cell_ID$cell) - length(s2.filtered@cell_ID$cell) # 200 deleted spots

### Visualization plot
s2.spots.plot <- spatPlot2D(s2, cell_color = ("lightgrey"), point_size = 2,
                         select_cells = s2.filtered@cell_ID$cell, # Kept spots
                         other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                         title = "Deleted spots (sample 2)"); s2.spots.plot
ggsave("./project/outcomes/LPS2_deleted.png", plot = s2.spots.plot, scale = 2.5, width = 1920, height = 1080, units = "px")


####--SAMPLE 3--####

# Obtain data
s3 <- createGiottoVisiumObject(
  visium_dir = sample.path[3], #LPS3
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = createGiottoInstructions(save_plot = F, show_plot = F))

## Filter over tissue genes
in.tissue <- pDataDT(s3)[in_tissue == 1]$cell_ID 
if (sum((in.tissue == pDataDT(s3)$cell_ID) == F) == 0) {
  cat("All genes are over tissue")
  rm(in.tissue)
} else {
  s1 <- subsetGiotto(s3, cell_ids = in.tissue)
  rm(in.tissue)
} # All genes are over tissue

s3.matrix <- s3@expression[["cell"]][["rna"]][["raw"]]@exprMat # Extract expression matrix from Giotto object


# Quality control 
## Create cell metadata
reads.depth <- apply(s3.matrix, 2, sum) # Reads per cell
genes.per.cell <- 21965 - apply(s3.matrix == 0, 2, sum) # Detected genes per cell
length(genes.per.cell); length(reads.depth) # Check objects length

is_mito <- grepl("(^MT-)|(^Mt-)|(^mt-)", s3@feat_ID$rna) # Identify mitochondrial genes
if (sum(is_mito == T) != 0) {
  s3@feat_ID$is_mito <- is_mito
  s3.mito <- subsetGiotto(s3, feat_ids = s3@feat_ID$rna[is_mito == T])
  cat(sum(is_mito == T), "mitochondrial genes found:", s3.mito@feat_ID$rna)
  rm(is_mito)
} else {
  cat('No mitochondrial genes found')
  rm(is_mito)
} # 13 mitochondrial genes found

expr.sum <- apply(s3.matrix, 2, sum) # Sum of gene expression values per cell
mito.sum <- apply(s3.matrix[s3.mito@feat_ID$rna, ], 2, sum) # Sum of mito gene expression values per cell
mito.prop <- round(mito.sum/expr.sum, 3) # Mito genes proportion

s3.cell.meta <- data.frame(nReads = reads.depth, nGenes = genes.per.cell, mitoProp = mito.prop) # Write data
rownames(s3.cell.meta) <- s3@cell_ID$cell # Add cell names
rm(reads.depth, genes.per.cell, expr.sum, mito.sum, mito.prop, s3.mito) # Rm variables


## Calcs
### Mito proportion threshold
which(s3.cell.meta$mitoProp > 0.28) # No cells exceed 0.28 mito genes proportion

## Libsize vs Detected genes correlation
with(s3.cell.meta, cor.test(nReads, nGenes)) # Check correlation (0.9584)


## Plots
### Detected genes
s3.detected <- filterDistributions(s3, detection = "cells", nr_bins = 50,
                                   method = "threshold", expression_threshold = 1, 
                                   default_save_name = "detected_genes")

### Genes detection along cells
s3.det.along <- filterDistributions(s3, detection = "feats", nr_bins = 50,
                                    method = "threshold", 
                                    default_save_name = "library_size")

### Libsize vs Detected genes
s3.libvdet <- ggplot(s3.cell.meta, aes(nReads, nGenes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = "Corr = 0.9584", x = 60000, y = 500)

### Threshold evaluation
s3.thresholds <- filterCombinations(s3, expression_thresholds = 1,
                                    feat_det_in_min_cells = c(50, 100, 50, 100),
                                    min_det_feats_per_cell = c(2500, 2500, 4000, 4000), 
                                    show_plot = F,
                                    default_save_name = "thresholds")
s3.thresholds[["ggplot"]][["theme"]][["legend.position"]] <- "none" # Remove legend

s3.metric.plots <- ggarrange(s3.detected, s3.det.along, s3.libvdet, s3.thresholds$ggplot) # Combine all plots
rm(s3.detected, s3.det.along, s3.libvdet, s3.thresholds); s3.metric.plots
ggsave("./project/outcomes/LPS3_metrics.png", plot = s3.metric.plots, scale = 2, width = 1920, height = 1080, units = "px")

## Filter sample
s3.filtered <- filterGiotto(s3, expression_values = "raw", expression_threshold = 1,
                            feat_det_in_min_cells = 50,
                            min_det_feats_per_cell = 2500)
save(s3.filtered, file = "./project/material/filtered_samples/s3_filtered.R") # load("./project/filtered_samples/s3_filtered.R")

## Deleted spots
length(s3@cell_ID$cell) - length(s3.filtered@cell_ID$cell) # 194 deleted spots

### Visualization plot
s3.spots.plot <- spatPlot2D(s3, cell_color = ("lightgrey"), point_size = 2,
                            select_cells = s3.filtered@cell_ID$cell, # Kept spots
                            other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                            title = "Deleted spots (sample 3)"); s3.spots.plot
ggsave("./project/outcomes/LPS3_deleted.png", plot = s3.spots.plot, scale = 2.5, width = 1920, height = 1080, units = "px")


####--SAMPLE 4--####

# Obtain data
s4 <- createGiottoVisiumObject(
  visium_dir = sample.path[4], #Saline1
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = createGiottoInstructions(save_plot = F, show_plot = F))

## Filter over tissue genes
in.tissue <- pDataDT(s4)[in_tissue == 1]$cell_ID 
if (sum((in.tissue == pDataDT(s4)$cell_ID) == F) == 0) {
  cat("All genes are over tissue")
  rm(in.tissue)
} else {
  s1 <- subsetGiotto(s4, cell_ids = in.tissue)
  rm(in.tissue)
} # All genes are over tissue

s4.matrix <- s4@expression[["cell"]][["rna"]][["raw"]]@exprMat # Extract expression matrix from Giotto object


# Quality control 
## Create cell metadata
reads.depth <- apply(s4.matrix, 2, sum) # Reads per cell
genes.per.cell <- 21352 - apply(s4.matrix == 0, 2, sum) # Detected genes per cell
length(genes.per.cell); length(reads.depth) # Check objects length

is_mito <- grepl("(^MT-)|(^Mt-)|(^mt-)", s4@feat_ID$rna) # Identify mitochondrial genes
if (sum(is_mito == T) != 0) {
  s4@feat_ID$is_mito <- is_mito
  s4.mito <- subsetGiotto(s4, feat_ids = s4@feat_ID$rna[is_mito == T])
  cat(sum(is_mito == T), "mitochondrial genes found:", s4.mito@feat_ID$rna)
  rm(is_mito)
} else {
  cat('No mitochondrial genes found')
  rm(is_mito)
} # 13 mitochondrial genes found

expr.sum <- apply(s4.matrix, 2, sum) # Sum of gene expression values per cell
mito.sum <- apply(s4.matrix[s4.mito@feat_ID$rna, ], 2, sum) # Sum of mito gene expression values per cell
mito.prop <- round(mito.sum/expr.sum, 3) # Mito genes proportion

s4.cell.meta <- data.frame(nReads = reads.depth, nGenes = genes.per.cell, mitoProp = mito.prop) # Write data
rownames(s4.cell.meta) <- s4@cell_ID$cell # Add cell names
rm(reads.depth, genes.per.cell, expr.sum, mito.sum, mito.prop, s4.mito) # Rm variables


## Calcs
### Mito proportion threshold
which(s4.cell.meta$mitoProp > 0.28) # No cells exceed 0.28 mito genes proportion

## Libsize vs Detected genes correlation
with(s4.cell.meta, cor.test(nReads, nGenes)) # Check correlation (0.9574)


## Plots
### Detected genes
s4.detected <- filterDistributions(s4, detection = "cells", nr_bins = 50,
                                   method = "threshold", expression_threshold = 1, 
                                   default_save_name = "detected_genes")

### Genes detection along cells
s4.det.along <- filterDistributions(s4, detection = "feats", nr_bins = 50,
                                    method = "threshold", 
                                    default_save_name = "library_size")

### Libsize vs Detected genes
s4.libvdet <- ggplot(s4.cell.meta, aes(nReads, nGenes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = "Corr = 0.9574", x = 60000, y = 500)

### Threshold evaluation
s4.thresholds <- filterCombinations(s4, expression_thresholds = 1,
                                    feat_det_in_min_cells = c(50, 100, 50, 100),
                                    min_det_feats_per_cell = c(3000, 3000, 4000, 4000), 
                                    show_plot = F,
                                    default_save_name = "thresholds")
s4.thresholds[["ggplot"]][["theme"]][["legend.position"]] <- "none" # Remove legend

s4.metric.plots <- ggarrange(s4.detected, s4.det.along, s4.libvdet, s4.thresholds$ggplot) # Combine all plots
rm(s4.detected, s4.det.along, s4.libvdet, s4.thresholds); s4.metric.plots
ggsave("./project/outcomes/Saline1_metrics.png", plot = s4.metric.plots, scale = 2, width = 1920, height = 1080, units = "px")

## Filter sample
s4.filtered <- filterGiotto(s4, expression_values = "raw", expression_threshold = 1,
                            feat_det_in_min_cells = 50,
                            min_det_feats_per_cell = 3000)
save(s4.filtered, file = "./project/material/filtered_samples/s4_filtered.R") # load("./project/filtered_samples/s4_filtered.R")

## Deleted spots
length(s4@cell_ID$cell) - length(s4.filtered@cell_ID$cell) # 141 deleted spots

### Visualization plot
s4.spots.plot <- spatPlot2D(s4, cell_color = ("lightgrey"), point_size = 2,
                            select_cells = s4.filtered@cell_ID$cell, # Kept spots
                            other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                            title = "Deleted spots (sample 4)"); s4.spots.plot
ggsave("./project/outcomes/Saline1_deleted.png", plot = s4.spots.plot, scale = 2.5, width = 1920, height = 1080, units = "px")


####--SAMPLE 5--####

# Obtain data
s5 <- createGiottoVisiumObject(
  visium_dir = sample.path[5], #Saline2
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = createGiottoInstructions(save_plot = F, show_plot = F))

## Filter over tissue genes
in.tissue <- pDataDT(s5)[in_tissue == 1]$cell_ID 
if (sum((in.tissue == pDataDT(s5)$cell_ID) == F) == 0) {
  cat("All genes are over tissue")
  rm(in.tissue)
} else {
  s1 <- subsetGiotto(s5, cell_ids = in.tissue)
  rm(in.tissue)
} # All genes are over tissue

s5.matrix <- s5@expression[["cell"]][["rna"]][["raw"]]@exprMat # Extract expression matrix from Giotto object


# Quality control 
## Create cell metadata
reads.depth <- apply(s5.matrix, 2, sum) # Reads per cell
genes.per.cell <- 20975 - apply(s5.matrix == 0, 2, sum) # Detected genes per cell
length(genes.per.cell); length(reads.depth) # Check objects length

is_mito <- grepl("(^MT-)|(^Mt-)|(^mt-)", s5@feat_ID$rna) # Identify mitochondrial genes
if (sum(is_mito == T) != 0) {
  s5@feat_ID$is_mito <- is_mito
  s5.mito <- subsetGiotto(s5, feat_ids = s5@feat_ID$rna[is_mito == T])
  cat(sum(is_mito == T), "mitochondrial genes found:", s5.mito@feat_ID$rna)
  rm(is_mito)
} else {
  cat('No mitochondrial genes found')
  rm(is_mito)
} # 13 mitochondrial genes found

expr.sum <- apply(s5.matrix, 2, sum) # Sum of gene expression values per cell
mito.sum <- apply(s5.matrix[s5.mito@feat_ID$rna, ], 2, sum) # Sum of mito gene expression values per cell
mito.prop <- round(mito.sum/expr.sum, 3) # Mito genes proportion

s5.cell.meta <- data.frame(nReads = reads.depth, nGenes = genes.per.cell, mitoProp = mito.prop) # Write data
rownames(s5.cell.meta) <- s5@cell_ID$cell # Add cell names
rm(reads.depth, genes.per.cell, expr.sum, mito.sum, mito.prop, s5.mito) # Rm variables


## Calcs
### Mito proportion threshold
which(s5.cell.meta$mitoProp > 0.28) # No cells exceed 0.28 mito genes proportion

## Libsize vs Detected genes correlation
with(s5.cell.meta, cor.test(nReads, nGenes)) # Check correlation (0.93322)


## Plots
### Detected genes
s5.detected <- filterDistributions(s5, detection = "cells", nr_bins = 50,
                                   method = "threshold", expression_threshold = 1, 
                                   default_save_name = "detected_genes")

### Genes detection along cells
s5.det.along <- filterDistributions(s5, detection = "feats", nr_bins = 50,
                                    method = "threshold", 
                                    default_save_name = "library_size")

### Libsize vs Detected genes
s5.libvdet <- ggplot(s5.cell.meta, aes(nReads, nGenes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = "Corr = 0.93322", x = 60000, y = 500)

### Threshold evaluation
s5.thresholds <- filterCombinations(s5, expression_thresholds = 1,
                                    feat_det_in_min_cells = c(50, 100, 50, 100),
                                    min_det_feats_per_cell = c(150, 150, 250, 250), 
                                    show_plot = F,
                                    default_save_name = "thresholds")
s5.thresholds[["ggplot"]][["theme"]][["legend.position"]] <- "none" # Remove legend

s5.metric.plots <- ggarrange(s5.detected, s5.det.along, s5.libvdet, s5.thresholds$ggplot) # Combine all plots
rm(s5.detected, s5.det.along, s5.libvdet, s5.thresholds); s5.metric.plots
ggsave("./project/outcomes/Saline2_metrics.png", plot = s5.metric.plots, scale = 2, width = 1920, height = 1080, units = "px")

## Filter sample
s5.filtered <- filterGiotto(s5, expression_values = "raw", expression_threshold = 1,
                            feat_det_in_min_cells = 50,
                            min_det_feats_per_cell = 150)
save(s5.filtered, file = "./project/material/filtered_samples/s5_filtered.R") # load("./project/filtered_samples/s5_filtered.R")

## Deleted spots
length(s5@cell_ID$cell) - length(s5.filtered@cell_ID$cell) # 148 deleted spots

### Visualization plot
s5.spots.plot <- spatPlot2D(s5, cell_color = ("lightgrey"), point_size = 2,
                            select_cells = s5.filtered@cell_ID$cell, # Kept spots
                            other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                            title = "Deleted spots (sample 5)"); s5.spots.plot
ggsave("./project/outcomes/Saline2_deleted.png", plot = s5.spots.plot, scale = 2.5, width = 1920, height = 1080, units = "px")


####--SAMPLE 6--####

# Obtain data
s6 <- createGiottoVisiumObject(
  visium_dir = sample.path[6], #Saline6
  gene_column_index = 2, # Use gene symbols
  expr_data = "filter", # Use filtered data
  png_name = "tissue_lowres_image.png", # Select lowres image
  instructions = createGiottoInstructions(save_plot = F, show_plot = F))

## Filter over tissue genes
in.tissue <- pDataDT(s6)[in_tissue == 1]$cell_ID 
if (sum((in.tissue == pDataDT(s6)$cell_ID) == F) == 0) {
  cat("All genes are over tissue")
  rm(in.tissue)
} else {
  s1 <- subsetGiotto(s6, cell_ids = in.tissue)
  rm(in.tissue)
} # All genes are over tissue

s6.matrix <- s6@expression[["cell"]][["rna"]][["raw"]]@exprMat # Extract expression matrix from Giotto object


# Quality control 
## Create cell metadata
reads.depth <- apply(s6.matrix, 2, sum) # Reads per cell
genes.per.cell <- 21714 - apply(s6.matrix == 0, 2, sum) # Detected genes per cell
length(genes.per.cell); length(reads.depth) # Check objects length

is_mito <- grepl("(^MT-)|(^Mt-)|(^mt-)", s6@feat_ID$rna) # Identify mitochondrial genes
if (sum(is_mito == T) != 0) {
  s6@feat_ID$is_mito <- is_mito
  s6.mito <- subsetGiotto(s6, feat_ids = s6@feat_ID$rna[is_mito == T])
  cat(sum(is_mito == T), "mitochondrial genes found:", s6.mito@feat_ID$rna)
  rm(is_mito)
} else {
  cat('No mitochondrial genes found')
  rm(is_mito)
} # 13 mitochondrial genes found

expr.sum <- apply(s6.matrix, 2, sum) # Sum of gene expression values per cell
mito.sum <- apply(s6.matrix[s6.mito@feat_ID$rna, ], 2, sum) # Sum of mito gene expression values per cell
mito.prop <- round(mito.sum/expr.sum, 3) # Mito genes proportion

s6.cell.meta <- data.frame(nReads = reads.depth, nGenes = genes.per.cell, mitoProp = mito.prop) # Write data
rownames(s6.cell.meta) <- s6@cell_ID$cell # Add cell names
rm(reads.depth, genes.per.cell, expr.sum, mito.sum, mito.prop, s6.mito) # Rm variables


## Calcs
### Mito proportion threshold
which(s6.cell.meta$mitoProp > 0.28) # No cells exceed 0.28 mito genes proportion

## Libsize vs Detected genes correlation
with(s6.cell.meta, cor.test(nReads, nGenes)) # Check correlation (0.9272)


## Plots
### Detected genes
s6.detected <- filterDistributions(s6, detection = "cells", nr_bins = 50,
                                   method = "threshold", expression_threshold = 1, 
                                   default_save_name = "detected_genes")

### Genes detection along cells
s6.det.along <- filterDistributions(s6, detection = "feats", nr_bins = 50,
                                    method = "threshold", 
                                    default_save_name = "library_size")

### Libsize vs Detected genes
s6.libvdet <- ggplot(s6.cell.meta, aes(nReads, nGenes)) + geom_point(aes(alpha = 0.3)) +
  theme_classic() + theme(legend.position = "none") + 
  xlab("total reads per spot") + ylab("detected genes per spot") +
  annotate("text", label = "Corr = 0.9272", x = 95000, y = 500)

### Threshold evaluation
s6.thresholds <- filterCombinations(s6, expression_thresholds = 1,
                                    feat_det_in_min_cells = c(50, 100, 50, 100),
                                    min_det_feats_per_cell = c(1500, 1500, 2000, 2000), 
                                    show_plot = F,
                                    default_save_name = "thresholds")
s6.thresholds[["ggplot"]][["theme"]][["legend.position"]] <- "none" # Remove legend

s6.metric.plots <- ggarrange(s6.detected, s6.det.along, s6.libvdet, s6.thresholds$ggplot) # Combine all plots
rm(s6.detected, s6.det.along, s6.libvdet, s6.thresholds); s6.metric.plots
ggsave("./project/outcomes/Saline3_metrics.png", plot = s6.metric.plots, scale = 2, width = 1920, height = 1080, units = "px")

## Filter sample
s6.filtered <- filterGiotto(s6, expression_values = "raw", expression_threshold = 1,
                            feat_det_in_min_cells = 50,
                            min_det_feats_per_cell = 1500)
save(s6.filtered, file = "./project/material/filtered_samples/s6_filtered.R") # load("./project/filtered_samples/s6_filtered.R")

## Deleted spots
length(s6@cell_ID$cell) - length(s6.filtered@cell_ID$cell) # 142 deleted spots

### Visualization plot
s6.spots.plot <- spatPlot2D(s6, cell_color = ("lightgrey"), point_size = 2,
                            select_cells = s6.filtered@cell_ID$cell, # Kept spots
                            other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                            title = "Deleted spots (sample 6)"); s6.spots.plot
ggsave("./project/outcomes/Saline3_deleted.png", plot = s6.spots.plot, scale = 2.5, width = 1920, height = 1080, units = "px")

