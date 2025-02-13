####--LIBRARIES--####

library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(spacexr) #devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
library(nicheDE) #devtools::install_github("kaishumason/NicheDE")
library(spatialGE) #devtools::install_github("fridleylab/spatialGE")
library(Giotto) #pak::pkg_install("drieslab/Giotto")


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


# Visualize spots
s1.spots <- spatPlot(s1, point_alpha = 0.6, default_save_name = "spots") # Save for later


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
save(s1.filtered, file = "./project/filtered_samples/s1_filtered.R") # load("./project/filtered_samples/s1_filtered.R")

## Deleted spots
length(s1@cell_ID$cell) - length(s1.filtered@cell_ID$cell) # 243 deleted spots

### Visualization plot
spots.plot <- spatPlot2D(s1, cell_color = ("lightgrey"), point_size = 2,
                         select_cells = s1.filtered@cell_ID$cell, # Select kept spots
                         other_cell_color = "red3", other_point_size = 2, # Mark deleted spots
                         title = "Deleted spots (sample 1)"); spots.plot
ggsave("./project/outcomes/LPS1_deleted.png", plot = spots.plot, scale = 2.5, width = 1920, height = 1080, units = "px")


###--SAMPLE 2--####