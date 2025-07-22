####--LIBRARIES--####
library(tidyverse)
library(ggpubr)
library(Giotto) #pak::pkg_install("drieslab/Giotto")
# library(spacexr) #devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
# library(nicheDE) #devtools::install_github("kaishumason/NicheDE")
# library(spatialGE) #devtools::install_github("fridleylab/spatialGE")
library(ggrepel)

setwd("~/TFM")

####--DATA--####
sample <- loadGiotto(path_to_folder = "./project/material/Giotto/HMRF_sample",
                     python_path = "C:/ProgramData/anaconda3/python.exe")

# Lesion core between conditions (Inactive VS Active)
LC <- subsetGiotto(sample, 
                   cell_ids = pDataDT(sample)[type %in% c("MSCA", "MSCI") & domain == "LC"]$cell_ID)

DEA <- findMarkers(LC,
                   expression_values = "normalized",
                   cluster_column = "type",
                   method = "scran")
DEA <- DEA[[2]]
minim <- min(DEA$FDR[DEA$FDR > 0])
DEA$FDR <- ifelse(DEA$FDR == 0, minim + runif(sum(DEA$FDR == 0), min = 0.1 * minim, max = -0.1 * minim),
                  DEA$FDR)
DEA$logFDR <- -log10(DEA$FDR)
DEG <- DEA[which(FDR < 0.05 & abs(logFC.MSCA) > 0.5)]; nrow(DEG)
colors <- ifelse(DEA$FDR < 0.05 & DEA$summary.logFC < -0.5, "under",
                 ifelse(DEA$FDR < 0.05 & DEA$summary.logFC > 0.5, "over", "nosig"))
names(colors) <- DEA$feats
DEA$sig <- colors
DEA$label <- ifelse(DEA$sig %in% c("under", "over"), DEA$feats, NA)


DEA.plot <- ggplot(DEA, aes(x = summary.logFC, y = logFDR, color = sig)) +
  geom_point(alpha = 0.8) +
  geom_text_repel(aes(label = label), col = "black", size = 3, max.overlaps = 15, alpha = 0.8) +
  scale_color_manual(values = c("under" = "#FFC98B", "over" = "#CCA7FF", "nosig" = "grey")) +
  theme_bw() + theme(legend.position = "") +
  labs(x = expression(Log[2]~Fold~Change), y = expression(-~Log[10]~FDR)) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")

ggsave(plot = DEA.plot, filename = "./project/outcomes/vis/MSCAvsMSCI.png", 
       width = 1920, height = 1080, units = "px", scale = 2)

# MSCA between domains (Lesion core vs )
MSCA <- subsetGiotto(sample, 
                     cell_ids = pDataDT(sample)[type == "MSCA" & domain %in% c("LC", "PPWM")]$cell_ID)

DEA <- findMarkers(MSCA,
                   expression_values = "normalized",
                   cluster_column = "domain",
                   method = "scran")
DEA <- DEA[[1]]
minim <- min(DEA$FDR[DEA$FDR > 0])
DEA$FDR <- ifelse(DEA$FDR == 0, minim + runif(sum(DEA$FDR == 0), min = 0.1 * minim, max = -0.1 * minim),
                  DEA$FDR)
DEA$logFDR <- -log10(DEA$FDR)
nrow(DEA[which(FDR < 0.05 & abs(logFC.PPWM) > 0.5)])
colors <- ifelse(DEA$FDR < 0.05 & DEA$summary.logFC < -0.5, "under",
                 ifelse(DEA$FDR < 0.05 & DEA$summary.logFC > 0.5, "over", "nosig"))
names(colors) <- DEA$feats
DEA$sig <- colors
DEA$label <- ifelse(DEA$sig %in% c("under", "over"), DEA$feats, NA)


DEA.plot2 <- ggplot(DEA, aes(x = summary.logFC, y = logFDR, color = sig)) +
  geom_point(alpha = 0.8) +
  geom_text_repel(aes(label = label), col = "black", size = 3, max.overlaps = 15, alpha = 0.8) +
  scale_color_manual(values = c("under" = "#FFC98B", "over" = "#CCA7FF", "nosig" = "grey")) +
  theme_bw() + theme(legend.position = "") +
  labs(x = expression(Log[2]~Fold~Change), y = expression(-~Log[10]~FDR)) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")

ggsave(plot = DEA.plot2, filename = "./project/outcomes/vis/LCvsPPWM_Giotto.png", 
       width = 1920, height = 1080, units = "px", scale = 2)
