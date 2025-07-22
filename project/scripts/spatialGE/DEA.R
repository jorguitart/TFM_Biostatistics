####--LIBRARIES--####
library(tidyverse)
library(ggpubr)
# library(Giotto) #pak::pkg_install("drieslab/Giotto")
# library(spacexr) #devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
# library(nicheDE) #devtools::install_github("kaishumason/NicheDE")
library(spatialGE) #devtools::install_github("fridleylab/spatialGE")
library(ggrepel)

setwd("~/TFM")

####--DATA--####
load("./project/material/spatialGE/typesGEobject.RData")

DEA <- STdiff(typesGE.obj, samples = "MSCA", annot = "domain", clusters = c("LC", "PPWM"), 
              topgenes = 11534, cores = 24)

DEA <- DEA$MSCA
minim <- min(DEA$adj_p_val[DEA$adj_p_val > 0])
DEA$FDR <- ifelse(DEA$adj_p_val == 0, minim, DEA$adj_p_val)
DEA$logFDR <- -log10(DEA$adj_p_val)
colors <- ifelse(DEA$FDR < 0.05 & DEA$avg_log2fc < -0.5, "under",
                 ifelse(DEA$FDR < 0.05 & DEA$avg_log2fc > 0.5, "over", "nosig"))
names(colors) <- DEA$gene
DEA$sig <- colors
DEA$label <- ifelse(DEA$sig %in% c("under", "over"), DEA$gene, NA)

save(DEA, file = "./project/material/spatialGE/DEA.RData")

DEA.plot <- ggplot(DEA[DEA$cluster_1 == "LC", ], aes(x = avg_log2fc, y = logFDR, color = sig)) +
  geom_point(alpha = 0.8) +
  geom_text_repel(aes(label = label), col = "black", size = 3, max.overlaps = 15, alpha = 0.8) +
  scale_color_manual(values = c("under" = "#FFC98B", "over" = "#CCA7FF", "nosig" = "grey")) +
  theme_bw() + theme(legend.position = "") +
  labs(x = expression(Log[2]~Fold~Change), y = expression(-~Log[10]~FDR)) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")
DEA.plot

ggsave(plot = DEA.plot2, filename = "./project/outcomes/vis/LCvsPPWM_sGE.png", 
       width = 1920, height = 1080, units = "px", scale = 2)
