####--LIBRARIES--####
library(tidyverse)
library(ggpubr)
library(Giotto)
library(spatialGE)
library(ggrepel)
library(patchwork)

####--DATA--####
Giotto <- loadGiotto("./project/material/Giotto/HMRF_sample", python_path = "C:/ProgramData/anaconda3/python.exe")
load("./project/material/Giotto/DEA.RData"); DEA.Giotto <- DEA; rm(DEA)
load("./project/material/spatialGE/newsGEobject.RData"); spatialGE <- newsGE.obj; rm(newsGE.obj)
load("./project/material/spatialGE/DEAt.RData"); DEA.spatialGE.t <- DEA.t; rm(DEA.t)
load("./project/material/spatialGE/DEAw.RData"); DEA.spatialGE.w <- DEA.w; rm(DEA.w)

####--VISUALIZATION--####
# Giotto
## HVF
calculateHVF(Giotto, expression_values = "normalized",
             method = "cov_loess", show_plot = T)

## Dim reduction
screePlot(Giotto)

plotUMAP(Giotto, cell_color = "type", point_alpha = 0.3, point_size = 1.2, 
         cell_color_code = getRainbowColors(3, slim = c(0.1, 1), vlim = 0.8, seed = 123),
         show_center_label = F, title = "")

## Clustering
spatPlot(Giotto, group_by = "list_ID", cell_color = "leiden_clus")

## Spatial domains
spatPlot(Giotto, group_by = "list_ID", cell_color = "HMRF_B20")

## Regions
spatPlot(Giotto, group_by = "list_ID", cell_color = "domain")

## DEA
DEA.Giotto.plot <- ggplot(DEA.Giotto, aes(x = summary.logFC, y = logFDR, color = sig)) +
  geom_point(alpha = 0.8) +
  geom_text_repel(aes(label = label), col = "black", size = 3, max.overlaps = 15, alpha = 0.8) +
  scale_color_manual(values = c("under" = "#FFC98B", "over" = "#CCA7FF", "nosig" = "grey")) +
  theme_bw() + theme(legend.position = "") +
  labs(x = expression(Log[2]~Fold~Change), y = expression(-~Log[10]~FDR)) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed"); DEA.Giotto.plot

# spatialGE
## Clusters (sp. domains)
STplot(spatialGE, plot_meta = "cluster", ptsize = 1.2)

## Domains (regions)
STplot(spatialGE, plot_meta = "domain", ptsize = 1.2)

## DEA
DEA.spatialGE <- DEA.spatialGE.t$MSCA

minim <- min(DEA.spatialGE$adj_p_val[DEA.spatialGE$adj_p_val > 0], na.rm = T)
DEA.spatialGE$FDR <- ifelse(DEA.spatialGE$adj_p_val == 0, minim, DEA.spatialGE$adj_p_val)
DEA.spatialGE$logFDR <- -log10(DEA.spatialGE$FDR)
DEG <- DEA.spatialGE$gene[which(DEA.spatialGE$FDR < 0.05 & abs(DEA.spatialGE$avg_log2fc) > 0.5)]; length(DEG)
colors <- ifelse(DEA.spatialGE$FDR < 0.05 & DEA.spatialGE$avg_log2fc < -0.5, "under",
                 ifelse(DEA.spatialGE$FDR < 0.05 & DEA.spatialGE$avg_log2fc > 0.5, "over", "nosig"))
names(colors) <- DEA.spatialGE$gene
DEA.spatialGE$sig <- colors
DEA.spatialGE$label <- ifelse(DEA.spatialGE$sig %in% c("under", "over"), DEA.spatialGE$gene, NA)

DEA.spatialGE.plot <- ggplot(DEA.spatialGE[DEA.spatialGE$cluster_1 == "LC", ], 
                             aes(x = avg_log2fc, y = logFDR, color = sig)) +
  geom_point(alpha = 0.8) +
  geom_text_repel(aes(label = label), col = "black", size = 3, max.overlaps = 15, alpha = 0.8) +
  scale_color_manual(values = c("under" = "#FFC98B", "over" = "#CCA7FF", "nosig" = "grey")) +
  theme_bw() + theme(legend.position = "") +
  labs(x = expression(Log[2]~Fold~Change), y = expression(-~Log[10]~FDR)) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed"); DEA.spatialGE.plot

ggsave(plot = DEA.spatialGE.plot, filename = "./project/outcomes/vis/LCvsPPWM_spatialGE.png", 
       width = 1920, height = 1080, units = "px", scale = 2)
