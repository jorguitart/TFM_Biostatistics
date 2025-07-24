library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

df <- DEA.spatialGE

# Separar genes sobreexpresados e infraexpresados
genes_over <- df$gene[df$sig == "over"]
genes_under <- df$gene[df$sig == "under"]

# Convertir símbolos de genes a IDs de Entrez
gene_ids_over <- bitr(genes_over, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
gene_ids_under <- bitr(genes_under, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

# Enriquecimiento GO para genes sobreexpresados
ego_over <- enrichGO(gene         = gene_ids_over$ENTREZID,
                     OrgDb        = org.Hs.eg.db,
                     keyType      = "ENTREZID",
                     ont          = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05)

# Enriquecimiento GO para genes infraexpresados
ego_under <- enrichGO(gene         = gene_ids_under$ENTREZID,
                      OrgDb        = org.Hs.eg.db,
                      keyType      = "ENTREZID",
                      ont          = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05)

# Dotplot
dotplot(ego_over, showCategory = 15, title = "GO enrichment - Overexpressed genes")
dotplot(ego_under, showCategory = 15, title = "GO enrichment - Underexpressed genes")

# Network plot
cnetplot(ego_over, categorySize = "pvalue", showCategory = 5)
cnetplot(ego_under, categorySize = "pvalue", showCategory = 5)
