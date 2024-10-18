Finding markers for each cell cluster


# Change the default assay
dat<-readRDS("Cluster1.rds")

DefaultAssay(object = dat) <- "RNA"

markers_genes <- FindAllMarkers(dat, log2FC.threshold = 0.2, test.use = "wilcox",min.pct = 0.1, min.diff.pct = 0.2, only.pos = TRUE, max.cells.per.ident = 50, assay = "RNA")

markers_genes %>%
    group_by(cluster) %>%
    top_n(-25, p_val_adj) -> top25
top25

markers_genes %>%
    group_by(cluster) %>%
    top_n(-5, p_val_adj) -> top5
top5


mypar(3, 9, mar = c(4, 6, 3, 1))
for (i in unique(top25$cluster)) {
    barplot(sort(setNames(top25$avg_log2FC, top25$gene)[top25$cluster == i], F),
        horiz = T, las = 1, main = paste0(i, " vs. rest"), border = "white", yaxs = "i")
    abline(v = c(0, 0.25), lty = c(1, 2))
}

results_top25<-top25 %>% group_by(cluster,gene)%>% summarize(N=n())
write.csv(results_top25,"Top25_Genes_Clusters")

DoHeatmap(dat, features = as.character(unique(top5$gene)), group.by = "seurat_clusters", assay = "RNA")
 
dat <- ScaleData(dat, features = as.character(unique(top5$gene)), assay = "RNA")
DoHeatmap(dat, features = as.character(unique(top5$gene)), group.by = "seurat_clusters", assay = "RNA")


# take top 3 genes per cluster/
top5 %>%
    group_by(cluster) %>%
    top_n(-2, p_val) -> top2


# set pt.size to zero if you do not want all the points to hide the violin
# shapes, or to a small value like 0.1
VlnPlot(dat, features = as.character(unique(top2$gene)), ncol = 5, group.by = "seurat_clusters", assay = "RNA", pt.size = 0)
    
cons0 <- FindConservedMarkers(dat, ident.1 = 0, only.pos = TRUE, logfc.threshold = 0.25, grouping.var = "sample")
cons1 <- FindConservedMarkers(dat, ident.1 = 1, only.pos = TRUE, logfc.threshold = 0.25, grouping.var = "sample" )
cons2 <- FindConservedMarkers(dat, ident.1 = 2, only.pos = TRUE, logfc.threshold = 0.25, grouping.var = "sample")
cons3 <- FindConservedMarkers(dat, ident.1 = 3, only.pos = TRUE, logfc.threshold = 0.25, grouping.var = "sample")
cons4 <- FindConservedMarkers(dat, ident.1 = 4, only.pos = TRUE, logfc.threshold = 0.25, grouping.var = "sample")
cons5 <- FindConservedMarkers(dat, ident.1 = 5, only.pos = TRUE, logfc.threshold = 0.25, grouping.var = "sample")
cons6 <- FindConservedMarkers(dat, ident.1 = 6, only.pos = TRUE, logfc.threshold = 0.25, grouping.var = "sample")
cons7 <- FindConservedMarkers(dat, ident.1 = 7, only.pos = TRUE, logfc.threshold = 0.25, grouping.var = "sample")
cons8 <- FindConservedMarkers(dat, ident.1 = 8, only.pos = TRUE, logfc.threshold = 0.25, grouping.var = "sample")
cons9 <- FindConservedMarkers(dat, ident.1 = 9, only.pos = TRUE, logfc.threshold = 0.25, grouping.var = "sample")
cons10 <- FindConservedMarkers(dat, ident.1 = 10, only.pos = TRUE, logfc.threshold = 0.25, grouping.var = "sample")
cons11 <- FindConservedMarkers(dat, ident.1 = 11, only.pos = TRUE, logfc.threshold = 0.25, grouping.var = "sample")
cons12 <- FindConservedMarkers(dat, ident.1 = 12, only.pos = TRUE, logfc.threshold = 0.25, grouping.var = "sample")
cons13 <- FindConservedMarkers(dat, ident.1 = 13, only.pos = TRUE, logfc.threshold = 0.25, grouping.var = "sample")

write.csv(cons0, file = "Cons0.csv")
write.csv(cons1, file = "Cons1.csv")
write.csv(cons2, file = "Cons2.csv")
write.csv(cons3, file = "Cons3.csv")
write.csv(cons4, file = "Cons4.csv")
write.csv(cons5, file = "Cons5.csv")
write.csv(cons6, file = "Cons6.csv")
write.csv(cons7, file = "Cons7.csv")
write.csv(cons8, file = "Cons8.csv")
write.csv(cons9, file = "Cons9.csv")
write.csv(cons10, file = "Cons10.csv")
write.csv(cons11, file = "Cons11.csv")
write.csv(cons12, file = "Cons12.csv")
write.csv(cons13, file = "Cons13.csv")
