#Load R libraries
library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(hdf5r)
library(presto)
library(GenomeInfoDb)
library(stringr)
library(ggplot2)
library(RColorBrewer)


library(BSgenome.Ggallus.UCSC.galGal4)
gtf <- rtracklayer::import('/proj/sllstore2017078/private/martinj/single_cell_pilot/annotation/Gallus_gallus.GRCg6a.104_filtered.gtf', format="gtf")
gene.coords <- gtf[gtf$type == 'gene' ]


dat<-readRDS("final_Merged.QC.rds")

#Run MACS2
DefaultAssay(dat)<-"ATAC"
peaks <- CallPeaks(dat, macs2.path = "/sw/bioinfo/MACS/2.1.0/rackham/bin/macs2")


# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")

DefaultAssay(dat) <- "ATAC"
saveRDS(peaks,file="combined.peakset.rds")

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(dat),
  features = peaks,
  cells = colnames(dat)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
dat[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = dat@assays$ATAC@fragments,
  annotation = gene.coords
)

#Add colors
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#pie(rep(1,n), col=sample(col_vector, n))

#my_cols = brewer.pal(8,"Dark2 ")
alpha_val=0.33
#RNA Processing
DefaultAssay(dat) <- "RNA"
dat <- SCTransform(dat)
dat <- RunPCA(dat)
dat<- RunUMAP(
  object = dat,
  reduction.name="rna_umap",
  reduction="pca",
  assay = "SCT",
  verbose = TRUE,
  dims=1:50
)
p1<-DimPlot(dat,reduction="rna_umap",group.by="sample",cols=alpha(col_vector,alpha_val))+ggtitle("RNA UMAP")

#p1<-DimPlot(dat,reduction="rna_umap",cols=alpha(my_cols,alpha_val))+ggtitle("RNA UMAP")



#DNA Accessibility processing
DefaultAssay(dat) <- "peaks"
dat <- FindTopFeatures(dat, min.cutoff = 5)
dat <- RunTFIDF(dat)
dat <- RunSVD(dat)
dat<- RunUMAP(
  object = dat,
  reduction.name="atac_umap",
  reduction="lsi",
  assay = "peaks",
  verbose = TRUE,
  dims=2:40
)
p2<-DimPlot(dat,reduction="atac_umap",group.by="sample",cols=alpha(col_vector,alpha_val))+ggtitle("ATAC UMAP")


# build a joint neighbor graph using both assays
dat <- FindMultiModalNeighbors(
  object = dat,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40), #I think the ATAC UMAP does a better job integrating samples, maybe skip dim 1 for RNA also?
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# build a joint UMAP visualization
dat <- RunUMAP(
  object = dat,
  nn.name = "weighted.nn",
  reduction.name="multimodal_umap",
  assay = "RNA",
  verbose = TRUE
)
p3<-DimPlot(dat,reduction="multimodal_umap",group.by="sample",cols=alpha(col_vector,alpha_val))+ggtitle("Multimodal UMAP")

#Cluster on multimodal graph
dat <- FindClusters(dat, resolution = 0.8, verbose = FALSE,graph="wknn")
#p4<-DimPlot(dat,reduction="multimodal_umap",group.by="seurat_clusters")+ggtitle("Multimodal UMAP Clusters")
p4<-DimPlot(dat,reduction="multimodal_umap",group.by="seurat_clusters",label = TRUE,cols=alpha(col_vector,alpha_val))+ggtitle("Multimodal UMAP Clusters")
#Finally Plot results
plt<-(p1 | p2)/(p3 | p4)
ggsave(plt,file="Cluster1_multimodal.umap.pdf")
#system("slack -F phase1_multimodal.umap.pdf ryan_todo")
saveRDS(dat,file="Cluster1.rds")




