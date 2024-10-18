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
library(DoubletFinder)

library(BSgenome.Ggallus.UCSC.galGal4)
gtf <- rtracklayer::import('/proj/sllstore2017078/private/vikas/single_cell_pilot/annotation/Gallus_gallus.GRCg6a.104_filtered.gtf', format="gtf")
gene.coords <- gtf[gtf$type == 'gene' ]

my_list <- list()

#sink("./output.txt",append=T)
setupseurat<-function(i){
  #CHANGE
  setwd(paste0("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/",i,"/outs"))
  #sink("./output.txt",append=T)
  counts <- Read10X_h5("filtered_feature_bc_matrix.h5") #count data
  fragpath <- "atac_fragments.tsv.gz" #atac fragments
  metadata_cellranger<-read.csv("per_barcode_metrics.csv") #metadata
  row.names(metadata_cellranger)<-metadata_cellranger$barcode

  # create a Seurat object containing the RNA adata
  dat <- CreateSeuratObject(
    counts = counts$`Gene Expression`,
    assay = "RNA"
  )

  # create ATAC assay and add it to the object
  dat[["ATAC"]] <- CreateChromatinAssay(
    counts = counts$Peaks,
    sep = c(":", "-"),
    fragments = fragpath,
    annotation = gene.coords
  )

  DefaultAssay(dat) <- "ATAC"

  dat <- NucleosomeSignal(dat)
  dat <- TSSEnrichment(dat)
  dat<-AddMetaData(dat,metadata=metadata_cellranger)

  #Important Mitochondrial contamination
  #grep("^MT-|^ND1|^ND3|^ND4|^ND4L|^ND5|^ND6|^ATP6$|^ATP8$|^CYTB|^COII|^COX3",rownames(dat@assays$RNA@counts),value = TRUE)
  PercentageFeatureSet(dat,pattern="^MT-|^ND1|^ND3|^ND4|^ND4L|^ND5|^ND6|^ATP6$|^ATP8$|^CYTB|^COII|^COX3",assay = "RNA") -> dat$percent.mt
  
  # Change active assay
  DefaultAssay(dat) <- "RNA" 
  
  # Get the output of the file
  dat
  counts_val = dat[['RNA']]@counts
  ncells_Total<-ncol(counts_val)
  #write.table(ncells,file=Cell_Count.txt,sep = ",",append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE)
  # Need to add the plot script

    subset_data<- subset(dat,
    subset = nFeature_RNA > 500 &
    nFeature_RNA < 10000 &
    percent.mt < 20
  )

  # Get the output of subset data
  subset_data 
  counts_val = subset_data[['RNA']]@counts
  ncells_filt1<- ncol(counts_val)
  #write.table(ncells,file=Cell_Count.txt,sep = ",",append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE)




  # Doublet finding and removal
   subset_data = NormalizeData(subset_data,assay = "RNA")
   subset_data = FindVariableFeatures(subset_data, verbose = F)
   subset_data = ScaleData(subset_data, vars.to.regress = c("nFeature_RNA", "percent.mt"), verbose = F)
   subset_data = RunPCA(subset_data, verbose = F, npcs = 20)
   subset_data = FindNeighbors(subset_data, dims = 1:20)
   subset_data = FindClusters(subset_data)
   subset_data = RunUMAP(subset_data, dims = 1:10, verbose = F)
   ## pK Identification  ---------------------------------------------------------------------------------------
   sweep.res.list_pbmc <- paramSweep_v3(subset_data, PCs = 1:20, sct = FALSE)
   sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
   bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
   #ggplot(bcmvn_pbmc, aes(pK, BCmetric, group = 1)) + geom_point() + geom_line()
   
   # select the pK that corresponds to max bcmvn to optimize doublet detection
   pK <- bcmvn_pbmc %>% filter(BCmetric == max(BCmetric)) %>% select(pK) 
   pK <- as.numeric(as.character(pK[[1]]))


   # Change the percentage
   annotations <- subset_data@meta.data$seurat_clusters
   homotypic.prop <- modelHomotypic(annotations)
   nExp_poi <- round(0.10*nrow(subset_data))   # add expected percentage
   nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
   subset_data <- doubletFinder_v3(subset_data, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
   

   # Plot the results
   DF.name = colnames(subset_data@meta.data)[grepl("DF.classification", colnames(subset_data@meta.data))]
   #cowplot::plot_grid(ncol = 2, DimPlot(subset_data, group.by = "orig.ident") + NoAxes(), DimPlot(subset_data, group.by = DF.name) + NoAxes())
   #VlnPlot(subset_data, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)
   # Another way to plot
   #DimPlot(subset_data, reduction = 'umap', group.by = "DF.classifications_0.25_0.09_1182") 
 
   #table(subset_data@meta.data$DF.classifications_0.25_0.09_1182)
   subset_data = subset_data[, subset_data@meta.data[, DF.name] == "Singlet"]
   
   #Print the object again
   subset_data
   counts_val = subset_data[['RNA']]@counts
   ncells_filt2<-ncol(counts_val)
   
   All_mt<-c(ncells_Total,ncells_filt1,ncells_filt2)  
   write.csv(All_mt, file = "Final_Count.csv")
   my_list[[i]] <- All_mt 
   #write.table(ncells,file=Cell_Count.txt,sep = ",",append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE) 
   #sink()
   saveRDS(subset_data,file=paste0(i,".QC.1.rds"))
}




#generate all seurat objects
# Try with few 
#lapply(c(paste0("id_",c(32,37,40,44,46))),setupseurat)
#write.csv(my_list, file = "Final_Count.csv")
#sink()
lapply(c(paste0("id_",c(5,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46))),setupseurat)


