# Many different methods for Annotation has been used

# Input file
dat<-readRDS("Cluster1.rds")
DefaultAssay(dat) <- "RNA"

# Automated SingleR method
library(SingleR)
ref<-celldex::MouseRNAseqData()
results<-SingleR(test=as.SingleCellExperiment(dat),ref=ref,labels=ref$label.main)
dat$singlr_labels<-results$labels
ref<-celldex::HumanPrimaryCellAtlasData()
results<-SingleR(test=as.SingleCellExperiment(dat),ref=ref,labels=ref$label.fine)
big_table<-dat %>% as_tibble()
big_table %>%  filter (seurat_clusters==0) %>% select(singlr_labels) %>% group_by(singlr_labels) %>% summarize(N=n())
big_table %>%  group_by(seurat_clusters, singlr_labels) %>% summarize(N=n()

#Add colors
alpha_val=0.33
n <- 1000
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#pie(rep(1,n), col=sample(col_vector, n))
p4<-DimPlot(dat,reduction="multimodal_umap",group.by="singlr_labels",label=TRUE,cols=alpha(col_vector,alpha_val))+ggtitle("Multimodal UMAP Clusters")



# Second method but not very promising 
library(scRNAseq)
library(ExperimentHub)
eh<-ExperimentHub()
query(hub, "TENxBrainData")
query(eh,"TabulaMurisData")
eh[['EH1617']]
lung_ref<-eh[['EH1617']]
query(eh,"TENxBrainData")
lung_ref<-eh[['EH1042']]  # brain 

lung_ref<-lung_ref[,lung_ref$tissue=='Lung']
lung_ref<-lung_ref[,!is.na(lung_ref$cell_ontology_class)]
library(scuttle)
lung_ref<-logNormCounts(lung_ref)
results<-SingleR(test=as.SingleCellExperiment(dat),ref=lung_ref,labels=lung_ref$cell_ontology_class)
dat$singlr_labels<-results$labels
p4<-DimPlot(dat,reduction="multimodal_umap",group.by="singlr_labels",label=TRUE,cols=alpha(col_vector,alpha_val))+ggtitle("Multimodal UMAP Clusters")


# Third method using celldex
library(SingleR)
library(celldex)
ref <- celldex::MouseRNAseqData()
class(ref)
table(ref$label.main)
seu_int_SingleR <- SingleR::SingleR(test = Seurat::GetAssayData(dat, slot = "data"),
                                ref = ref,
                                labels = ref$label.main)
singleR_labels <- seu_int_SingleR$labels
head(seu_int_SingleR)

SingleR::plotScoreHeatmap(seu_int_SingleR)
SingleR::plotDeltaDistribution(seu_int_SingleR)
SingleR::plotDeltaDistribution(seu_int_SingleR)
singleR_labels <- seu_int_SingleR$labels
t <- table(singleR_labels)
other <- names(t)[t < 10]
singleR_labels[singleR_labels %in% other] <- "none"
dat$SingleR_annot <- singleR_labels
dittoSeq::dittoDimPlot(dat, "SingleR_annot", size = 0.7)



# another method to make annotation 
library(scCATCH)
obj <- createscCATCH(data = dat[['RNA']]@data, cluster = as.character(Idents(dat)))
obj <- findmarkergene(object = obj, species = "Human", marker = cellmatch, tissue = "Brain")
obj <- findcelltype(object = obj)
write.csv(obj@celltype[1:6],"cell_marker.csv")
obj <- findmarkergene(object = obj,
                      species = "Human",
                      marker = cellmatch,
                      tissue = c("Brain","Dorsolateral prefrontal cortex","Embryonic brain","Embryonic prefrontal cortex","Fetal brain","Hippocampus","Inferior colliculus","Midbrain","Sympathetic ganglion"))
obj <- findmarkergene(object = obj, species = "Mouse", marker = cellmatch,tissue = "Brain") 
obj <- findmarkergene(object = obj, species = "Human", marker = cellmatch, tissue = "Brain",use_method = "1")
obj <- findmarkergene(object = obj, species = "Human", marker = cellmatch, tissue = "Brain",use_method = "2")



