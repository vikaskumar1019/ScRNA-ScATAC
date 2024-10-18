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
gtf <- rtracklayer::import('/proj/sllstore2017078/private/vikas/single_cell_pilot/annotation/Gallus_gallus.GRCg6a.104_filtered.gtf', format="gtf")
gene.coords <- gtf[gtf$type == 'gene' ]



# CHANGE combine seurat objects and add metadata column for sample
#samp_1<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_1/outs/id_1.QC.1.rds"); samp_1$sample<-"sample_1"
#samp_2<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_2/outs/id_2.QC.1.rds"); samp_2$sample<-"sample_2"
#samp_3<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_3/outs/id_3.QC.1.rds"); samp_3$sample<-"sample_3"
#samp_4<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_4/outs/id_4.QC.1.rds"); samp_4$sample<-"sample_4"
samp_5<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_5/outs/id_5.QC.1.rds"); samp_5$sample<-"sample_5"
#samp_6<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_6/outs/id_6.QC.1.rds"); samp_6$sample<-"sample_6"
samp_7<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_7/outs/id_7.QC.1.rds"); samp_7$sample<-"sample_7"
#samp_8<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_8/outs/id_8.QC.1.rds"); samp_8$sample<-"sample_8"
#samp_9<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_9/outs/id_9.QC.1.rds"); samp_9$sample<-"sample_9"
samp_10<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_10/outs/id_10.QC.1.rds"); samp_10$sample<-"sample_10"
#samp_11<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_11/outs/id_11.QC.1.rds"); samp_11$sample<-"sample_11"
samp_12<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_12/outs/id_12.QC.1.rds"); samp_12$sample<-"sample_12"
samp_13<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_13/outs/id_13.QC.1.rds"); samp_13$sample<-"sample_13"
samp_14<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_14/outs/id_14.QC.1.rds"); samp_14$sample<-"sample_14"
samp_15<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_15/outs/id_15.QC.1.rds"); samp_15$sample<-"sample_15"
samp_16<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_16/outs/id_16.QC.1.rds"); samp_16$sample<-"sample_16"
samp_17<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_17/outs/id_17.QC.1.rds"); samp_17$sample<-"sample_17"
#samp_18<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_18/outs/id_18.QC.1.rds"); samp_18$sample<-"sample_18"
samp_19<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_19/outs/id_19.QC.1.rds"); samp_19$sample<-"sample_19"
samp_20<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_20/outs/id_20.QC.1.rds"); samp_20$sample<-"sample_20"
samp_21<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_21/outs/id_21.QC.1.rds"); samp_21$sample<-"sample_21"
samp_22<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_22/outs/id_22.QC.1.rds"); samp_22$sample<-"sample_22"
samp_23<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_23/outs/id_23.QC.1.rds"); samp_23$sample<-"sample_23"
samp_24<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_24/outs/id_24.QC.1.rds"); samp_24$sample<-"sample_24"
samp_25<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_25/outs/id_25.QC.1.rds"); samp_25$sample<-"sample_25"
#samp_26<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_26/outs/id_26.QC.1.rds"); samp_26$sample<-"sample_26"
#samp_27<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_27/outs/id_27.QC.1.rds"); samp_27$sample<-"sample_27"
samp_28<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_28/outs/id_28.QC.1.rds"); samp_28$sample<-"sample_28"
samp_29<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_29/outs/id_29.QC.1.rds"); samp_29$sample<-"sample_29"
samp_30<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_30/outs/id_30.QC.1.rds"); samp_30$sample<-"sample_30"
samp_31<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_31/outs/id_31.QC.1.rds"); samp_31$sample<-"sample_31"
samp_32<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_32/outs/id_32.QC.1.rds"); samp_32$sample<-"sample_32"
samp_33<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_33/outs/id_33.QC.1.rds"); samp_33$sample<-"sample_33"
samp_34<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_34/outs/id_34.QC.1.rds"); samp_34$sample<-"sample_34"
samp_35<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_35/outs/id_35.QC.1.rds"); samp_35$sample<-"sample_35"
samp_36<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_36/outs/id_36.QC.1.rds"); samp_36$sample<-"sample_36"
samp_37<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_37/outs/id_37.QC.1.rds"); samp_37$sample<-"sample_37"
samp_38<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_38/outs/id_38.QC.1.rds"); samp_38$sample<-"sample_38"
samp_39<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_39/outs/id_39.QC.1.rds"); samp_39$sample<-"sample_39"
samp_40<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_40/outs/id_40.QC.1.rds"); samp_40$sample<-"sample_40"
samp_41<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_41/outs/id_41.QC.1.rds"); samp_41$sample<-"sample_41"
samp_42<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_42/outs/id_42.QC.1.rds"); samp_42$sample<-"sample_42"
#samp_43<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_43/outs/id_43.QC.1.rds"); samp_43$sample<-"sample_43"
samp_44<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_44/outs/id_44.QC.1.rds"); samp_44$sample<-"sample_44"
samp_45<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_45/outs/id_45.QC.1.rds"); samp_45$sample<-"sample_45"
samp_46<-readRDS("/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_46/outs/id_46.QC.1.rds"); samp_46$sample<-"sample_46"


#dat <- merge(samp_1, y = c(samp_2,samp_3,samp_4,samp_5,samp_6,samp_7,samp_8,samp_10,samp_11,samp_12,samp_13,samp_14,samp_15,samp_16,samp_17,samp_18,samp_19,samp_20,samp_21,samp_22,samp_23,samp_24,samp_25,samp_26,samp_27,samp_28,samp_29,samp_30,samp_31,samp_32,samp_33,samp_34,samp_35,samp_36,samp_37,samp_38,samp_39,samp_40,samp_41,samp_42,samp_43,samp_44,samp_45,samp_46), add.cell.ids = c("samp_1","samp_2","samp_3","samp_4","samp_5","samp_6","samp_7","samp_8","samp_10","samp_11","samp_12","samp_13","samp_14","samp_15","samp_16","samp_17","samp_18","samp_19","samp_20","samp_21","samp_22","samp_23","samp_24","samp_25","samp_26","samp_27","samp_28","samp_29","samp_30","samp_31","samp_32","samp_33","samp_34","samp_35","samp_36","samp_37","samp_38","samp_39","samp_40","samp_41","samp_42","samp_43","samp_44","samp_45","samp_46"), project = "primary_phase1")

dat <- merge(samp_5, y = c(samp_7,samp_10,samp_12,samp_13,samp_14,samp_15,samp_16,samp_17,samp_19,samp_20,samp_21,samp_22,samp_23,samp_24,samp_25,samp_28,samp_29,samp_30,samp_31,samp_32,samp_33,samp_34,samp_35,samp_36,samp_37,samp_38,samp_39,samp_40,samp_41,samp_42,samp_44,samp_45,samp_46), add.cell.ids = c("samp_5","samp_7","samp_10","samp_12","samp_13","samp_14","samp_15","samp_16","samp_17","samp_19","samp_20","samp_21","samp_22","samp_23","samp_24","samp_25","samp_28","samp_29","samp_30","samp_31","samp_32","samp_33","samp_34","samp_35","samp_36","samp_37","samp_38","samp_39","samp_40","samp_41","samp_42","samp_44","samp_45","samp_46"), project = "primary_phase1")

saveRDS(dat,file="/proj/sllstore2017078/private/vikas/single_cell_pilot/cellranger_counts/check/id_11/remove/one/two/Check_NormDup_Rem/Final_Run/Merge/final_Merged.QC.rds")
