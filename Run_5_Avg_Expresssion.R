library(readr)

# Get the Average expression
dat$Samp_Cluster<-paste (dat$sample, dat$seurat_clusters,sep = "_")
avg_expr <- AverageExpression(dat, assays="RNA", group.by="Samp_Cluster")
write.csv(avg_expr, file="new_check.csv")
Clust0<-read_csv("new_check.csv")

#Clust0_renamed<-rename(Clust0, Gene_0 = 1)
#vk<-Clust0_renamed %>% select(ends_with('_0'))
#write.csv(vk, file="Final_Clust0.csv",row.names=FALSE)
#names(vk)<-c('id','10','12','13','14','15','16','17','19','20','21','22','23','24','25','28','29','30','31','32','33','34','35','36','37','38','39','40','41','42','44','45','46','5','7')
#first_col<-Clust0[,1]

library("tidyseurat")

Clust0<-read_csv("new_check.csv")
#Clust0_renamed<-dplyr::rename(Clust0, Gene_0 = 1)
vk<-Clust0_renamed %>% select(ends_with('_0'))
write.csv(vk, file="Final_Clust0.csv",row.names=FALSE)

Clust0_renamed<-dplyr::rename(Clust0, Gene_1 = 1)
vk<-Clust0_renamed %>% select(ends_with('_1'))
write.csv(vk, file="Final_Clust1.csv",row.names=FALSE)

Clust0_renamed<-dplyr::rename(Clust0, Gene_2 = 1)
vk<-Clust0_renamed %>% select(ends_with('_2'))
write.csv(vk, file="Final_Clust2.csv",row.names=FALSE)

Clust0_renamed<-dplyr::rename(Clust0, Gene_3 = 1)
vk<-Clust0_renamed %>% select(ends_with('_3'))
write.csv(vk, file="Final_Clust3.csv",row.names=FALSE)

Clust0_renamed<-dplyr::rename(Clust0, Gene_4 = 1)
vk<-Clust0_renamed %>% select(ends_with('_4'))
write.csv(vk, file="Final_Clust4.csv",row.names=FALSE)

Clust0_renamed<-dplyr::rename(Clust0, Gene_5 = 1)
vk<-Clust0_renamed %>% select(ends_with('_5'))
write.csv(vk, file="Final_Clust5.csv",row.names=FALSE)

Clust0_renamed<-dplyr::rename(Clust0, Gene_6 = 1)
vk<-Clust0_renamed %>% select(ends_with('_6'))
write.csv(vk, file="Final_Clust6.csv",row.names=FALSE)

Clust0_renamed<-dplyr::rename(Clust0, Gene_7 = 1)
vk<-Clust0_renamed %>% select(ends_with('_7'))
write.csv(vk, file="Final_Clust7.csv",row.names=FALSE)

Clust0_renamed<-dplyr::rename(Clust0, Gene_8 = 1)
vk<-Clust0_renamed %>% select(ends_with('_8'))
write.csv(vk, file="Final_Clust8.csv",row.names=FALSE)

Clust0_renamed<-dplyr::rename(Clust0, Gene_9 = 1)
vk<-Clust0_renamed %>% select(ends_with('_9'))
write.csv(vk, file="Final_Clust9.csv",row.names=FALSE)

Clust0_renamed<-dplyr::rename(Clust0, Gene_10 = 1)
vk<-Clust0_renamed %>% select(ends_with('_10'))
write.csv(vk, file="Final_Clust10.csv",row.names=FALSE)

Clust0_renamed<-dplyr::rename(Clust0, Gene_11 = 1)
vk<-Clust0_renamed %>% select(ends_with('_11'))
write.csv(vk, file="Final_Clust11.csv",row.names=FALSE)

Clust0_renamed<-dplyr::rename(Clust0, Gene_12 = 1)
vk<-Clust0_renamed %>% select(ends_with('_12'))
write.csv(vk, file="Final_Clust12.csv",row.names=FALSE)

Clust0_renamed<-dplyr::rename(Clust0, Gene_13 = 1)
vk<-Clust0_renamed %>% select(ends_with('_13'))
write.csv(vk, file="Final_Clust13.csv",row.names=FALSE)

Clust0_renamed<-dplyr::rename(Clust0, Gene_14 = 1)
vk<-Clust0_renamed %>% select(ends_with('_14'))
write.csv(vk, file="Final_Clust14.csv",row.names=FALSE)

Clust0_renamed<-dplyr::rename(Clust0, Gene_15 = 1)
vk<-Clust0_renamed %>% select(ends_with('_15'))
write.csv(vk, file="Final_Clust15.csv",row.names=FALSE)

Clust0_renamed<-dplyr::rename(Clust0, Gene_16 = 1)
vk<-Clust0_renamed %>% select(ends_with('_16'))
write.csv(vk, file="Final_Clust16.csv",row.names=FALSE)

Clust0_renamed<-dplyr::rename(Clust0, Gene_17 = 1)
vk<-Clust0_renamed %>% select(ends_with('_17'))
write.csv(vk, file="Final_Clust17.csv",row.names=FALSE)

Clust0_renamed<-dplyr::rename(Clust0, Gene_18 = 1)
vk<-Clust0_renamed %>% select(ends_with('_18'))
write.csv(vk, file="Final_Clust18.csv",row.names=FALSE)

Clust0_renamed<-dplyr::rename(Clust0, Gene_19 = 1)
vk<-Clust0_renamed %>% select(ends_with('_19'))
write.csv(vk, file="Final_Clust19.csv",row.names=FALSE)

Clust0_renamed<-dplyr::rename(Clust0, Gene_20 = 1)
vk<-Clust0_renamed %>% select(ends_with('_20'))
write.csv(vk, file="Final_Clust20.csv",row.names=FALSE)

Clust0_renamed<-dplyr::rename(Clust0, Gene_21 = 1)
vk<-Clust0_renamed %>% select(ends_with('_21'))
write.csv(vk, file="Final_Clust21.csv",row.names=FALSE)

Clust0_renamed<-dplyr::rename(Clust0, Gene_22 = 1)
vk<-Clust0_renamed %>% select(ends_with('_22'))
write.csv(vk, file="Final_Clust22.csv",row.names=FALSE)

Clust0_renamed<-dplyr::rename(Clust0, Gene_23 = 1)
vk<-Clust0_renamed %>% select(ends_with('_23'))
write.csv(vk, file="Final_Clust23.csv",row.names=FALSE)

Clust0_renamed<-dplyr::rename(Clust0, Gene_24 = 1)
vk<-Clust0_renamed %>% select(ends_with('_24'))
write.csv(vk, file="Final_Clust24.csv",row.names=FALSE)

Clust0_renamed<-dplyr::rename(Clust0, Gene_25 = 1)
vk<-Clust0_renamed %>% select(ends_with('_25'))
write.csv(vk, file="Final_Clust25.csv",row.names=FALSE)

Clust0_renamed<-dplyr::rename(Clust0, Gene_26 = 1)
vk<-Clust0_renamed %>% select(ends_with('_26'))
write.csv(vk, file="Final_Clust26.csv",row.names=FALSE)

Clust0_renamed<-dplyr::rename(Clust0, Gene_27 = 1)
vk<-Clust0_renamed %>% select(ends_with('_27'))
write.csv(vk, file="Final_Clust27.csv",row.names=FALSE)

Clust0_renamed<-dplyr::rename(Clust0, Gene_28 = 1)
vk<-Clust0_renamed %>% select(ends_with('_28'))
write.csv(vk, file="Final_Clust28.csv",row.names=FALSE)

Clust0_renamed<-dplyr::rename(Clust0, Gene_29 = 1)
vk<-Clust0_renamed %>% select(ends_with('_29'))
write.csv(vk, file="Final_Clust29.csv",row.names=FALSE)

Clust0_renamed<-dplyr::rename(Clust0, Gene_30 = 1)
vk<-Clust0_renamed %>% select(ends_with('_30'))
write.csv(vk, file="Final_Clust30.csv",row.names=FALSE)

Clust0_renamed<-dplyr::rename(Clust0, Gene_31 = 1)
vk<-Clust0_renamed %>% select(ends_with('_31'))
write.csv(vk, file="Final_Clust31.csv",row.names=FALSE)

Clust0_renamed<-dplyr::rename(Clust0, Gene_32 = 1)
vk<-Clust0_renamed %>% select(ends_with('_32'))
write.csv(vk, file="Final_Clust32.csv",row.names=FALSE)

Clust0_renamed<-dplyr::rename(Clust0, Gene_33 = 1)
vk<-Clust0_renamed %>% select(ends_with('_33'))
write.csv(vk, file="Final_Clust33.csv",row.names=FALSE)

Clust0_renamed<-dplyr::rename(Clust0, Gene_34 = 1)
vk<-Clust0_renamed %>% select(ends_with('_34'))
write.csv(vk, file="Final_Clust34.csv",row.names=FALSE)

Clust0_renamed<-dplyr::rename(Clust0, Gene_35 = 1)
vk<-Clust0_renamed %>% select(ends_with('_35'))
write.csv(vk, file="Final_Clust35.csv",row.names=FALSE)

Clust0_renamed<-dplyr::rename(Clust0, Gene_36 = 1)
vk<-Clust0_renamed %>% select(ends_with('_36'))
write.csv(vk, file="Final_Clust36.csv",row.names=FALSE)
