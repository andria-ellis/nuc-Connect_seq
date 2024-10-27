sessionInfo()
.libPaths(.libPaths()[-1])
set.seed(33)

# This code block processes single-cell RNA sequencing of neurons. Neurons are re-aligned and clustered. We re-assign specific clusters to excitatory ("Glu") or inhibitory ("GABA") cell types, and plots expression of marker genes across clusters. These excitatory and inhibirtory groups and saved seperately for further exploration.

#Input files: cds_neurons_withchen.rds

#Output files: cds_glu_chen.rds , cds_gaba_chen.rds

# Figures created: Supplementary figure 2

#loading libraries
library(RColorBrewer)
library(dplyr)
library(plyr)
library(ggplot2)
library(monocle3)
library(reshape2)
library(scales)
library(Matrix)
library(stringr)
library(repr)
library(Seurat)

#Load neuron data
cds=readRDS("cds_neurons_withchen.rds")

#Supplementary Figure 2C

cds$log.umi=log(cds$UMI)
plot_cells(cds,color_cells_by="log.umi",label_cell_groups=FALSE)+
  theme(legend.position="None",
        text=element_text(size=20))+
  scale_x_continuous(limits=c(-10, 5))

#Supplementary Figure 2B

plot_cells(cds,color_cells_by="lab",label_cell_groups=FALSE)+
  scale_color_manual(values = c("steelblue3","navy"))+
  theme(legend.position="None",
        text=element_text(size=20))+
  scale_x_continuous(limits=c(-10, 5))

cds$cell_label="none"
cds[,cds$chen.cluster=="Glu1"]$cell_label="yes"
cds[,cds$chen.cluster=="Glu2"]$cell_label="yes"
cds[,cds$chen.cluster=="Glu3"]$cell_label="yes"
cds[,cds$chen.cluster=="Glu4"]$cell_label="yes"
cds[,cds$chen.cluster=="Glu5"]$cell_label="yes"
cds[,cds$chen.cluster=="Glu6"]$cell_label="yes"
cds[,cds$chen.cluster=="Glu7"]$cell_label="yes"
cds[,cds$chen.cluster=="Glu8"]$cell_label="yes"
cds[,cds$chen.cluster=="Glu9"]$cell_label="yes"
cds[,cds$chen.cluster=="Glu10"]$cell_label="yes"
cds[,cds$chen.cluster=="Glu11"]$cell_label="yes"
cds[,cds$chen.cluster=="Glu12"]$cell_label="yes"
cds[,cds$chen.cluster=="Glu13"]$cell_label="yes"
cds[,cds$chen.cluster=="Glu14"]$cell_label="yes"
cds[,cds$chen.cluster=="Glu15"]$cell_label="yes"
cds[,cds$chen.cluster=="GABA1"]$cell_label="yes"
cds[,cds$chen.cluster=="GABA2"]$cell_label="yes"
cds[,cds$chen.cluster=="GABA3"]$cell_label="yes"
cds[,cds$chen.cluster=="GABA4"]$cell_label="yes"
cds[,cds$chen.cluster=="GABA5"]$cell_label="yes"
cds[,cds$chen.cluster=="GABA6"]$cell_label="yes"
cds[,cds$chen.cluster=="GABA7"]$cell_label="yes"
cds[,cds$chen.cluster=="GABA8"]$cell_label="yes"
cds[,cds$chen.cluster=="GABA9"]$cell_label="yes"
cds[,cds$chen.cluster=="GABA10"]$cell_label="yes"
cds[,cds$chen.cluster=="GABA11"]$cell_label="yes"
cds[,cds$chen.cluster=="GABA12"]$cell_label="yes"
cds[,cds$chen.cluster=="GABA13"]$cell_label="yes"
cds[,cds$chen.cluster=="GABA14"]$cell_label="yes"
cds[,cds$chen.cluster=="GABA15"]$cell_label="yes"
cds[,cds$chen.cluster=="GABA16"]$cell_label="yes"
cds[,cds$chen.cluster=="GABA17"]$cell_label="yes"
cds[,cds$chen.cluster=="GABA18"]$cell_label="yes"

#Supplementary figure 2D
plot_cells(cds,color_cells_by="cell_label",label_cell_groups = FALSE)+
  scale_colour_manual(values = c("grey","orangered")) +
  theme(legend.position="None",
        text=element_text(size=20))+
  scale_x_continuous(limits=c(-10, 5))

#Supplementary fiure 2A
cds$log.prv=log(cds$PRV)
plot_cells(cds,color_cells_by="log.prv",label_cell_groups=FALSE)+
  theme(legend.position="None",
        text=element_text(size=20))+
  scale_x_continuous(limits=c(-10, 5))

cds$typed="TRUE"
cds[,cds$chen.cluster=="none"]$typed="FALSE"
cds[,cds$chen.cluster=="zothers"]$typed="FALSE"

plot_cells(cds,color_cells_by="typed")+
  scale_color_manual(values = c("steelblue3","navy"))+
  text=element_text(size=20)+
  scale_x_continuous(limits=c(-10, 5))
#remove cells that have less than 5000 UMI

cds5000=cds[,-which(cds$UMI<5000)]

#remove cells that have a log(PRV)>=7.5

cds5000$logprv=log(cds5000$PRV)

cdstemp=cds5000[,-which(cds5000$logprv>=7.5)]

#downsample my data for the clustering

cds_chen=cdstemp[,cdstemp$lab=="chen"]
cds_buck=cdstemp[,cdstemp$lab=="buck"]

buck.cells=data.frame(cds_buck$cell_sample)
sampled.buck=sample_n(buck.cells, dim(cds_chen)[2])
colnames(buck.cells)=c("cells")
colnames(sampled.buck)=c("cells")

ind_list=match(sampled.buck$cells,buck.cells$cells)
cds_buck=cds_buck[,ind_list]

expression.matrix=cbind(exprs(cds_chen),exprs(cds_buck))
cell.names=rbind(data.frame(pData(cds_chen)),data.frame(pData(cds_buck)))
feature.names=join(data.frame(fData(cds_chen)),data.frame(fData(cds_buck)),
                   by = "gene_id",type="left",match="all")

row.names(expression.matrix) <- feature.names$gene_id
row.names(feature.names) <- feature.names$gene_id
row.names(cell.names) <- cell.names$cell_sample
colnames(expression.matrix) <- cell.names$cell_sample

cds <- new_cell_data_set(expression.matrix,
                         cell_metadata = cell.names,
                         gene_metadata = feature.names)

#Align to the chen data using the "lab" variable as the alignment group 

coldata_df=colData(cds) %>% as.data.frame()

count_mat = assay(cds)

coldata_df$Cell = cds$cell_sample
cds_seurat = 
  CreateSeuratObject(counts = count_mat,
                     project = "buck",
                     assay = "RNA",
                     meta.data = coldata_df)
cds.list <- 
  SplitObject(cds_seurat, 
              split.by = "lab")


features <- SelectIntegrationFeatures(object.list = cds.list)

cds.list <- lapply(X = cds.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = T)
  x <- RunPCA(x, features = features, verbose = T, npcs = 30)
  
})

anchors <- FindIntegrationAnchors(object.list = cds.list, 
                                  reduction = "rpca", 
                                  dims = 1:30,
                                  verbose = T,
                                  k.anchor = 30)
cds.integrated <- IntegrateData(anchorset = anchors, 
                                dims = 1:30,
                                verbose = T)

cds.integrated <- ScaleData(cds.integrated, verbose = T)

cds.integrated <- RunPCA(cds.integrated, verbose = T, npcs = 30)

cds.integrated <- RunUMAP(cds.integrated, dims = 1:30)

cds.integrated <- FindNeighbors(cds.integrated, reduction = "pca", dims = 1:30)

cds.integrated <- FindClusters(cds.integrated, resolution = 0.5)

# Visualization
p1 <- DimPlot(cds.integrated, reduction = "umap", split.by = "sample", group.by="sample", ncol=4)
p2 <- DimPlot(cds.integrated, reduction = "umap", label = TRUE, repel = TRUE)

# Pull out the UMAP embeddings from the Seurat CDS object
seurat_umap_embeddings = 
  as.data.frame(Embeddings(cds.integrated, reduction = "umap"))

seurat_umap_embeddings$Cell = rownames(seurat_umap_embeddings)

cdstemp$umap1=seurat_umap_embeddings$UMAP_1
cdstemp$umap2=seurat_umap_embeddings$UMAP_2

df=as.data.frame(colData(cdstemp))
saved_coords =as.matrix(data.frame(df$umap1,df$umap2))


reducedDims(cdstemp)[["UMAP"]] = saved_coords

cdsfinal1 <- cluster_cells(cdstemp,resolution=1e-1)
colData(cdsfinal1)$cluster<- as.character(clusters(cdsfinal1))

g1=plot_cells(cdsfinal1)
ggsave(g1,file="resolution_high.png")

cdsfinal2 <- cluster_cells(cdstemp,resolution=1e-2)

g2=plot_cells(cdsfinal2)
ggsave(g2,file="resolution_med.png")

cdsfinal3 <- cluster_cells(cdstemp,resolution=1e-3)


nm_list=c("Gad1","Gad2","Slc17a6","Slc17a7","Slc17a8","Slc32a1","Fos")
p2=plot_genes_by_group(cdsfinal2,
                       nm_list,
                       group_cells_by="cluster",
                       ordering_type="maximal_on_diag",
                       max.size=3)
ggsave(p2,height=4,width=20,file="nm_markers.png")

#Marker levels allow us to label by clusters manually

colData(cdsfinal2)$assigned_cell_type = dplyr::recode(colData(cdsfinal2)$cluster,
                                                      "1"="Glu",
                                                      "2"="Glu",
                                                      "3"="Gaba",
                                                      "4"="Gaba",
                                                      "5"="Glu",
                                                      "6"="Gaba",
                                                      "7"="Glu",
                                                      "8"="Gaba",
                                                      "9"="Glu",
                                                      "10"="Gaba",
                                                      "11"="Gaba",
                                                      "12"="Glu",
                                                      "13"="Glu",
                                                      "14"="Gaba",
                                                      "15"="Gaba",
                                                      "16"="Gaba",
                                                      "17"="Gaba",
                                                      "18"="Glu",
                                                      "19"="Glu",
                                                      "20"="Gaba",
                                                      "21"="Gaba",
                                                      "22"="Glu",
                                                      "23"="Gaba",
                                                      "24"="Gaba",
                                                      "25"="Gaba",
                                                      "26"="Gaba",
                                                      "27"="Gaba",
                                                      "28"="Gaba",
                                                      "29"="Glu",
                                                      "30"="Gaba",
                                                      "31"="Gaba",
                                                      "32"="Glu",
                                                      "33"="Glu",
                                                      "34"="Gaba",
                                                      "35"="Gaba",
                                                      "36"="Gaba",
                                                      "37"="Gaba",
                                                      "38"="Glu",
                                                      "39"="Gaba",
                                                      "40"="Gaba",
                                                      "41"="Glu",
                                                      "42"="Gaba",
                                                      "43"="Glu",
                                                      "44"="Gaba",
                                                      "45"="Glu",
                                                      "46"="Gaba",
                                                      "47"="Gaba",
                                                      "48"="Gaba",
                                                      "49"="Gaba",
                                                      "50"="Glu",
                                                      "51"="Glu",
                                                      "52"="Gaba",
                                                      "53"="Gaba",
                                                      "54"="Gaba",
                                                      "55"="Glu",
                                                      "56"="Glu",
                                                      "57"="Glu",
                                                      "58"="Glu",
                                                      "59"="Glu",
                                                      "60"="Glu",
                                                      "61"="Gaba",
                                                      "62"="Gaba",
                                                      "63"="None",
                                                      "64"="Gaba",
                                                      "65"="None",
                                                      "66"="Gaba",
                                                      "67"="Glu",
                                                      "68"="Glu",
                                                      "69"="Gaba",
                                                      "70"="Fos")

nm_cell_type=plot_cells(cdsfinal2,color_cells_by="assigned_cell_type")

plot_cells(cdsfinal2,color_cells_by="assigned_cell_type")+
  text=element_text(size=20)+
  scale_x_continuous(limits=c(-10, 5))

cdsfinal2$temp=TRUE
plot_cells(cdsfinal2,color_cells_by="temp",label_cell_groups=FALSE)+
  scale_color_manual(values = c("grey"))

plot_cells(cdsfinal2,gene="Slc16a6",label_cell_groups=FALSE,cell_size=1)

ggsave(nm_cell_type,file="cell_type_nm.png")

cds_gaba=cdsfinal2[,cdsfinal2$assigned_cell_type=="Gaba"]
cds_glu=cdsfinal2[,cdsfinal2$assigned_cell_type=="Glu"]

saveRDS(cdsfinal1,file="cdsfinal_res1.rds")
saveRDS(cdsfinal2,file="cdsfinal_res2.rds")
saveRDS(cdsfinal3,file="cdsfinal_res3.rds")

#############################finding gaba and glu cells######################################

cds=cluster_cells(cdsfinal1)
colData(cds)$cluster<- as.character(clusters(cds))

cds1=cds[,cds$cluster=="1"]
cds2=cds[,cds$cluster=="2"]
cds3=cds[,cds$cluster=="3"]


###########seurat for cluster 1 ##############
coldata_df=colData(cds1) %>% as.data.frame()

count_mat = assay(cds1)

coldata_df$Cell = cds1$cell_sample
cds_seurat = 
  CreateSeuratObject(counts = count_mat,
                     project = "buck",
                     assay = "RNA",
                     meta.data = coldata_df)
cds.list <- 
  SplitObject(cds_seurat, 
              split.by = "lab")


features <- SelectIntegrationFeatures(object.list = cds.list)

cds.list <- lapply(X = cds.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = T)
  x <- RunPCA(x, features = features, verbose = T, npcs = 30)
  
})

anchors <- FindIntegrationAnchors(object.list = cds.list, 
                                  reduction = "rpca", 
                                  dims = 1:30,
                                  verbose = T,
                                  k.anchor = 30)
cds.integrated <- IntegrateData(anchorset = anchors, 
                                dims = 1:30,
                                verbose = T)

cds.integrated <- ScaleData(cds.integrated, verbose = T)

cds.integrated <- RunPCA(cds.integrated, verbose = T, npcs = 30)

cds.integrated <- RunUMAP(cds.integrated, dims = 1:30)

cds.integrated <- FindNeighbors(cds.integrated, reduction = "pca", dims = 1:30)

cds.integrated <- FindClusters(cds.integrated, resolution = 0.5)

# Visualization
p1 <- DimPlot(cds.integrated, reduction = "umap", split.by = "sample", group.by="sample", ncol=4)
p2 <- DimPlot(cds.integrated, reduction = "umap", label = TRUE, repel = TRUE)

# Pull out the UMAP embeddings from the Seurat CDS object
seurat_umap_embeddings = 
  as.data.frame(Embeddings(cds.integrated, reduction = "umap"))

seurat_umap_embeddings$Cell = rownames(seurat_umap_embeddings)

cds1$umap1=seurat_umap_embeddings$UMAP_1
cds1$umap2=seurat_umap_embeddings$UMAP_2

df=as.data.frame(colData(cds1))
saved_coords =as.matrix(data.frame(df$umap1,df$umap2))


reducedDims(cds1)[["UMAP"]] = saved_coords

cds1 <- cluster_cells(cds1,resolution=1e-2)
saveRDS(cds1,file="cluster1.rds")

###########seurat for cluster 2 ##############
coldata_df=colData(cds2) %>% as.data.frame()

count_mat = assay(cds2)

coldata_df$Cell = cds2$cell_sample
cds_seurat = 
  CreateSeuratObject(counts = count_mat,
                     project = "buck",
                     assay = "RNA",
                     meta.data = coldata_df)
cds.list <- 
  SplitObject(cds_seurat, 
              split.by = "lab")


features <- SelectIntegrationFeatures(object.list = cds.list)

cds.list <- lapply(X = cds.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = T)
  x <- RunPCA(x, features = features, verbose = T, npcs = 30)
  
})

anchors <- FindIntegrationAnchors(object.list = cds.list, 
                                  reduction = "rpca", 
                                  dims = 1:30,
                                  verbose = T,
                                  k.anchor = 30)
cds.integrated <- IntegrateData(anchorset = anchors, 
                                dims = 1:30,
                                verbose = T)

cds.integrated <- ScaleData(cds.integrated, verbose = T)

cds.integrated <- RunPCA(cds.integrated, verbose = T, npcs = 30)

cds.integrated <- RunUMAP(cds.integrated, dims = 1:30)

cds.integrated <- FindNeighbors(cds.integrated, reduction = "pca", dims = 1:30)

cds.integrated <- FindClusters(cds.integrated, resolution = 0.5)

# Visualization
p1 <- DimPlot(cds.integrated, reduction = "umap", split.by = "sample", group.by="sample", ncol=4)
p2 <- DimPlot(cds.integrated, reduction = "umap", label = TRUE, repel = TRUE)

# Pull out the UMAP embeddings from the Seurat CDS object
seurat_umap_embeddings = 
  as.data.frame(Embeddings(cds.integrated, reduction = "umap"))

seurat_umap_embeddings$Cell = rownames(seurat_umap_embeddings)

cds2$umap1=seurat_umap_embeddings$UMAP_1
cds2$umap2=seurat_umap_embeddings$UMAP_2

df=as.data.frame(colData(cds2))
saved_coords =as.matrix(data.frame(df$umap1,df$umap2))


reducedDims(cds2)[["UMAP"]] = saved_coords

cds2 <- cluster_cells(cds2,resolution=1e-2)
saveRDS(cds2,file="cluster2.rds")

###########seurat for cluster 3 ##############
coldata_df=colData(cds3) %>% as.data.frame()

count_mat = assay(cds3)

coldata_df$Cell = cds3$cell_sample
cds_seurat = 
  CreateSeuratObject(counts = count_mat,
                     project = "buck",
                     assay = "RNA",
                     meta.data = coldata_df)
cds.list <- 
  SplitObject(cds_seurat, 
              split.by = "lab")


features <- SelectIntegrationFeatures(object.list = cds.list)

cds.list <- lapply(X = cds.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = T)
  x <- RunPCA(x, features = features, verbose = T, npcs = 30)
  
})

anchors <- FindIntegrationAnchors(object.list = cds.list, 
                                  reduction = "rpca", 
                                  dims = 1:30,
                                  verbose = T,
                                  k.anchor = 30)
cds.integrated <- IntegrateData(anchorset = anchors, 
                                dims = 1:30,
                                verbose = T)

cds.integrated <- ScaleData(cds.integrated, verbose = T)

cds.integrated <- RunPCA(cds.integrated, verbose = T, npcs = 30)

cds.integrated <- RunUMAP(cds.integrated, dims = 1:30)

cds.integrated <- FindNeighbors(cds.integrated, reduction = "pca", dims = 1:30)

cds.integrated <- FindClusters(cds.integrated, resolution = 0.5)

# Visualization
p1 <- DimPlot(cds.integrated, reduction = "umap", split.by = "sample", group.by="sample", ncol=4)
p2 <- DimPlot(cds.integrated, reduction = "umap", label = TRUE, repel = TRUE)

# Pull out the UMAP embeddings from the Seurat CDS object
seurat_umap_embeddings = 
  as.data.frame(Embeddings(cds.integrated, reduction = "umap"))

seurat_umap_embeddings$Cell = rownames(seurat_umap_embeddings)

cds3$umap1=seurat_umap_embeddings$UMAP_1
cds3$umap2=seurat_umap_embeddings$UMAP_2

df=as.data.frame(colData(cds3))
saved_coords =as.matrix(data.frame(df$umap1,df$umap2))


reducedDims(cds3)[["UMAP"]] = saved_coords

cds3 <- cluster_cells(cds3,resolution=1e-2)
saveRDS(cds3,file="cluster3.rds")


##################################seurat for gaba cells######################################

coldata_df=colData(cds_gaba) %>% as.data.frame()

count_mat = assay(cds_gaba)

coldata_df$Cell = cds_gaba$cell_sample
cds_seurat = 
  CreateSeuratObject(counts = count_mat,
                     project = "buck",
                     assay = "RNA",
                     meta.data = coldata_df)
cds.list <- 
  SplitObject(cds_seurat, 
              split.by = "lab")


features <- SelectIntegrationFeatures(object.list = cds.list)

cds.list <- lapply(X = cds.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = T)
  x <- RunPCA(x, features = features, verbose = T, npcs = 20)
  
})

anchors <- FindIntegrationAnchors(object.list = cds.list, 
                                  reduction = "rpca", 
                                  dims = 1:20,
                                  verbose = T,
                                  k.anchor = 20)
cds.integrated <- IntegrateData(anchorset = anchors, 
                                dims = 1:20,
                                verbose = T)

cds.integrated <- ScaleData(cds.integrated, verbose = T)

cds.integrated <- RunPCA(cds.integrated, verbose = T, npcs = 20)

cds.integrated <- RunUMAP(cds.integrated, dims = 1:20)

cds.integrated <- FindNeighbors(cds.integrated, reduction = "pca", dims = 1:20)

cds.integrated <- FindClusters(cds.integrated, resolution = 0.5)

# Visualization
p1 <- DimPlot(cds.integrated, reduction = "umap", split.by = "sample", group.by="sample", ncol=4)
p2 <- DimPlot(cds.integrated, reduction = "umap", label = TRUE, repel = TRUE)

# Pull out the UMAP embeddings from the Seurat CDS object
seurat_umap_embeddings = 
  as.data.frame(Embeddings(cds.integrated, reduction = "umap"))

seurat_umap_embeddings$Cell = rownames(seurat_umap_embeddings)

cds_gaba$umap1=seurat_umap_embeddings$UMAP_1
cds_gaba$umap2=seurat_umap_embeddings$UMAP_2

df=as.data.frame(colData(cds_gaba))
saved_coords =as.matrix(data.frame(df$umap1,df$umap2))


reducedDims(cdsgaba)[["UMAP"]] = saved_coords

cds_gaba <- cluster_cells(cdsgaba,resolution=1e-2)
saveRDS(cds_gaba,file="cds_gaba_chen.rds")

################################seurat for glu cells######################################

coldata_df=colData(cds_glu) %>% as.data.frame()

count_mat = assay(cds_glu)

coldata_df$Cell = cds_glu$cell_sample
cds_seurat = 
  CreateSeuratObject(counts = count_mat,
                     project = "buck",
                     assay = "RNA",
                     meta.data = coldata_df)
cds.list <- 
  SplitObject(cds_seurat, 
              split.by = "lab")


features <- SelectIntegrationFeatures(object.list = cds.list)

cds.list <- lapply(X = cds.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = T)
  x <- RunPCA(x, features = features, verbose = T, npcs = 20)
  
})

anchors <- FindIntegrationAnchors(object.list = cds.list, 
                                  reduction = "rpca", 
                                  dims = 1:20,
                                  verbose = T,
                                  k.anchor = 20)
cds.integrated <- IntegrateData(anchorset = anchors, 
                                dims = 1:20,
                                verbose = T)

cds.integrated <- ScaleData(cds.integrated, verbose = T)

cds.integrated <- RunPCA(cds.integrated, verbose = T, npcs = 20)

cds.integrated <- RunUMAP(cds.integrated, dims = 1:20)

cds.integrated <- FindNeighbors(cds.integrated, reduction = "pca", dims = 1:20)

cds.integrated <- FindClusters(cds.integrated, resolution = 0.5)

# Visualization
p1 <- DimPlot(cds.integrated, reduction = "umap", split.by = "sample", group.by="sample", ncol=4)
p2 <- DimPlot(cds.integrated, reduction = "umap", label = TRUE, repel = TRUE)

# Pull out the UMAP embeddings from the Seurat CDS object
seurat_umap_embeddings = 
  as.data.frame(Embeddings(cds.integrated, reduction = "umap"))

seurat_umap_embeddings$Cell = rownames(seurat_umap_embeddings)

cds_glu$umap1=seurat_umap_embeddings$UMAP_1
cds_glu$umap2=seurat_umap_embeddings$UMAP_2

df=as.data.frame(colData(cds_glu))
saved_coords =as.matrix(data.frame(df$umap1,df$umap2))


reducedDims(cds_glu)[["UMAP"]] = saved_coords

cds_glu <- cluster_cells(cds_glu,resolution=1e-2)
saveRDS(cds_glu,file="cds_glu_chen.rds")