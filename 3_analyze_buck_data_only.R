sessionInfo()
set.seed(33)

# The purpose of this code block is to plot clusters by cell type. Genes indicative of certain cell types, such as oligodendrocytes, astrocytes, microglia, and neurons, are visualized with distinct color gradients. The code also annotates clusters with assigned cell types.

# Input files: cds_filtered_nochen.rds

#Figures created:   Figure 1
#                   Figure 2
#                   Supplementary Figure 1


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

#Read in the filtered CDS with just our lab's data

cds=readRDS("Downloads/cds_filtered_nochen.rds")

fData(cds)=fData(cds)[,10:12]

pData(cds)$condition=pData(cds)$sample_name
pData(cds)=pData(cds)[,-3]

#Re-align the CDS to remove the cluster with the chen data

cds_buck <- preprocess_cds(cds,num_dim=45)

cds_buck <- align_cds(cds_buck, alignment_group = "sample")

cds_buck <- reduce_dimension(cds_buck)

cds_buck <- cluster_cells(cds_buck,resolution=1e-4)

#The following plots make up the marker gene plots for Figure 2, these help inform about which clusters and different cell types and how many neurons we collected.

#oligodendrycytes
p1=plot_cells(cds_buck,genes=c("Olig1"),label_groups_by_cluster=FALSE,label_cell_groups = FALSE,cell_size=1) +
  scale_colour_gradientn(colours = c("#CCCCCC","#CCCCCC","#000333")) +
  theme(legend.position="None",
        axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        text=element_text(size=50))

ggsave(p1,filename="buck.olig1.png")

p2=plot_cells(cds_buck,genes=c("Top2a"),label_groups_by_cluster=FALSE,label_cell_groups = FALSE,cell_size=1) +
  scale_colour_gradientn(colours = c("#CCCCCC","#CCCCCC","#000333")) +
  theme(legend.position="None",
        axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        text=element_text(size=50))

ggsave(p2,filename="buck.top2a.png")

p3=plot_cells(cds_buck,genes=c("Pdgfra"),label_groups_by_cluster=FALSE,label_cell_groups = FALSE,cell_size=1) +
  scale_colour_gradientn(colours = c("#CCCCCC","#CCCCCC","#000333")) +
  theme(legend.position="None",
        axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        text=element_text(size=50))

ggsave(p3,filename="buck.pdgfra.png")

#astrocytes
p4=plot_cells(cds_buck,genes=c("Aldh1l1"),label_groups_by_cluster=FALSE,label_cell_groups = FALSE,cell_size=1) +
  scale_colour_gradientn(colours = c("#CCCCCC","#CCCCCC","#000333")) +
  theme(legend.position="None",
        axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        text=element_text(size=50))

ggsave(p4,filename="buck.aldh1l1.png")

p5=plot_cells(cds_buck,genes=c("Gfap"),label_groups_by_cluster=FALSE,label_cell_groups = FALSE,cell_size=1) +
  scale_colour_gradientn(colours = c("#CCCCCC","#CCCCCC","#000333"))+
  theme(legend.position="None",
        axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        text=element_text(size=50))

ggsave(p5,filename="buck.gfap.png")

p6=plot_cells(cds_buck,genes=c("Aldoc"),label_groups_by_cluster=FALSE,label_cell_groups = FALSE,cell_size=1) +
  scale_colour_gradientn(colours = c("#CCCCCC","#CCCCCC","#000333"))+
  theme(legend.position="None",
        axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        text=element_text(size=50))

ggsave(p6,filename="buck.aldoc.png")

#microglia
p7=plot_cells(cds_buck,genes=c("Aif1"),label_groups_by_cluster=FALSE,label_cell_groups = FALSE,cell_size=1) +
  scale_colour_gradientn(colours = c("#CCCCCC","#CCCCCC","#000333")) +
  theme(legend.position="None",
        axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        text=element_text(size=50))

ggsave(p7,filename="buck.aif1.png")

p8=plot_cells(cds_buck,genes=c("Hexb"),label_groups_by_cluster=FALSE,label_cell_groups = FALSE,cell_size=1) +
  scale_colour_gradientn(colours = c("#CCCCCC","#CCCCCC","#000333")) +
  theme(legend.position="None",
        axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        text=element_text(size=50))

ggsave(p8,filename="buck.hexb.png")

p9=plot_cells(cds_buck,genes=c("Mafb"),label_groups_by_cluster=FALSE,label_cell_groups = FALSE,cell_size=1) +
  scale_colour_gradientn(colours = c("#CCCCCC","#CCCCCC","#000333")) +
  theme(legend.position="None",
        axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        text=element_text(size=50))

ggsave(p9,filename="buck.mafb.png")

#neurons

p10=plot_cells(cds_buck,genes=c("Syt1"),label_groups_by_cluster=FALSE,label_cell_groups = FALSE,cell_size=1) +
  scale_colour_gradientn(colours = c("#CCCCCC","#CCCCCC","#000333")) +
  theme(legend.position="None",
        axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        text=element_text(size=50))

ggsave(p10,filename="buck.syt1.png")

p11=plot_cells(cds_buck,genes=c("Map2"),label_groups_by_cluster=FALSE,label_cell_groups = FALSE,cell_size=1) +
  scale_colour_gradientn(colours = c("#CCCCCC","#CCCCCC","#000333")) +
  theme(legend.position="None",
        axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        text=element_text(size=50))

ggsave(p11,filename="buck.map2.png")

plot_cells(cds_buck,genes=c("Snap25"),label_groups_by_cluster=FALSE,label_cell_groups = FALSE,cell_size=1) +
  scale_colour_gradientn(colours = c("#CCCCCC","#CCCCCC","#000333")) +
  theme(legend.position="None",
        axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        text=element_text(size=50))


colData(cds_buck)$cluster<- as.character(clusters(cds_buck))

colData(cds_buck)$assigned_cell_type = dplyr::recode(colData(cds_buck)$cluster,
                                                     "1"="Neuron",
                                                     "2"="Microglia",
                                                     "3"="Neuron",
                                                     "4"="Oligodendrocyte",
                                                     "5"="Microglia",
                                                     "6"="Neuron",
                                                     "7"="Neuron",
                                                     "8"="Neuron",
                                                     "9"="Astrocyte",
                                                     "10"="Microglia",
                                                     "11"="Neuron",
                                                     "12"="O.P.C.",
                                                     "13"="Unknown",
                                                     "14"="P.O.P.C.",
                                                     "15"="Neuron",
                                                     "16"="Unknown",
                                                     "17"="Unknown",
                                                     "18"="Unknown",
                                                     "19"="Unknown",
                                                     "20"="Neuron",
                                                     "21"="Neuron")


# Plot Figure 1B with assigned cell types
plot_cells(cds_buck,color_cells_by="assigned_cell_type",cell_size=1,label_cell_groups = FALSE)+theme(legend.position="none")+
  scale_colour_manual(values = c("orange2","green4",
                                 "dodgerblue2","orangered","darkmagenta","orchid","midnightblue"))

# Plot miniature for Fig1A
plot_cells(cds_buck,color_cells_by="assigned_cell_type",cell_size=1,label_cell_groups = FALSE)+theme(legend.position="none")+
  scale_colour_manual(values = c("orange2","green4",
                                 "dodgerblue2","orangered","darkmagenta","orchid","midnightblue"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

#Load garnett mouse brain classifier to test our cell type assignments against.

library(garnett)

classifier <- readRDS("Downloads/mmBrain_20191017.RDS")

pbmc_cds <- classify_cells(cds_buck, classifier,
                           db = org.Hs.eg.db,
                           cluster_extend = TRUE)

mm_neuro_cds=pbmc_cds
head(pData(mm_neuro_cds))

df=data.frame(mm_neuro_cds$umap1,mm_neuro_cds$umap2,mm_neuro_cds$cluster_ext_type)
colnames(df)=c("umap1","umap2","type")


df=as.data.frame(colData(mm_neuro_cds))
saved_coords =as.matrix(data.frame(df$umap1,df$umap2))

reducedDims(mm_neuro_cds)[["UMAP"]] = saved_coords

cds <- cluster_cells(mm_neuro_cds)

cds$exp=cds$sample_name
cds$sample_name=NULL

plot_cells(cds)

# Supplementary Figure 1
plot_cells(cds,color_cells_by="cluster_ext_type",label_cell_groups=FALSE)+theme(legend.position="bottom")