sessionInfo()

# The purpose of this code block is to customize the UMAP plots to display clusters by cell type and lab of origin. Genes indicative of certain cell types, such as oligodendrocytes, astrocytes, microglia, and neurons, are visualized with distinct color gradients. The code also annotates clusters with assigned cell types, compares cell-type composition between Chen and Buck samples, and plots the proportions of each cell type across both datasets. Lastly, it isolates the neuron subset and saves this specific dataset for further analysis, allowing targeted exploration of neuronal gene expression patterns.

# Input files: integrated_cds.chen.rds

# Output files: cds_neurons_withchen.rds



#Figures created: Figure 3 and 4

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

# Read in the integrated CDS
cds=readRDS("integrated_cds.chen.rds")


plot_cells(cds,color_cells_by="assigned_cell_type", label_cell_groups = FALSE,cell_size = 1)+
  scale_colour_manual(values = c("orange2","seagreen3","dodgerblue2","orangered","darkmagenta","orchid","royalblue4")) +
  theme(legend.position="None",
        axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        text=element_text(size=16))

plot_cells(cds,color_cells_by="lab", label_cell_groups = FALSE,cell_size = 1)+
  scale_colour_manual(values = c("steelblue3","navy"))+
  theme(legend.position="None",
        axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        text=element_text(size=16))

cds$cell="temp"

plot_cells(cds,color_cells_by="cell",label_cell_groups = FALSE)+
  scale_colour_manual(values = c("grey"))

## The following plots are used to show known marker labels for various cell types that were collected. We combine garnett labeling with markers to make determinations on the final cell types of each cluster. Figure 4 is made from the following marker label plots. 

#oligodendrycytes
p1=plot_cells(cds,genes=c("Olig1"),label_groups_by_cluster=FALSE,label_cell_groups = FALSE,cell_size=1) +
  scale_colour_gradientn(colours = c("#CCCCCC","#CCCCCC","#000333")) +
  theme(legend.position="None",
        axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        text=element_text(size=50))

ggsave(p1,filename="chen.olig1.png")

p2=plot_cells(cds,genes=c("Top2a"),label_groups_by_cluster=FALSE,label_cell_groups = FALSE,cell_size=1) +
  scale_colour_gradientn(colours = c("#CCCCCC","#CCCCCC","#000333")) +
  theme(legend.position="None",
        axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        text=element_text(size=50))

ggsave(p2,filename="chen.top2a.png")

p3=plot_cells(cds,genes=c("Pdgfra"),label_groups_by_cluster=FALSE,label_cell_groups = FALSE,cell_size=1) +
  scale_colour_gradientn(colours = c("#CCCCCC","#CCCCCC","#000333")) +
  theme(legend.position="None",
        axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        text=element_text(size=50))

ggsave(p3,filename="chen.pdgfra.png")

#astrocytes
p4=plot_cells(cds,genes=c("Aldh1l1"),label_groups_by_cluster=FALSE,label_cell_groups = FALSE,cell_size=1) +
  scale_colour_gradientn(colours = c("#CCCCCC","#CCCCCC","#000333")) +
  theme(legend.position="None",
        axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        text=element_text(size=50))

ggsave(p4,filename="chen.aldh1l1.png")

p5=plot_cells(cds,genes=c("Gfap"),label_groups_by_cluster=FALSE,label_cell_groups = FALSE,cell_size=1) +
  scale_colour_gradientn(colours = c("#CCCCCC","#CCCCCC","#000333"))+
  theme(legend.position="None",
        axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        text=element_text(size=50))

ggsave(p5,filename="chen.gfap.png")

p6=plot_cells(cds,genes=c("Aldoc"),label_groups_by_cluster=FALSE,label_cell_groups = FALSE,cell_size=1) +
  scale_colour_gradientn(colours = c("#CCCCCC","#CCCCCC","#000333"))+
  theme(legend.position="None",
        axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        text=element_text(size=50))

ggsave(p6,filename="chen.aldoc.png")

#microglia
p7=plot_cells(cds,genes=c("Aif1"),label_groups_by_cluster=FALSE,label_cell_groups = FALSE,cell_size=1) +
  scale_colour_gradientn(colours = c("#CCCCCC","#CCCCCC","#000333")) +
  theme(legend.position="None",
        axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        text=element_text(size=50))

ggsave(p7,filename="chen.aif1.png")

p8=plot_cells(cds,genes=c("Hexb"),label_groups_by_cluster=FALSE,label_cell_groups = FALSE,cell_size=1) +
  scale_colour_gradientn(colours = c("#CCCCCC","#CCCCCC","#000333")) +
  theme(legend.position="None",
        axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        text=element_text(size=50))

ggsave(p8,filename="chen.hexb.png")

p9=plot_cells(cds,genes=c("Mafb"),label_groups_by_cluster=FALSE,label_cell_groups = FALSE,cell_size=1) +
  scale_colour_gradientn(colours = c("#CCCCCC","#CCCCCC","#000333")) +
  theme(legend.position="None",
        axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        text=element_text(size=50))

ggsave(p9,filename="chen.mafb.png")

#neurons

p10=plot_cells(cds,genes=c("Syt1"),label_groups_by_cluster=FALSE,label_cell_groups = FALSE,cell_size=1) +
  scale_colour_gradientn(colours = c("#CCCCCC","#CCCCCC","#000333")) +
  theme(legend.position="None",
        axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        text=element_text(size=50))

ggsave(p10,filename="chen.syt1.png")

p11=plot_cells(cds,genes=c("Map2"),label_groups_by_cluster=FALSE,label_cell_groups = FALSE,cell_size=1) +
  scale_colour_gradientn(colours = c("#CCCCCC","#CCCCCC","#000333")) +
  theme(legend.position="None",
        axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        text=element_text(size=50))

ggsave(p11,filename="chen.map2.png")

p12=plot_cells(cds,genes=c("Snap25"),label_groups_by_cluster=FALSE,label_cell_groups = FALSE,cell_size=1) +
  scale_colour_gradientn(colours = c("#CCCCCC","#CCCCCC","#000333")) +
  theme(legend.position="None",
        axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        text=element_text(size=50))

ggsave(p12,filename="chen.snap25.png")

p13=plot_cells(cds,genes=c("Igfbp7"),label_groups_by_cluster=FALSE,label_cell_groups = FALSE,cell_size=1) +
  scale_colour_gradientn(colours = c("#CCCCCC","#CCCCCC","#000333")) +
  theme(legend.position="None",
        axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        text=element_text(size=50))

ggsave(p13,filename="chen.igfbp7.png")

p14=plot_cells(cds,genes=c("Rax"),label_groups_by_cluster=FALSE,label_cell_groups = FALSE,cell_size=1) +
  scale_colour_gradientn(colours = c("#CCCCCC","#CCCCCC","#000333")) +
  theme(legend.position="None",
        axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        text=element_text(size=50))

ggsave(p14,filename="chen.rax.png")

p14=plot_cells(cds,genes=c("Ccdc153"),label_groups_by_cluster=FALSE,label_cell_groups = FALSE,cell_size=1) +
  scale_colour_gradientn(colours = c("#CCCCCC","#CCCCCC","#000333")) +
  theme(legend.position="None",
        axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        text=element_text(size=50))

ggsave(p14,filename="chen.ccdc153.png")


plot_cells(cds,genes=c("Fyn"),label_groups_by_cluster=FALSE,label_cell_groups = FALSE,cell_size=1) +
  scale_colour_gradientn(colours = c("#CCCCCC","#CCCCCC","#000333")) +
  theme(legend.position="None",
        axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        text=element_text(size=50))

ggsave(p15,filename="chen.fyn.png")

p16=plot_cells(cds,genes=c("Ctss"),label_groups_by_cluster=FALSE,label_cell_groups = FALSE,cell_size=1) +
  scale_colour_gradientn(colours = c("#CCCCCC","#CCCCCC","#000333")) +
  theme(legend.position="None",
        axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        text=element_text(size=50))

ggsave(p16,filename="chen.ctss.png")


colData(cds)$cluster<- as.character(clusters(cds))


#Final cluster assignments, based on assigned cell type and marker labels
colData(cds)$assigned_cell_type = dplyr::recode(colData(cds)$cluster,
                                                "1"="Microglia",
                                                "2"="Oligodendrocyte",
                                                "3"="Neuron",
                                                "4"="Neuron",
                                                "5"="Neuron",
                                                "6"="O.P.C.",
                                                "7"="Astrocyte",
                                                "8"="Epithelial",
                                                "9"="Neuron",
                                                "10"="Microglia",
                                                "11"="Oligodendrocyte",
                                                "12"="Neuron",
                                                "13"="Tanycyte",
                                                "14"="Neuron",
                                                "15"="Neuron",
                                                "16"="Neuron",
                                                "17"="Ependycyte",
                                                "18"="P.O.P.C.",
                                                "19"="Neuron",
                                                "20"="Epithelial",
                                                "21"="Neuron",
                                                "22"="Microglia",
                                                "23"="Microglia",
                                                "24"="Epithelial",
                                                "25"="Neuron",
                                                "26"="Epithelial",
                                                "27"="Neuron",
                                                "28"="O.P.C.",
                                                "29"="Neuron",
                                                "30"="Neuron",
                                                "31"="Neuron")


# Visualization of final cell types: Figure 3A
plot_cells(cds,color_cells_by="assigned_cell_type",cell_size=1,label_cell_groups = FALSE)+theme(legend.position="none")+
  scale_colour_manual(values = c("orange2","coral","sienna","green4",
                                 "dodgerblue2","orangered","darkmagenta","orchid","deeppink")) 

# Distribution of cells per lab in the UMAP plot: Figure 3B
plot_cells(cds,color_cells_by="lab",cell_size=1,label_cell_groups = FALSE)+theme(legend.position="none")+
  scale_colour_manual(values = c("steelblue3","navy")) 


c=which(cds$lab=="chen")
b=which(cds$lab=="buck")
chendf=data.frame(table(cds[,c]$assigned_cell_type))
buckdf=data.frame(table(cds[,b]$assigned_cell_type))

chendf$percent=chendf$Freq/length(c)*100
buckdf$percent=buckdf$Freq/length(b)*100

colnames(chendf)=c("Cell_Type","Chen_Number","Chen_Percent")
colnames(buckdf)=c("Cell_Type","Buck_Number","Buck_Percent")

chendf$Buck_Number=buckdf$Buck_Number
chendf$Buck_Percent=buckdf$Buck_Percent

newdf=data.frame(chendf$Cell_Type,chendf$Buck_Percent,chendf$Chen_Percent)

dat=melt(newdf)

#Percentage of each cell type by lab: Figure 3C
p<-ggplot(dat, aes(x=chendf.Cell_Type, y=value, fill=variable)) +
  geom_bar(stat="identity",position=position_dodge())+ 
  theme_minimal() + 
  scale_fill_manual(values = c("steelblue3","navy"))+
  scale_y_continuous(name="Percentage (%)", limits=c(0, 60))+
  theme(legend.position="None",
        axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        text=element_text(size=20),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

p +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


cds_neurons=cds[,cds$assigned_cell_type=="Neurons"]

#Finsl CDS contains only neuron cells
saveRDS(cds_neurons,file="cds_neurons_withchen.rds")