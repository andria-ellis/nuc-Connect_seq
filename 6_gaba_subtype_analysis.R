set.seed(33)

#This code analyzes and visualizes data for GABAergic cell subtypes, focusing on their classification, marker gene expression, and distribution across experimental conditions. It starts by loading cell data, correcting subtype classifications, and removing unwanted cells. Key marker genes are identified and visualized, particularly Fos, to highlight cells with high activity. Finally, neuropeptide gene expression patterns and the relative proportions of each subtype across conditions are plotted, offering insights into the molecular diversity and distribution of GABA cells in different experimental settings.

# Input files: gaba_final.rds , final_GLU.rds

# Output files: cds_gaba2,rds

# Figures created: Figure 6, Figure 7A, Figure 7D, Figure 8C
#                  Supplementary Figure 4B, Supplementary Figure 5B, Supplementary Figure 3B


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

#Collect all final GLU and GABA cells 

cds_gaba=readRDS("Downloads/gaba_final.rds")
cds_glu=readRDS("Downloads/final_GLU.rds")

# There was one subset of glu cells (Glu13) that were found in the gaba cells
# Move GLU13 from Gaba to Glu

cds_glu13=cds_gaba[,cds_gaba$sub_type=="Glu13"]
cds_gaba=cds_gaba[,!cds_gaba$sub_type=="Glu13"]

#####Recombine Glu cds including Glu13

expression.matrix=cbind(exprs(cds_glu),exprs(cds_glu13))
cell.names=rbind(data.frame(pData(cds_glu)),data.frame(pData(cds_glu13)))
feature.names=join(data.frame(fData(cds_glu)),data.frame(fData(cds_glu13)),by = "gene_id",type="left",match="all")

row.names(expression.matrix) <- feature.names$gene_id
row.names(feature.names) <- feature.names$gene_id
row.names(cell.names) <- cell.names$cell_sample
colnames(expression.matrix) <- cell.names$cell_sample

cds_glu <- new_cell_data_set(expression.matrix,
                             cell_metadata = cell.names,
                             gene_metadata = feature.names)

###Remove unwanted cells

unique(cds_gaba$sub_type)

cds_gaba=cds_gaba[,cds_gaba$sub_type!="MO"]
cds_gaba=cds_gaba[,cds_gaba$sub_type!="choline transporter"]

cds_gaba[,cds_gaba$sub_type=="None"]$sub_type="Mixed/Unknown"

cds_gaba89=cds_gaba[,cds_gaba$sub_type=="GABA8/9"]
cds_gaba89[,cds_gaba89$cluster==9]$sub_type="GABA8"
cds_gaba89[,cds_gaba89$cluster==3]$sub_type="GABA9"

cds_gaba[,cds_gaba$sub_type=="GABA8/9"]$sub_type=cds_gaba89$sub_type

cds_gaba[,cds_gaba$sub_type=="Chol"]$sub_type="GABA19"
cds_gaba[,cds_gaba$sub_type=="Sox6"]$sub_type="GABA20"
cds_gaba=cds_gaba[,cds_gaba$sub_type!="GABA20"]

unique(cds_gaba$sub_type)

cds_gaba2=cds_gaba


plot_cells(cds_gaba2)

#Figure 6A: GABA subtype assignments

plot_cells(cds_gaba2[,cds_gaba2$lab=="buck"], color_cells_by="sub_type",cell_size=.8,label_cell_groups=FALSE)+
  theme(legend.position="left")+
  scale_color_manual(values=c('#666633', '#663333','#004488', '#DDAA33', '#BB5566','#FF6600', 
                              '#332288', '#DDCC77', '#117733', '#AA4499', '#33BBEE', '#EE3377', 
                              '#CC3311', '#009988', '#882255', '#FFCC00', '#999933', '#DDDDDD'))

#Supplementary Figure 3B

plot_cells(cds_gaba2[,cds_gaba2$lab=="chen"], color_cells_by="sub_type",cell_size=.8,label_cell_groups=FALSE)+
  theme(legend.position="left")+
  scale_color_manual(values=c('#666633', '#663333','#004488', '#DDAA33', '#BB5566','#FF6600', 
                              '#332288', '#DDCC77', '#117733', '#AA4499', '#33BBEE', '#EE3377', 
                              '#CC3311', '#009988', '#882255', '#FFCC00', '#999933', '#DDDDDD'))

#Figure 7D : GABA cells clustering by experiment type

plot_cells(cds_gaba2[,cds_gaba2$lab=="buck"], color_cells_by="exp",cell_size=.8,label_cell_groups=FALSE)+
  theme(legend.position="none")+
  scale_color_manual(values = c("royalblue","orangered"))

#### High Fos expression filtering

ind=which(rowData(cds_gaba2)$gene_short_name=="Fos")
cds_gaba2$fos=exprs(cds_gaba2)[ind,]*colData(cds_gaba2)$Size_Factor

cds_gaba2$high_fos=FALSE
cds_gaba2[,cds_gaba2$fos>=5]$high_fos=TRUE
plot_cells(cds_gaba2[,cds_gaba2$lab=="buck"],color_cells_by="high_fos",cell_size=.8)+
  scale_color_manual(values = c("lightgray","navy"))+
  theme(legend.position="None",
        text=element_text(size=16))

ind=which(rowData(cds_gaba2)$gene_short_name=="Fos")
cds_gaba2$fos=exprs(cds_gaba2)[ind,]*colData(cds_gaba2)$Size_Factor

cds_gaba2$high_fos=FALSE
cds_gaba2[,cds_gaba2$fos>=5]$high_fos=TRUE


#Figure 8C: High fos expression in GABA cells

plot_cells(cds_gaba2[,cds_gaba2$lab=="buck"],color_cells_by="high_fos",cell_size=.8)+
  scale_color_manual(values = c("lightgray","navy"))+
  theme(legend.position="None",
        text=element_text(size=16))


#Finding top marker genes in GABA cells

marker_test_res <- top_markers(cds_gaba2, group_cells_by="sub_type", 
                               reference_cells=1000, cores=8)
top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(3, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

np_table=read.table("Downloads/neuropeptide_list_2020.txt",sep="\t",header=TRUE)

new_np_list=c("Sst","Gal","Tac1","Ghrh","Nucb2","Penk","Pnoc","Avp","Nts","Cartpt","Nms")

chen_list=c("Pvalb","Npas1","Bcl11b","Lhx8","Pax6","Trh","Vipr2","Vip","Prok2",
            "Ghrh","Crabp1","Slc18a2","Gal","Cbln4","Cox6a2","Slc6a3","Lhx1","Klhl1")

np_list=np_table$Genes

np_list=c("Adcyap1","Adipoq","Agrp","Agt","Avp","Calca",
          "Cartpt","Cck","Crh","Cort","Gal","Galp","Gast","Ghrh",
          "Ghrl","Grp","Hcrt","Kiss1","Kng1","Nmb","Nms","Nmu","Npb",
          "Npff","Nppa","Nppc","Npy","Nts","Nucb2","Oxt","Pdyn",
          "Penk","Pmch","Pnoc","Pomc","Qrfp","Rln1",
          "Sct","Sst","Tac1","Tac2","Trh","Vip")

#Supplementary Figure 4B: Neuropeptide expression in gaba subtypes

plot_genes_by_group(cds_gaba,
                    np_list,
                    group_cells_by="sub_type",
                    ordering_type="none",
                    max.size=5)+
  scale_x_discrete(limits = c("GABA1","GABA2","GABA3","GABA5","GABA6","GABA7","GABA8","GABA9",
                              "GABA10","GABA11","GABA12","GABA13","GABA14","GABA16","GABA17","GABA18","GABA19","GABA20"))

#Figure 6B : Chen marker genes on gaba subclusters
plot_genes_by_group(cds_gaba,
                    chen_list,
                    group_cells_by="sub_type",
                    ordering_type="none",
                    max.size=5)+
  scale_x_discrete(limits = c("GABA1","GABA2","GABA3","GABA5","GABA6","GABA7","GABA8","GABA9",
                              "GABA10","GABA11","GABA12","GABA13","GABA14","GABA16","GABA17","GABA18","GABA19"))+
  theme(axis.text.x = element_text(angle = 90))+
  labs(x = "Subtype")+
  theme(legend.position = "none")

#Supplementary Figure 5B: Top marker expression in gaba subtypes

plot_genes_by_group(cds_gaba,
                    top_specific_marker_ids,
                    group_cells_by="sub_type",
                    ordering_type="maximal_on_diag",
                    max.size=3)+
  scale_x_discrete(limits = c("GABA1","GABA10","GABA11","GABA12","GABA13","GABA14","GABA16","GABA17","GABA18","GABA19",
                              "GABA2","GABA20","GABA3","GABA5","GABA6","GABA7","GABA8","GABA9"))+
  theme(axis.text.x = element_text(angle = 90))

saveRDS(cds_gaba2, file = "cds_gaba2.rds")

cds_gaba_subtypes=cds_gaba2[,cds_gaba2$sub_type!="Mixed/Unknown"]


plot_cells(cds_gaba_subtypes, color_cells_by="sub_type",cell_size=1)+theme(legend.position="left")


marker_test_res <- top_markers(cds_gaba_subtypes, group_cells_by="sub_type", 
                               reference_cells=1000, cores=8)
top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(3, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

np_table=read.table("Downloads/neuropeptide_list_2020.txt",sep="\t",header=TRUE)

new_np_list=c("Sst","Gal","Tac1","Ghrh","Nucb2","Penk","Pnoc","Avp","Nts","Cartpt","Nms")

chen_marker_list=c("Sox6","Pvalb","Npas1","Bcl11b","Lhx8","Pax6","Trh","Vipr2","Avp","Vip","Prok2","Ghrh","Crabp1","Slc18a2","Gal","Cbln4","Cox6a2","Slc6a3","Klhl1","Sst","Penk","Pnoc","Cartpt","Nts","Tac1","Nms","Coch","Chat","Lhx6","Lhx8","Prima1","Rarb")

plot_genes_by_group(cds_gaba_subtypes,
                    chen_marker_list,
                    group_cells_by="sub_type",
                    ordering_type="maximal_on_diag",
                    max.size=5)

### Figure 6C GABA subtype percentage

df=data.frame(table(cds_gaba$lab,cds_gaba$sub_type))

df$num_cells=0

a=sum(df[df$Var1=="buck",]$Freq)
b=sum(df[df$Var1=="chen",]$Freq)

df[df$Var1=="buck",]$num_cells=a
df[df$Var1=="chen",]$num_cells=b

df$percent=df$Freq/df$num_cells*100

ggplot(data=df, aes(x=Var2, y=percent, fill=Var1)) +
  theme_bw()+
  ylab("Percentage (%)") +
  geom_bar(stat="identity", position=position_dodge())+
  theme(axis.text.x = element_text(angle = 90,vjust=.5))+
  scale_x_discrete(limits = c("GABA1","GABA2","GABA3","GABA5","GABA6","GABA7","GABA8","GABA9",
                              "GABA10","GABA11","GABA12","GABA13","GABA14","GABA16","GABA17","GABA18","GABA19"))+
  scale_fill_manual(values = c("steelblue3","navy"))+
  theme(legend.position="None",
        axis.title.x=element_blank(),
        text=element_text(size=16))+
  scale_y_continuous(limits=c(0, 15))



ggplot(data=df, aes(x=Var2, y=percent, fill=Var1)) +
  theme_bw()+
  ylab("Percentage (%)") +
  geom_bar(stat="identity", position=position_dodge())+
  theme(axis.text.x = element_text(angle = 90,hjust=0.5))+
  scale_x_discrete(limits = c("Glu1","Glu2","Glu3","Glu4","Glu5","Glu6","Glu7","Glu8",
                              "Glu9","Glu10","Glu11","Glu12","Glu13","Glu14","Glu15","Glu16","Glu17","Glu18"))+
  scale_fill_manual(values = c("steelblue3","navy"))+
  theme(legend.position="None",
        axis.title.x=element_blank(),
        text=element_text(size=16))+
  scale_y_continuous(limits=c(0, 30))