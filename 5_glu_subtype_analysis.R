set.seed(33)

#This code processes and visualizes single-cell gene expression data, focusing on GABAergic (GABA) and glutamatergic (GLU) neurons. The script uses Seurat and monocle3 to perform cell clustering, dimensionality reduction, and subtype-specific marker gene analysis. Initially, GLU13 cells identified in the GABA dataset are moved to the GLU dataset to maintain subtype integrity. The combined GLU dataset undergoes preprocessing, dimensionality reduction (UMAP), and filtering to remove non-relevant subtypes. Specific cells are clustered, annotated, and visualized, with focus on high-Fos-expressing cells. Key marker genes for subtype identification and neuropeptide genes are plotted by group, highlighting expression patterns across GLU subtypes.

# Input files: gaba_final.rds , final_GLU.rds

# Output files: cds_glu2.rds

# Figures: Figure 5 and Figure 7B, Figure 7C , Figure 8B
#           Supplmentary Figure 4A, Supplmentary Figure 5A, Supplementary Figure 3B

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

cds_gaba=readRDS("gaba_final.rds")
cds_glu=readRDS("final_GLU.rds")

cds_glu$xcoord=cds_glu$umap1
cds_glu$ycoord=cds_glu$umap2

# There was one subset of glu cells (Glu13) that were found in the gaba cells
# Move GLU13 from Gaba to Glu

cds_glu13=cds_gaba[,cds_gaba$sub_type=="Glu13"]

cds_glu13$xcoord=cds_glu13$umap1
cds_glu13$ycoord=cds_glu13$umap2

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

cds_glu=preprocess_cds(cds_glu)
cds_glu=reduce_dimension(cds_glu)

###Remove unwanted cells

unique(cds_glu$sub_type)

cds_glu=cds_glu[,!cds_glu$sub_type=="Ctss"]
cds_glu=cds_glu[,!cds_glu$sub_type=="Non-neuronal"]
cds_glu=cds_glu[,!cds_glu$sub_type=="Hista"]

cds_glu[,cds_glu$sub_type=="Clu4"]$sub_type="Mixed"
cds_glu[,cds_glu$sub_type=="Clu5"]$sub_type="Mixed"
cds_glu[,cds_glu$sub_type=="Dopa"]$sub_type="Glu16"
cds_glu[,cds_glu$sub_type=="Sub_Mixed"]$sub_type="Mixed"

unique(cds_glu$sub_type)

cds_glu2=cds_glu

saved_coords =as.matrix(data.frame(cds_glu$xcoord,cds_glu$ycoord))

reducedDims(cds_glu2)[["UMAP"]] = saved_coords

cds_glu2 <- cluster_cells(cds_glu2,resolution=1e-3)

colData(cds_glu2)$cluster<- as.character(clusters(cds_glu2))

plot_cells(cds_glu2)

cds_glu_buck=cds_glu2[,cds_glu2$lab=="buck"]
cds_glu_chen=cds_glu2[,cds_glu2$lab=="chen"]

plot_cells(cds_glu2, color_cells_by="sub_type")
cds_glu_temp=cds_glu2[,cds_glu2$sub_type!="Unknown"]

ind=which(rowData(cds_glu2)$gene_short_name=="Fos")
cds_glu2$fos=exprs(cds_glu2)[ind,]*colData(cds_glu2)$Size_Factor

cds_glu2$high_fos=FALSE
cds_glu2[,cds_glu2$fos>=5]$high_fos=TRUE

saveRDS(cds_glu2, file = "cds_glu2.rds")


#Figure 5A: Glu cells colored by subtype
plot_cells(cds_glu_buck, color_cells_by="sub_type",cell_size=.8,label_cell_groups=FALSE)+
  theme(legend.position="left")+
  scale_color_manual(values=c('#666633', '#663333','#004488', '#DDAA33', '#BB5566','#FF6600', 
                              '#332288', '#DDCC77', '#117733', '#882255','#88CCEE','#EE3377', 
                              '#CC3311', '#009988',  '#AA4499', '#DDDDDD'))

#Supplmentary Figure 3B

plot_cells(cds_glu_chen, color_cells_by="sub_type",cell_size=.8,label_cell_groups=FALSE)+
  theme(legend.position="left")+
  scale_color_manual(values=c('#666633', '#663333','#004488', '#DDAA33', '#BB5566','#FF6600', 
                              '#332288', '#DDCC77', '#117733', '#882255','#88CCEE','#EE3377', 
                              '#CC3311', '#009988',  '#AA4499', '#DDDDDD','#DDDDDD'))

# Figure 7C: Glu cell clustering by experiment group
plot_cells(cds_glu_buck, color_cells_by="exp",cell_size=.8,label_cell_groups=FALSE)+
  theme(legend.position="none")+
  scale_color_manual(values = c("royalblue","orangered"))

####
ind=which(rowData(cds_glu_buck)$gene_short_name=="Fos")
cds_glu_buck$fos=exprs(cds_glu_buck)[ind,]*colData(cds_glu_buck)$Size_Factor

cds_glu_buck$high_fos=FALSE
cds_glu_buck[,cds_glu_buck$fos>=5]$high_fos=TRUE


#Figure 8B: High Fos Glu Cells
plot_cells(cds_glu_buck[,cds_glu_buck$lab=="buck"],color_cells_by="high_fos",cell_size=.8)+
  scale_color_manual(values = c("lightgray","navy"))+
  theme(legend.position="None",
        text=element_text(size=16))


#Marker Genes for chen subtypes

chen_marker_list=c("Crh","Fezf2","Sln","Samd3","Shox2","Foxb1",
                   "Tac1","Fezf1","Cartpt","Gng8","Trh","Kiss1","Vgll2","Pomc","Avp","Sim1","Sst","Prdm8","Ddc")

#Figure 5B : Expression of marker genes in subtype cells
plot_genes_by_group(cds_glu_temp,
                    chen_marker_list,
                    group_cells_by="sub_type",
                    ordering_type="none",
                    max.size=5)+
  scale_x_discrete(limits = c("Glu1","Glu2","Glu3","Glu4","Glu5","Glu6","Glu7","Glu8",
                              "Glu9","Glu10","Glu11","Glu12","Glu13","Glu14","Glu15","Glu16"))+
  theme(axis.text.x = element_text(angle = 90))+
  labs(x = "Subtype")+
  theme(legend.position = "none")

#Supplmentary Figure 5A : Top differentially expreseed genes in glu cell types
plot_genes_by_group(cds_glu_temp,
                    top_specific_marker_ids,
                    group_cells_by="sub_type",
                    ordering_type="maximal_on_diag",
                    max.size=5)+
  scale_x_discrete(limits = c("Glu1","Glu2","Glu3","Glu4","Glu5","Glu6","Glu7","Glu8",
                              "Glu9","Glu10","Glu11","Glu12","Glu13","Glu14","Glu15","Glu16"))


np_list=np_table$Genes

np_list=c("Adcyap1","Adipoq","Agrp","Agt","Avp","Calca",
          "Cartpt","Cck","Crh","Cort","Gal","Galp","Gast","Ghrh",
          "Ghrl","Grp","Hcrt","Kiss1","Kng1","Nmb","Nms","Nmu","Npb",
          "Npff","Nppa","Nppc","Npy","Nts","Nucb2","Oxt","Pdyn",
          "Penk","Pmch","Pnoc","Pomc","Qrfp","Rln1",
          "Sct","Sst","Tac1","Tac2","Trh","Vip")


# Supplmentary Figure 4A : Neuropeptide expression in glu cell types
plot_genes_by_group(cds_glu_temp,
                    np_list,
                    group_cells_by="sub_type",
                    ordering_type="none",
                    max.size=5)+
  scale_x_discrete(limits = c("Glu1","Glu2","Glu3","Glu4","Glu5","Glu6","Glu7","Glu8",
                              "Glu9","Glu10","Glu11","Glu12","Glu13","Glu14","Glu15","Glu16"))



marker_test_res <- top_markers(cds_glu, group_cells_by="sub_type", 
                               reference_cells=1000, cores=8)
top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(3, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

np_table=read.table("Downloads/neuropeptide_list_2020.txt",sep="\t",header=TRUE)

new_np_list=c("Crh","Gal","Tac1","Tac2","Hcrt","Pomc","Sst","Avp","Oxt","Cartpt","Nts","Kiss1","Trh")

plot_genes_by_group(cds_glu2,
                    np_table$Genes,
                    group_cells_by="sub_type",
                    ordering_type="maximal_on_diag",
                    max.size=5)

plot_genes_by_group(cds_glu2,
                    top_specific_marker_ids,
                    group_cells_by="sub_type",
                    ordering_type="maximal_on_diag",
                    max.size=3)

cds_glu_subtypes=cds_glu2[,cds_glu2$sub_type!="Mixed"]


library(data.table)

## Figure 5C : Percentage glu cell types in experiment vs in hypothalamus

df=data.frame(table(cds_glu$lab,cds_glu$sub_type))

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
  scale_x_discrete(limits = c("Glu1","Glu2","Glu3","Glu4","Glu5","Glu6","Glu7","Glu8",
                              "Glu9","Glu10","Glu11","Glu12","Glu13","Glu14","Glu15","Glu16"))+
  scale_fill_manual(values = c("steelblue3","navy"))+
  theme(legend.position="None",
        axis.title.x=element_blank(),
        text=element_text(size=16))+
  scale_y_continuous(limits=c(0, 30))


cds_glu_buck=cds_glu[,cds_glu$lab=="buck"]


cds=cds_glu_buck


ind=which(rowData(cds)$gene_short_name=="Fos")

cds$fos_level="ANone"
cds$fos=exprs(cds)[ind,]*colData(cds)$Size_Factor

cds[,cds$fos>0]$fos_level="BLow"

cds[,cds$fos>2]$fos_level="CMedium"

cds[,cds$fos>5]$fos_level="DHigh"

plot_cells(cds,color_cells_by="fos_level")

cds$binary=0
cds[,cds$exp=="Restraint"]$binary=1

df_glu=data.frame(cds$fos_level,cds$binary)

df_glu$cds.binary=as.factor(df_glu$cds.binary)
df_glu$cds.fos_level=as.factor(df_glu$cds.fos_level)
model <- glm(cds.binary ~.,family=binomial(link='logit'),data=df_glu)




###DEG testing

cds_glu_subtypes$glu4_testing="NO"

cds_glu_subtypes[,cds_glu_subtypes$sub_type=="Glu4"]$glu4_testing="YES"

cell_groups = c("YES", "NO")

sub_cds = cds_glu_subtypes
sub_cds = estimate_size_factors(sub_cds)
sub_cds = detect_genes(sub_cds)
sub_cds = sub_cds[Matrix::rowSums(SingleCellExperiment::counts(sub_cds) > 0) > 5,]
fits = fit_models(sub_cds, model_formula_str = "~glu4_testing", cores = 1) 
# the model formula string here is the style of group that you want to contrast
mod_coef = coefficient_table(fits)

mod_coef = mod_coef %>% dplyr::mutate(up_in = case_when(
  estimate < 0 ~ "NO",
  estimate > 0 ~ "YES"))

mod_coef %>% dplyr::filter(term != "(Intercept)") %>% dplyr::filter(q_value < 0.05 ) %>%  
  dplyr::arrange(up_in, q_value) %>%
  dplyr::select(up_in, gene_id, gene_short_name, q_value, estimate) %>%
  fwrite("glu4_testing_DEG_list.csv", row.names = FALSE, sep = ",", na = "NA")