set.seed(33)

# This R code processes and analyzes gene expression data for Glu and GABA cell datasets. It computes expression levels for the gene "Fos," categorizes cells based on Fos expression levels, and visualizes the results. A generalized linear model (GLM) is employed to analyze the relationship between binary experimental conditions (Restraint vs. Control) and Fos levels. The code also evaluates the expression of neuropeptides across different Fos expression levels, generating multiple plots to summarize the results.

#Input files: cds_glu2.rds,cds_gaba2.rds

#Figures created: Figure 8 A,B , Figure 7A,D,E

## You must have already run all previous code block scripts to make the GLU and GABA CDS's

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


# Load in final glu/gaba groups

cds_glu2 <- readRDS("cds_glu2.rds")
cds_gaba2 <- readRDS("cds_gaba2.rds")

pData(cds_glu2)=pData(cds_glu2)[ , -which(names(pData(cds_glu2)) %in% c("xcoord","ycoord"))]

colData(cds_gaba2)
expression.matrix=cbind(exprs(cds_glu2),exprs(cds_gaba2))
cell.names=rbind(data.frame(pData(cds_glu2)),data.frame(pData(cds_gaba2)))
feature.names=join(data.frame(fData(cds_glu2)),data.frame(fData(cds_gaba2)),by = "gene_id",type="left",match="all")

row.names(expression.matrix) <- feature.names$gene_id
row.names(feature.names) <- feature.names$gene_id
row.names(cell.names) <- cell.names$cell_sample
colnames(expression.matrix) <- cell.names$cell_sample

cds_new <- new_cell_data_set(expression.matrix,
                             cell_metadata = cell.names,
                             gene_metadata = feature.names)

cds_new=cds_new[,cds_new$lab=="buck"]

#Filtering based on fos level

cds=cds_new

ind=which(rowData(cds)$gene_short_name=="Fos")

cds$fos_level="ANone"
cds$fos=exprs(cds)[ind,]*colData(cds)$Size_Factor

cds[,cds$fos>0]$fos_level="BLow"
cds[,cds$fos>1]$fos_level="CLow"

cds[,cds$fos>2]$fos_level="DMedium"
cds[,cds$fos>3]$fos_level="EMedium"
cds[,cds$fos>4]$fos_level="FMedium"

cds[,cds$fos>5]$fos_level="GHigh"


cds_gaba2$fos=exprs(cds_gaba2)[ind,]*colData(cds_gaba2)$Size_Factor

cds_gaba2$high_fos=FALSE
cds_gaba2[,cds_gaba2$fos>=5]$high_fos=TRUE
plot_cells(cds_gaba2,color_cells_by="high_fos")


cds$binary=0
cds[,cds$exp=="Restraint"]$binary=1

df_glu=data.frame(cds$fos_level,cds$binary)

#GLM to model fos level across glue cells

df_glu$cds.binary=as.factor(df_glu$cds.binary)
df_glu$cds.fos_level=as.factor(df_glu$cds.fos_level)
model <- glm(cds.binary ~.,family=binomial(link='logit'),data=df_glu)

summary(model)

estimate=c(0.036,0.27,0.28,0.61,1.26)
fos=c(1,2,3,4,5)
cells=c(248,218,185,159,133)

df_temp=data.frame(estimate,fos,cells)

# Figure 8A
ggplot(data=df_temp, aes(x=fos, y=estimate)) +
  geom_bar(stat="identity") +theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        text=element_text(size=16))

ggplot(data=df_temp, aes(x=fos, y=cells)) +
  geom_bar(stat="identity")+theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        text=element_text(size=16))

#Filtering FOS expression for gaba cells

cds_glu_buck=cds_glu2[,cds_glu2$lab=="buck"]
cds_gaba_buck=cds_gaba2[,cds_gaba2$lab=="buck"]
cds_new_buck=cds_new[,cds_new$lab=="buck"]

cds=cds_gaba_buck

cds$fos=exprs(cds)[ind,]*colData(cds)$Size_Factor


cds_fos1=cds[,cds$fos>=1]
cds_fos2=cds[,cds$fos>=2]
cds_fos3=cds[,cds$fos>=3]
cds_fos4=cds[,cds$fos>=4]
cds_fos5=cds[,cds$fos>=5]

cds_fos5
cds_fos0=cds

plot_df0=as.data.frame(table(cds_fos0$exp,cds_fos0$sub_type))
plot_df5=as.data.frame(table(cds_fos5$exp,cds_fos5$sub_type))


#Figure 7A: Gaba subtypes by condition
ggplot(plot_df0, aes(fill=Var1, y=Freq, x=Var2)) + 
  theme_bw()+
  geom_bar(position="dodge", stat="identity")+
  scale_fill_manual(values = c("royalblue","orangered"))+
  ylim(0, 150)+
  theme(plot.title = element_text(hjust = 0.5,size=24,face="bold"), 
        text=element_text(size=24),
        legend.position="none",
        axis.title.x=element_blank())+
  scale_x_discrete(limits = c("GABA1","GABA2","GABA3","GABA5","GABA6","GABA7","GABA8","GABA9","GABA10","GABA11","GABA12","GABA13","GABA14","GABA15","GABA16","GABA17","GABA18","GABA19"))+
  theme(axis.text.x = element_text(angle = 90,vjust=.5))+
  labs(y = "Cell Counts")

# Figure 8D : High Fos gaba cells by condition
ggplot(plot_df5, aes(fill=Var1, y=Freq, x=Var2)) + 
  theme_bw()+
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values = c("royalblue","orangered"))+
  theme(plot.title = element_text(hjust = 0.5,size=24,face="bold"), 
        text=element_text(size=24),
        legend.position="none",
        axis.title.x=element_blank())+
  scale_x_discrete(limits = c("GABA2","GABA3","GABA5","GABA6","GABA9","GABA12","GABA13","GABA14","GABA16","GABA19","Mixed/Unknown"),
                   labels = c("GABA2","GABA3","GABA5","GABA6","GABA9","GABA12","GABA13","GABA14","GABA16","GABA19","Mixed"))+
  theme(axis.text.x = element_text(angle = 90,vjust=.5))+
  labs(y = "Cell Counts")

#Filtering high fos cells for glu subtypes

cds=cds_glu_buck

cds$fos=exprs(cds)[ind,]*colData(cds)$Size_Factor


cds_fos1=cds[,cds$fos>=1]
cds_fos2=cds[,cds$fos>=2]
cds_fos3=cds[,cds$fos>=3]
cds_fos4=cds[,cds$fos>=4]
cds_fos5=cds[,cds$fos>=5]

cds_fos5
cds_fos0=cds

plot_df0=as.data.frame(table(cds_fos0$exp,cds_fos0$sub_type))
plot_df5=as.data.frame(table(cds_fos5$exp,cds_fos5$sub_type))


#Figure 7B: Gaba subtypes by condition
ggplot(plot_df0, aes(fill=Var1, y=Freq, x=Var2)) + 
  theme_bw()+
  geom_bar(position="dodge", stat="identity")+
  scale_fill_manual(values = c("royalblue","orangered"))+
  ylim(0, 650)+
  theme(plot.title = element_text(hjust = 0.5,size=24,face="bold"), 
        text=element_text(size=24),
        legend.position="none",
        axis.title.x=element_blank())+
  scale_x_discrete(limits = c("Glu2","Glu3","Glu4","Glu5","Glu6","Glu7","Glu8","Glu9","Glu10","Glu11","Glu12","Glu13","Glu14","Glu15","Glu16"))+
  theme(axis.text.x = element_text(angle = 90,vjust=.5))+
  labs(y = "Cell Counts")

# Figure 8D: High fos gaba subtypes by condition
ggplot(plot_df5, aes(fill=Var1, y=Freq, x=Var2)) + 
  theme_bw()+
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values = c("royalblue","orangered"))+
  theme(plot.title = element_text(hjust = 0.5,size=24,face="bold"), 
        text=element_text(size=24),
        legend.position="none",
        axis.title.x=element_blank())+
  scale_x_discrete(limits = c("Glu3","Glu4","Glu7","Glu8","Glu14","Glu16","Mixed"))+
  theme(axis.text.x = element_text(angle = 90,vjust=.5))+
  labs(y = "Cell Counts")


## Run everything below this to produce FIG8E

np_table=read.table("Downloads/neuropeptide_list_2020.txt",sep="\t",header=TRUE)

np_table$Genes=as.character(np_table$Genes)
np_list=np_table$Genes

restraint=c()
total=c()
control=c()

cds_control=cds[,which(cds$exp=="Control")]
cds_restraint=cds[,which(cds$exp=="Restraint")]

CUTOFF=0

for (i in 1:length(np_list)){
  total_cells_express=length(which(exprs(cds)[which(fData(cds)$gene_short_name==np_list[i]),]>CUTOFF))
  control_cell_express=length(which(exprs(cds_control)[which(fData(cds_control)$gene_short_name==np_list[i]),]>CUTOFF))
  restraint_cell_express=length(which(exprs(cds_restraint)[which(fData(cds_restraint)$gene_short_name==np_list[i]),]>CUTOFF))
  
  control=c(control, control_cell_express)
  restraint=c(restraint,restraint_cell_express)
  total=c(total,total_cells_express)
}


control_percent=control/length(which(cds$exp=="Control"))*100
restraint_percent=restraint/length(which(cds$exp=="Restraint"))*100
total_percent=total/length(cds$exp)*100

counts=c(control_percent,restraint_percent,total_percent)
np=rep(np_list,3)

np_ind=length(np_list)
treatment=c(rep("control",np_ind),rep("restaint",np_ind),rep("total",np_ind))
fos=rep(0,3*np_ind)

df_fos0=data.frame(np,counts,treatment,fos)


cds=cds_fos1
cds_control=cds[,which(cds$exp=="Control")]
cds_restraint=cds[,which(cds$exp=="Restraint")]

restraint=c()
total=c()
control=c()

CUTOFF=0

for (i in 1:length(np_list)){
  total_cells_express=length(which(exprs(cds)[which(fData(cds)$gene_short_name==np_list[i]),]>CUTOFF))
  control_cell_express=length(which(exprs(cds_control)[which(fData(cds_control)$gene_short_name==np_list[i]),]>CUTOFF))
  restraint_cell_express=length(which(exprs(cds_restraint)[which(fData(cds_restraint)$gene_short_name==np_list[i]),]>CUTOFF))
  
  control=c(control, control_cell_express)
  restraint=c(restraint,restraint_cell_express)
  total=c(total,total_cells_express)
}

control_percent=control/length(which(cds$exp=="Control"))*100
restraint_percent=restraint/length(which(cds$exp=="Restraint"))*100
total_percent=total/length(cds$exp)*100

counts=c(control_percent,restraint_percent,total_percent)
np=rep(np_list,3)

np_ind=length(np_list)
treatment=c(rep("control",np_ind),rep("restaint",np_ind),rep("total",np_ind))
fos=rep(1,3*np_ind)

df_fos1=data.frame(np,counts,treatment,fos)


###############################################

cds=cds_fos2
cds_control=cds[,which(cds$exp=="Control")]
cds_restraint=cds[,which(cds$exp=="Restraint")]

restraint=c()
total=c()
control=c()

CUTOFF=0

for (i in 1:length(np_list)){
  total_cells_express=length(which(exprs(cds)[which(fData(cds)$gene_short_name==np_list[i]),]>CUTOFF))
  control_cell_express=length(which(exprs(cds_control)[which(fData(cds_control)$gene_short_name==np_list[i]),]>CUTOFF))
  restraint_cell_express=length(which(exprs(cds_restraint)[which(fData(cds_restraint)$gene_short_name==np_list[i]),]>CUTOFF))
  
  control=c(control, control_cell_express)
  restraint=c(restraint,restraint_cell_express)
  total=c(total,total_cells_express)
}

control_percent=control/length(which(cds$exp=="Control"))*100
restraint_percent=restraint/length(which(cds$exp=="Restraint"))*100
total_percent=total/length(cds$exp)*100

counts=c(control_percent,restraint_percent,total_percent)
np=rep(np_list,3)

np_ind=length(np_list)
treatment=c(rep("control",np_ind),rep("restaint",np_ind),rep("total",np_ind))
fos=rep(2,3*np_ind)

df_fos2=data.frame(np,counts,treatment,fos)

##################

cds=cds_fos3
cds_control=cds[,which(cds$exp=="Control")]
cds_restraint=cds[,which(cds$exp=="Restraint")]

restraint=c()
total=c()
control=c()

CUTOFF=0

for (i in 1:length(np_list)){
  total_cells_express=length(which(exprs(cds)[which(fData(cds)$gene_short_name==np_list[i]),]>CUTOFF))
  control_cell_express=length(which(exprs(cds_control)[which(fData(cds_control)$gene_short_name==np_list[i]),]>CUTOFF))
  restraint_cell_express=length(which(exprs(cds_restraint)[which(fData(cds_restraint)$gene_short_name==np_list[i]),]>CUTOFF))
  
  control=c(control, control_cell_express)
  restraint=c(restraint,restraint_cell_express)
  total=c(total,total_cells_express)
}

control_percent=control/length(which(cds$exp=="Control"))*100
restraint_percent=restraint/length(which(cds$exp=="Restraint"))*100
total_percent=total/length(cds$exp)*100

counts=c(control_percent,restraint_percent,total_percent)
np=rep(np_list,3)

np_ind=length(np_list)
treatment=c(rep("control",np_ind),rep("restaint",np_ind),rep("total",np_ind))
fos=rep(3,3*np_ind)

df_fos3=data.frame(np,counts,treatment,fos)

############

cds=cds_fos4
cds_control=cds[,which(cds$exp=="Control")]
cds_restraint=cds[,which(cds$exp=="Restraint")]

restraint=c()
total=c()
control=c()


for (i in 1:length(np_list)){
  total_cells_express=length(which(exprs(cds)[which(fData(cds)$gene_short_name==np_list[i]),]>CUTOFF))
  control_cell_express=length(which(exprs(cds_control)[which(fData(cds_control)$gene_short_name==np_list[i]),]>CUTOFF))
  restraint_cell_express=length(which(exprs(cds_restraint)[which(fData(cds_restraint)$gene_short_name==np_list[i]),]>CUTOFF))
  
  control=c(control, control_cell_express)
  restraint=c(restraint,restraint_cell_express)
  total=c(total,total_cells_express)
}

control_percent=control/length(which(cds$exp=="Control"))*100
restraint_percent=restraint/length(which(cds$exp=="Restraint"))*100
total_percent=total/length(cds$exp)*100

counts=c(control_percent,restraint_percent,total_percent)
np=rep(np_list,3)

np_ind=length(np_list)
treatment=c(rep("control",np_ind),rep("restaint",np_ind),rep("total",np_ind))
fos=rep(4,3*np_ind)

df_fos4=data.frame(np,counts,treatment,fos)

################

cds=cds_fos5
cds_control=cds[,which(cds$exp=="Control")]
cds_restraint=cds[,which(cds$exp=="Restraint")]

restraint=c()
total=c()
control=c()


for (i in 1:length(np_list)){
  total_cells_express=length(which(exprs(cds)[which(fData(cds)$gene_short_name==np_list[i]),]>CUTOFF))
  control_cell_express=length(which(exprs(cds_control)[which(fData(cds_control)$gene_short_name==np_list[i]),]>CUTOFF))
  restraint_cell_express=length(which(exprs(cds_restraint)[which(fData(cds_restraint)$gene_short_name==np_list[i]),]>CUTOFF))
  
  control=c(control, control_cell_express)
  restraint=c(restraint,restraint_cell_express)
  total=c(total,total_cells_express)
}

control_percent=control/length(which(cds$exp=="Control"))*100
restraint_percent=restraint/length(which(cds$exp=="Restraint"))*100
total_percent=total/length(cds$exp)*100

counts=c(control_percent,restraint_percent,total_percent)
np=rep(np_list,3)

np_ind=length(np_list)
treatment=c(rep("control",np_ind),rep("restaint",np_ind),rep("total",np_ind))
fos=rep(5,3*np_ind)

df_fos5=data.frame(np,counts,treatment,fos)


df=rbind(df_fos0,df_fos1,df_fos2,df_fos3,df_fos4,df_fos5)

library(patchwork)

monocle_theme_opts <- function()
{
  theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(panel.border = element_blank()) +
    theme(axis.line.x = element_line(size=0.25, color="black")) +
    theme(axis.line.y = element_line(size=0.25, color="black")) +
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
    theme(panel.background = element_rect(fill='white')) +
    theme(legend.key=element_blank())
}

df=df[df$treatment!="total",]

for (i in 1:length(np_list)){
  df_sub=df[which(df$np==np_list[i]),]
  p=ggplot(data=df_sub, aes(x=fos, y=counts, group=treatment)) +
    geom_line(aes(color=treatment),size=2)+
    #geom_point()+
    monocle_theme_opts()+
    ggtitle(c(np_list[i])) +
    xlab("Fos expression") + 
    ylab("% positive for NP")+
    scale_color_manual(values = c("royalblue","orangered"))+
    theme(plot.title = element_text(hjust = 0.5,size=24,face="bold"), 
          text=element_text(size=24),
          legend.position="none",
          axis.title.x=element_blank(),
          axis.title.y=element_blank())
  assign(paste("p",i,sep=""),p)
}

#glu
wrap_plots(p48,p45,p60)

#gaba
wrap_plots(p27,p30,p58)