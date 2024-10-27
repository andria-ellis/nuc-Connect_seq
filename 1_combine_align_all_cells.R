sessionInfo()

# The purpose of this code block is to process, integrate, and analyze single-cell gene expression data from control and restraint condition libraries, as well as a reference dataset from Chen et al., to create a unified Cell Data Set (CDS) object. After filtering low-quality cells and doublets, we use Seurat and Monocle3 to perform data integration, dimensionality reduction, and clustering. Finally, cell types are classified using a pre-trained model to support further functional interpretation of the dataset.

#Load in the 10X libraries that store the both the control and restraint conditions
#Files used in this code block are
# C1 - control library 1
# C2 - control library 2
# C4 - restraint library 1
# C5 - restraint library 2
# C6 - Restrain library 3
#Files from Chen et al. were downloaded from https://pmc.ncbi.nlm.nih.gov/articles/PMC5782816/#S23


# Output file : integrated_cds.chen.rds
#               integrated_cds_nochen.rds

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


#load in each library seperately and add experimental information

### C1
matrix_dir1 = "/net/trapnell/vol1/andria/10X/Naresh_data/C1/C1_withintrons/outs/filtered_feature_bc_matrix/"
barcode.path1 <- paste0(matrix_dir1, "barcodes.tsv")
features.path1 <- paste0(matrix_dir1, "features.tsv")
matrix.path1 <- paste0(matrix_dir1, "matrix.mtx")
mat1 <- readMM(file = matrix.path1)
feature.names1 = read.delim(features.path1, 
                            header = FALSE,
                            stringsAsFactors = FALSE)
barcode.names1 = read.delim(barcode.path1, 
                            header = FALSE,
                            stringsAsFactors = FALSE)
colnames(mat1) = barcode.names1$V1
rownames(mat1) = feature.names1$V1
barcode.names1$sample=1
barcode.names1$sample_name="Control"
barcode.names1$Doublet <- c1_doublet$V1
row.names(mat1) <- feature.names1$V1
row.names(feature.names1) <- feature.names1$V1
colnames(feature.names1) <- c("gene_id","gene_short_name","type")
row.names(barcode.names1) <- barcode.names1$V1
barcode.names1$cell_sample=barcode.names1$V1
barcode.names1$lab="buck"
barcode.names1$chen.cluster="none"

barcode.names1$PRV=mat1[which(feature.names1$gene_short_name=="PRV"),]

### C2
matrix_dir2 = "/net/trapnell/vol1/andria/10X/Naresh_data/C2/C2_withintrons/outs/filtered_feature_bc_matrix/"
barcode.path2 <- paste0(matrix_dir2, "barcodes.tsv")
features.path2 <- paste0(matrix_dir2, "features.tsv")
matrix.path2 <- paste0(matrix_dir2, "matrix.mtx")
mat2 <- readMM(file = matrix.path2)
feature.names2 = read.delim(features.path2, 
                            header = FALSE,
                            stringsAsFactors = FALSE)
barcode.names2 = read.delim(barcode.path2, 
                            header = FALSE,
                            stringsAsFactors = FALSE)
colnames(mat2) = barcode.names2$V1
rownames(mat2) = feature.names2$V1
barcode.names2$sample=2
barcode.names2$sample_name="Control"
barcode.names2$Doublet=c2_doublet$V1
row.names(mat2) <- feature.names2$V1
row.names(feature.names2) <- feature.names2$V1
colnames(feature.names2) <- c("gene_id","gene_short_name","type")
row.names(barcode.names2) <-  str_replace(barcode.names2$V1, "-1", "-2")
barcode.names2$cell_sample=str_replace(barcode.names2$V1, "-1", "-2")
barcode.names2$lab="buck"
barcode.names2$chen.cluster="none"

barcode.names2$PRV=mat2[which(feature.names2$gene_short_name=="PRV"),]

### C4 : restriant library
matrix_dir4 = "/net/trapnell/vol1/andria/10X/Naresh_data/C4/C4_withintrons/outs/filtered_feature_bc_matrix/"
barcode.path4 <- paste0(matrix_dir4, "barcodes.tsv")
features.path4 <- paste0(matrix_dir4, "features.tsv")
matrix.path4 <- paste0(matrix_dir4, "matrix.mtx")
mat4 <- readMM(file = matrix.path4)
feature.names4 = read.delim(features.path4, 
                            header = FALSE,
                            stringsAsFactors = FALSE)
barcode.names4 = read.delim(barcode.path4, 
                            header = FALSE,
                            stringsAsFactors = FALSE)
colnames(mat4) = barcode.names4$V1
rownames(mat4) = feature.names4$V1
barcode.names4$sample=4
barcode.names4$sample_name="Restraint"
barcode.names4$Doublet=c4_doublet$V1
row.names(mat4) <- feature.names4$V1
row.names(feature.names4) <- feature.names4$V1
colnames(feature.names4) <- c("gene_id","gene_short_name","type")
row.names(barcode.names4) <-  str_replace(barcode.names4$V1, "-1", "-4")
barcode.names4$cell_sample=str_replace(barcode.names4$V1, "-1", "-4")
barcode.names4$lab="buck"
barcode.names4$chen.cluster="none"

barcode.names4$PRV=mat4[which(feature.names4$gene_short_name=="PRV"),]

### C5  : restriant library
matrix_dir5 = "/net/trapnell/vol1/andria/10X/Naresh_data/C5/C5_withintrons/outs/filtered_feature_bc_matrix/"
barcode.path5 <- paste0(matrix_dir5, "barcodes.tsv")
features.path5 <- paste0(matrix_dir5, "features.tsv")
matrix.path5 <- paste0(matrix_dir5, "matrix.mtx")
mat5 <- readMM(file = matrix.path5)
feature.names5 = read.delim(features.path5, 
                            header = FALSE,
                            stringsAsFactors = FALSE)
barcode.names5 = read.delim(barcode.path5, 
                            header = FALSE,
                            stringsAsFactors = FALSE)
colnames(mat5) = barcode.names5$V1
rownames(mat5) = feature.names5$V1
barcode.names5$sample=5
barcode.names5$sample_name="Restraint"
barcode.names5$Doublet=c5_doublet$V1
row.names(mat5) <- feature.names5$V1
row.names(feature.names5) <- feature.names5$V1
colnames(feature.names5) <- c("gene_id","gene_short_name","type")
row.names(barcode.names5) <-  str_replace(barcode.names5$V1, "-1", "-5")
barcode.names5$cell_sample=str_replace(barcode.names5$V1, "-1", "-5")
barcode.names5$lab="buck"
barcode.names5$chen.cluster="none"

barcode.names5$PRV=mat5[which(feature.names5$gene_short_name=="PRV"),]

### C6 : restriant library
matrix_dir6 = "/net/trapnell/vol1/andria/10X/Naresh_data/C6/C6_withintrons/outs/filtered_feature_bc_matrix/"
barcode.path6 <- paste0(matrix_dir6, "barcodes.tsv")
features.path6 <- paste0(matrix_dir6, "features.tsv")
matrix.path6 <- paste0(matrix_dir6, "matrix.mtx")
mat6 <- readMM(file = matrix.path6)
feature.names6 = read.delim(features.path6, 
                            header = FALSE,
                            stringsAsFactors = FALSE)
barcode.names6 = read.delim(barcode.path6, 
                            header = FALSE,
                            stringsAsFactors = FALSE)
colnames(mat6) = barcode.names6$V1
rownames(mat6) = feature.names6$V1
barcode.names6$sample=6
barcode.names6$sample_name="Restraint"
barcode.names6$Doublet=c6_doublet$V1
row.names(mat6) <- feature.names6$V1
row.names(feature.names6) <- feature.names6$V1
colnames(feature.names6) <- c("gene_id","gene_short_name","type")
row.names(barcode.names6) <-  str_replace(barcode.names6$V1, "-1", "-6")
barcode.names6$cell_sample=str_replace(barcode.names6$V1, "-1", "-6")
barcode.names6$lab="buck"
barcode.names6$chen.cluster="none"
barcode.names6$PRV=mat6[which(feature.names6$gene_short_name=="PRV"),]


#Load in Chen data. expression data is stored in a matrix, and the cluster "cell type" assignments are stored in a table

matrix_dir = "/net/trapnell/vol1/andria/Buck_Lab/"
load(paste0(matrix_dir, "GSE87544_1443737Cells.Expresssion.Matrix.log_tpm+1_.renamed.RData"))
chen.clusters=read.csv("/net/trapnell/vol1/andria/Buck_Lab/GSE87544_1443737Cells.SVM.cluster.identity.renamed.csv")

dim(Expresssion_Matrix_unfiltered)
my_list=match(rownames(Expresssion_Matrix_unfiltered),feature.names1$gene_short_name)
mynums=1:23284

#filter gene lists so that both libraries are comparing the same genes
df1=data.frame(mynums,my_list)
colnames(df1)=c("C_gene_loc","B_gene_loc")
remove=which(is.na(df1$B_gene_loc))
df_gene_indices=df1[-remove,]
chen_ind_list=df_gene_indices$C_gene_loc
buck_ind_list=df_gene_indices$B_gene_loc

#filter all experimental libraries
mat1=mat1[buck_ind_list,]
mat2=mat2[buck_ind_list,]
mat4=mat4[buck_ind_list,]
mat5=mat5[buck_ind_list,]
mat6=mat6[buck_ind_list,]

feature.names1=feature.names1[buck_ind_list,]
feature.names2=feature.names2[buck_ind_list,]
feature.names4=feature.names4[buck_ind_list,]
feature.names5=feature.names5[buck_ind_list,]
feature.names6=feature.names6[buck_ind_list,]

#format Chen data to monocle cds format
mat.chen=Expresssion_Matrix_unfiltered[chen_ind_list,]
feature.names.chen.temp=rownames(Expresssion_Matrix_unfiltered)[chen_ind_list]
barcode.names.chen.temp=colnames(Expresssion_Matrix_unfiltered)

feature.names.chen=data.frame(feature.names1$gene_id,feature.names.chen.temp,feature.names1$gene_id)

colnames(feature.names.chen)=c("gene_id","gene_short_name","type")


condition=t(data.frame(str_split(barcode.names.chen.temp, "_"))[2,])
cell=t(data.frame(str_split(barcode.names.chen.temp, "_"))[1,])

barcode.names.chen=data.frame(condition,cell)
colnames(barcode.names.chen)=c("sample_name","V1")
barcode.names.chen$cell_sample=barcode.names.chen.temp
barcode.names.chen$lab="chen"
barcode.names.chen$PRV=0
barcode.names.chen$sample=7
barcode.names.chen$Doublet=0
barcode.names.chen$chen.cluster=chen.clusters$SVM_clusterID

#combine expression matrices
expression_matrix=cbind(mat1,mat2,mat4,mat5,mat6,as.matrix(mat.chen))


#combine cell info files
barcode.names=rbind(barcode.names1,barcode.names2,barcode.names4,barcode.names5,barcode.names6,
                    barcode.names.chen)

#combine gene info files
feature.names=join(feature.names1,feature.names2,by = "gene_id",type="left",match="all")
feature.names=join(feature.names,feature.names4,by = "gene_id",type="left",match="all")
feature.names=join(feature.names,feature.names5,by = "gene_id",type="left",match="all")
feature.names=join(feature.names,feature.names6,by = "gene_id",type="left",match="all")
feature.names=join(feature.names,feature.names.chen,by = "gene_id",type="left",match="all")


row.names(expression_matrix) <- feature.names$gene_id
row.names(feature.names) <- feature.names$gene_id

row.names(barcode.names) <- barcode.names$cell_sample
colnames(expression_matrix) <- barcode.names$cell_sample

as.matrix(expression_matrix)

#create merged cds
cds <- new_cell_data_set(as.matrix(expression_matrix),
                         cell_metadata = barcode.names,
                         gene_metadata = feature.names)

num_genes_expressed=colSums(expression_matrix != 0)
num_cells_expressed=rowSums(expression_matrix != 0)

cell_umi=colSums(expression_matrix)
pData(cds)$UMI=cell_umi
pData(cds)$num_genes=num_genes_expressed
fData(cds)$num_cells=num_cells_expressed

#filter low expressing cells (less than 500 genes)

cds_filter_low_cells=cds[,-(which(pData(cds)$num_genes<500))]

#filter doublets

cds_filter_doublets=cds_filter_low_cells[,-(which(pData(cds_filter_low_cells)$Doublet==1))]

cds=cds_filter_doublets

cds_buck=cds[,which(cds$lab=="buck")]
cds_chen=cds[,which(cds$lab=="chen")]

coldata_df=colData(cds) %>% as.data.frame()

count_mat = assay(cds)

coldata_df$Cell = cds$cell_sample

#create a seurat object to integrate the two libraries

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
                                  dims = 1:100,
                                  verbose = T,
                                  k.anchor = 20)

cds.integrated <- IntegrateData(anchorset = anchors, 
                                dims = 1:100,
                                verbose = T)

cds.integrated <- ScaleData(cds.integrated, verbose = T)

cds.integrated <- RunPCA(cds.integrated, verbose = T, npcs = 100)

cds.integrated <- RunUMAP(cds.integrated, dims = 1:100)

cds.integrated <- FindNeighbors(cds.integrated, reduction = "pca", dims = 1:30)

cds.integrated <- FindClusters(cds.integrated, resolution = 0.5)

# Pull out the UMAP embeddings from the Seurat CDS object
seurat_umap_embeddings = 
  as.data.frame(Embeddings(cds.integrated, reduction = "umap"))

seurat_umap_embeddings$Cell = rownames(seurat_umap_embeddings)

cds$umap1=seurat_umap_embeddings$UMAP_1
cds$umap2=seurat_umap_embeddings$UMAP_2

# Load garnett mouse neuron classifier

library(org.Mm.eg.db)
library(garnett)
mm_neuro_classifier <- readRDS("/net/trapnell/vol1/andria/Buck_Lab/mmBrain_20191017.RDS")

# Classify the cells

mm_neuro_cds <- classify_cells(cds, mm_neuro_classifier,
                               db = org.Mm.eg.db,
                               cluster_extend = TRUE)

head(pData(mm_neuro_cds))

# Save the cluster extended type for each cell
df=data.frame(mm_neuro_cds$umap1,mm_neuro_cds$umap2,mm_neuro_cds$cluster_ext_type)
colnames(df)=c("umap1","umap2","type")

df=as.data.frame(colData(mm_neuro_cds))
saved_coords =as.matrix(data.frame(df$umap1,df$umap2))

reducedDims(mm_neuro_cds)[["UMAP"]] = saved_coords

cds <- cluster_cells(mm_neuro_cds)

cds$exp=cds$sample_name
cds$sample_name=NULL

#Save the final cds as the transformed, integrated libraries
saveRDS(cds,file="integrated_cds.chen.rds")