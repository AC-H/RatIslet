rm(list = ls())
library('Seurat')
library('patchwork') 
library('dplyr')
library('scran') 
library('cowplot')
library(conos)
PATH <- '/home/ach/Doc/Proyecto_INGAP/sc-Vivoli/'
setwd(paste0(PATH,'filtered_feature_bc_matrix_O_Rep1'))
Oleate_0<-ReadMtx(mtx = "matrix.mtx.gz", features = "features.tsv.gz", cells = "barcodes.tsv.gz")
Oleate_1 <- CreateSeuratObject(counts = Oleate_0, project ="Oleate_Rep1")
Oleate_1[["percent.mt"]] <- PercentageFeatureSet(Oleate_1, pattern = "^Mt-")

setwd(paste0(PATH,'SeuratAnalysis_Oleate'))
# Visualize QC metrics as a violin plot
pdf("01_nGene_percent.mito_nUMI_VlnPlot_Rep1.pdf", height=10, width=20)
VlnPlot(Oleate_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
graphics.off()
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
pdf("02_nGene_percent.mito_nUMI_GenePlot_Rep1.pdf", height=10, width=20)
FeatureScatter(Oleate_1, feature1 = "nCount_RNA", feature2 = "percent.mt")
graphics.off()
# Filtering step. Exclude cells with min & max nFeatures, percent.mito and nCount_RNA
Oleate_R1 <- subset(Oleate_1, subset = nFeature_RNA > 2000 & nFeature_RNA < 5000 & nCount_RNA > 5000 & nCount_RNA < 50000 & percent.mt<10)

PATH <- '/home/ach/Doc/Proyecto_INGAP/sc-Vivoli/'
setwd(paste0(PATH,'filtered_feature_bc_matrix_O_Rep2'))
Oleate_0<-ReadMtx(mtx = "matrix.mtx.gz", features = "features.tsv.gz", cells = "barcodes.tsv.gz")
Oleate_2 <- CreateSeuratObject(counts = Oleate_0, project ="Oleate_Rep2")
Oleate_2[["percent.mt"]] <- PercentageFeatureSet(Oleate_2, pattern = "^Mt-")

setwd(paste0(PATH,'SeuratAnalysis_Oleate'))
# Visualize QC metrics as a violin plot
pdf("01_nGene_percent.mito_nUMI_VlnPlot_Rep2.pdf", height=10, width=20)
VlnPlot(Oleate_2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
graphics.off()
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
pdf("02_nGene_percent.mito_nUMI_GenePlot_Rep2.pdf", height=10, width=20)
FeatureScatter(Oleate_2, feature1 = "nCount_RNA", feature2 = "percent.mt")
graphics.off()
# Filtering step. Exclude cells with min & max nFeatures, percent.mito and nCount_RNA
Oleate_R2 <- subset(Oleate_2, subset = nFeature_RNA > 2000 & nFeature_RNA < 5000 & nCount_RNA > 10000 & nCount_RNA < 40000 & percent.mt>1 & percent.mt<10)

PATH <- '/home/ach/Doc/Proyecto_INGAP/sc-Vivoli/'
setwd(paste0(PATH,'filtered_feature_bc_matrix_O_Rep4'))
Oleate_0<-ReadMtx(mtx = "matrix.mtx.gz", features = "features.tsv.gz", cells = "barcodes.tsv.gz")
Oleate_4 <- CreateSeuratObject(counts = Oleate_0, project ="Oleate_Rep4")
Oleate_4[["percent.mt"]] <- PercentageFeatureSet(Oleate_4, pattern = "^Mt-")

setwd(paste0(PATH,'SeuratAnalysis_Oleate'))
# Visualize QC metrics as a violin plot
pdf("01_nGene_percent.mito_nUMI_VlnPlot_Rep4.pdf", height=10, width=20)
VlnPlot(Oleate_4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
graphics.off()
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
pdf("02_nGene_percent.mito_nUMI_GenePlot_Rep2.pdf", height=10, width=20)
FeatureScatter(Oleate_4, feature1 = "nCount_RNA", feature2 = "percent.mt")
graphics.off()
# Filtering step. Exclude cells with min & max nFeatures, percent.mito and nCount_RNA
Oleate_R4 <- subset(Oleate_4, subset = nFeature_RNA > 2000 & nFeature_RNA < 5000 & nCount_RNA > 5000 & nCount_RNA < 50000 & percent.mt> 1 & percent.mt<8)

PATH <- '/media/grupo_srs/Datos/Ana/Proyecto_INGAP/sc-Vivoli/Vivoli_Oleate_Rep3_chrMT/outs/'
setwd(paste0(PATH,'filtered_feature_bc_matrix_O_Rep3'))
Oleate_0<-ReadMtx(mtx = "matrix.mtx.gz", features = "features.tsv.gz", cells = "barcodes.tsv.gz")
Oleate_3 <- CreateSeuratObject(counts = Oleate_0, project ="Oleate_Rep3")
Oleate_3[["percent.mt"]] <- PercentageFeatureSet(Oleate_3, pattern = "^Mt-")

setwd('/media/grupo_srs/Datos/Ana/Proyecto_INGAP/sc-Vivoli/SeuratAnalysis_Vivoli')
# Visualize QC metrics as a violin plot
pdf("01_nGene_percent.mito_nUMI_VlnPlot_Rep3.pdf", height=10, width=20)
VlnPlot(Oleate_3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
graphics.off()
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
pdf("02_nGene_percent.mito_nUMI_GenePlot_Rep3.pdf", height=10, width=20)
FeatureScatter(Oleate_3, feature1 = "nCount_RNA", feature2 = "percent.mt")
graphics.off()
# Filtering step. Exclude cells with min & max nFeatures, percent.mito and nCount_RNA
Oleate_R3 <- subset(Oleate_3, subset = nFeature_RNA > 2000 & nFeature_RNA < 5000 & nCount_RNA > 5000 & nCount_RNA < 40000 & percent.mt>1 & percent.mt<10)

markers_islet=read.csv('/home/ach/Doc/Proyecto_Procr/1_Mouse_Islets/5_Integrated/1_Seurat/marcadores_RevFinal.csv')
markers<-markers_islet
markers_islet$endothelial<-markers_islet$B.cell.T.cell<-markers_islet$macrophage<-NULL  
markers_mouse=(unique(unlist(markers, use.names = FALSE)))
markers_short<-(unique(unlist(markers_islet, use.names = FALSE)))
markers_mouse<-c(markers_mouse, 'Mki67','Cdk1','Pcna','MSTRG.16201')


# Identification of highly variable features (feature selection)
#-----------------------Wide variable genes to interesting ones-----------------------------
Oleate_R1 <- NormalizeData(Oleate_R1, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
Oleate_R1 <- FindVariableFeatures(Oleate_R1, selection.method = "vst", nfeatures = 3000)
x <- VariableFeatures(Oleate_R1)
setdiff(markers_short,x)
feat=unique(c(markers_short,x))
feat=feat[feat!=''] #Remove empty elements
VariableFeatures(Oleate_R1)<-feat
length(feat)
#---------------------Run the standard workflow for visualization and clustering-----------
Oleate_R1 <- ScaleData(Oleate_R1, verbose = FALSE)
Oleate_R1 <- RunPCA(Oleate_R1, npcs = 50, verbose = FALSE)
saveRDS(Oleate_R1, "Oleate_R1.rds")

# Identification of highly variable features (feature selection)
#-----------------------Wide variable genes to interesting ones-----------------------------
Oleate_R2 <- NormalizeData(Oleate_R2, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
Oleate_R2 <- FindVariableFeatures(Oleate_R2, selection.method = "vst", nfeatures = 3000)
x <- VariableFeatures(Oleate_R2)
setdiff(markers_short,x)
feat=unique(c(markers_short,x))
feat=feat[feat!=''] #Remove empty elements
VariableFeatures(Oleate_R2)<-feat
length(feat)
#---------------------Run the standard workflow for visualization and clustering-----------
Oleate_R2 <- ScaleData(Oleate_R2, verbose = FALSE)
Oleate_R2 <- RunPCA(Oleate_R2, npcs = 50, verbose = FALSE)
saveRDS(Oleate_R2, "Oleate_R2.rds")

# Identification of highly variable features (feature selection)
#-----------------------Wide variable genes to interesting ones-----------------------------
Oleate_R3 <- NormalizeData(Oleate_R3, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
Oleate_R3 <- FindVariableFeatures(Oleate_R3, selection.method = "vst", nfeatures = 3000)
x <- VariableFeatures(Oleate_R3)
setdiff(markers_short,x)
feat=unique(c(markers_short,x))
feat=feat[feat!=''] #Remove empty elements
VariableFeatures(Oleate_R3)<-feat
length(feat)
#---------------------Run the standard workflow for visualization and clustering-----------
Oleate_R3 <- ScaleData(Oleate_R3, verbose = FALSE)
Oleate_R3 <- RunPCA(Oleate_R3, npcs = 50, verbose = FALSE)
saveRDS(Oleate_R3, "Oleate_R3.rds")

# Identification of highly variable features (feature selection)
#-----------------------Wide variable genes to interesting ones-----------------------------
Oleate_R4 <- NormalizeData(Oleate_R4, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
Oleate_R4 <- FindVariableFeatures(Oleate_R4, selection.method = "vst", nfeatures = 3000)
x <- VariableFeatures(Oleate_R4)
setdiff(markers_short,x)
feat=unique(c(markers_short,x))
feat=feat[feat!=''] #Remove empty elements
VariableFeatures(Oleate_R4)<-feat
length(feat)
#---------------------Run the standard workflow for visualization and clustering-----------
Oleate_R4 <- ScaleData(Oleate_R4, verbose = FALSE)
Oleate_R4 <- RunPCA(Oleate_R4, npcs = 50, verbose = FALSE)
saveRDS(Oleate_R4, "Oleate_R4.rds")


#-------------------------------Integration------------------------------------------------
Oleate<- list("Oleate_R1"= Oleate_R1, "Oleate_R2"= Oleate_R2, "Oleate_R3"=Oleate_R3, "Oleate_R4" = Oleate_R4)

Oleate.anchors <- FindIntegrationAnchors(object.list = Oleate, k.anchor = 5, dims = 1:20)

Oleate.integrated <- IntegrateData(anchorset = Oleate.anchors, dims = 1:20)
Oleate.integrated <- ScaleData(Oleate.integrated, verbose = FALSE)
Oleate.integrated <- RunPCA(Oleate.integrated, npcs = 20, verbose = FALSE)
#
#Embeding
Oleate.integrated <- RunUMAP(Oleate.integrated, reduction = "pca", dims = 1:30, verbose = FALSE)

Oleate.integrated=FindNeighbors(Oleate.integrated, reduction = "pca", dims = 1:20, k.param=20) 
Oleate.integrated=FindClusters(Oleate.integrated, resolution=0.4)

saveRDS(Oleate.integrated, file="Vivoli.integrated_dim1-20.rds")
#Nuevo analisis
Oleate.integrated_new<-Oleate.integrated
#Cargar analisis con nombres
#Oleate.integrated<-readRDS("/media/grupo_srs/Datos/Ana/Proyecto_Procr/1_Oleates/5_Integrated/1_Seurat/Oleate.integrated_dim1-30.rds")
table(Idents(Oleate.integrated))



#---------------------------------------------------------------------------------------------
# UMAP and Clustering
Oleate.integrated <- RunUMAP(Oleate.integrated, reduction = "pca", dims = 1:20)
Oleate.integrated <- FindNeighbors(Oleate.integrated, reduction = "pca", dims = 1:20)
Oleate.integrated <- FindClusters(Oleate.integrated, resolution = 0.4)
DefaultAssay(Oleate.integrated)<- "RNA" 
saveRDS(Oleate.integrated,file='Oleate_Oleate.rds')
Oleate_1<-readRDS('Oleate_Oleate.rds')
#-----------------------Exploring graph-----------------------------------------------------
pdf(paste0("Vivoli_Oleate_IsletMarkers.pdf"), height=6, width=16)
DimPlot(Oleate.integrated, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 4)
DotPlot(Oleate.integrated, features = markers_mouse, cols = c("white", "red"),col.min = -2.5, col.max = 2.5, dot.min = 0, cluster.idents=TRUE, dot.scale = 8) + RotatedAxis() 
graphics.off()

Cl0<-'Beta1'
Cl1<-'Beta2'
Cl2<-'Beta3'
Cl3<-'Beta4'
Cl4<-'Alpha'
Cl5<-'Beta5'
Cl6<-'Beta6'
Cl7<-'Mesenchymal'
Cl8<-'Delta'
Cl9<-'Immune'

Oleate_1_r<-RenameIdents(Oleate_1,`0` = Cl0, `1` = Cl1, `2` = Cl2, `3` = Cl3, `4` = Cl4, `5` = Cl5, `6` = Cl6, `7` = Cl7, `8` = Cl8, `9` = Cl9)
DefaultAssay(Oleate_1_r)<- "RNA" 
pdf(paste0("Vivoli_Oleate_IsletMarkers.pdf"), height=6, width=16)
DimPlot(Oleate_1_r, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 4)
DotPlot(Oleate_1_r, features = c(markers_mouse, 'MSTRG.7597','MSTRG.13389','MSTRG.4070','MSTRG.23603','MSTRG.9663','MSTRG.12772','MSTRG.7595') , cols = c("white", "red"),col.min = -2.5, col.max = 2.5, dot.min = 0, cluster.idents=TRUE, dot.scale = 8) + RotatedAxis() 
VlnPlot(Oleate_1_r, features = c('Mki67','Cdk1','Pcna','MSTRG.16201','MSTRG.7597','MSTRG.13389','MSTRG.4070','MSTRG.23603','MSTRG.9663','MSTRG.12772','MSTRG.7595'), pt.size = 0, combine = TRUE)  
graphics.off()



