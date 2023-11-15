rm(list = ls())
library('Seurat')
library('patchwork') 
library('dplyr')
library('scran') 
library('cowplot')

PATH <- '/home/ach/Doc/Proyecto_INGAP/sc-Vivoli/'
setwd(paste0(PATH,'filtered_feature_bc_matrix_P_Rep1'))
Pal_0<-ReadMtx(mtx = "matrix.mtx.gz", features = "features.tsv.gz", cells = "barcodes.tsv.gz")
Pal_1 <- CreateSeuratObject(counts = Pal_0, project ="Pal_Rep1")
Pal_1[["percent.mt"]] <- PercentageFeatureSet(Pal_1, pattern = "^Mt-")

setwd(paste0(PATH,'SeuratAnalysis_Palmitate'))
# Visualize QC metrics as a violin plot
pdf("01_nGene_percent.mito_nUMI_VlnPlot_Rep1.pdf", height=10, width=20)
VlnPlot(Pal_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
graphics.off()
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
pdf("02_nGene_percent.mito_nUMI_GenePlot_Rep1.pdf", height=10, width=20)
FeatureScatter(Pal_1, feature1 = "nCount_RNA", feature2 = "percent.mt")
graphics.off()
# Filtering step. Exclude cells with min & max nFeatures, percent.mito and nCount_RNA
Pal_R1 <- subset(Pal_1, subset = nFeature_RNA > 1500 & nFeature_RNA < 5000 & nCount_RNA > 5000 & nCount_RNA < 50000 & percent.mt<8)

PATH <- '/media/grupo_srs/Datos/Ana/Proyecto_INGAP/sc-Vivoli/Vivoli_Palmitate_Rep4_chrMT/outs/'
setwd(paste0(PATH,'filtered_feature_bc_matrix_P_Rep4'))
Pal_0<-ReadMtx(mtx = "matrix.mtx.gz", features = "features.tsv.gz", cells = "barcodes.tsv.gz")
Pal_4 <- CreateSeuratObject(counts = Pal_0, project ="Pal_Rep4")
Pal_4[["percent.mt"]] <- PercentageFeatureSet(Pal_4, pattern = "^Mt-")

setwd('/media/grupo_srs/Datos/Ana/Proyecto_INGAP/sc-Vivoli/SeuratAnalysis_Vivoli')
# Visualize QC metrics as a violin plot
pdf("01_nGene_percent.mito_nUMI_VlnPlot_Palm_Rep4.pdf", height=10, width=20)
VlnPlot(Pal_4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
graphics.off()
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
pdf("02_nGene_percent.mito_nUMI_GenePlot_Rep2.pdf", height=10, width=20)
FeatureScatter(Pal_4, feature1 = "nCount_RNA", feature2 = "percent.mt")
graphics.off()
# Filtering step. Exclude cells with min & max nFeatures, percent.mito and nCount_RNA
Pal_R4 <- subset(Pal_4, subset = nFeature_RNA > 1500 & nFeature_RNA < 5000 & nCount_RNA > 5000 & nCount_RNA < 40000 & percent.mt>3 & percent.mt<7)

PATH <- '/home/ach/Doc/Proyecto_INGAP/sc-Vivoli/'
setwd(paste0(PATH,'filtered_feature_bc_matrix_P_Rep3'))
Pal_0<-ReadMtx(mtx = "matrix.mtx.gz", features = "features.tsv.gz", cells = "barcodes.tsv.gz")
Pal_3 <- CreateSeuratObject(counts = Pal_0, project ="Pal_Rep3")
Pal_3[["percent.mt"]] <- PercentageFeatureSet(Pal_3, pattern = "^Mt-")

setwd(paste0(PATH,'SeuratAnalysis_Palmitate'))
# Visualize QC metrics as a violin plot
pdf("01_nGene_percent.mito_nUMI_VlnPlot_Rep3.pdf", height=10, width=20)
VlnPlot(Pal_3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
graphics.off()
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
pdf("02_nGene_percent.mito_nUMI_GenePlot_Rep3.pdf", height=10, width=20)
FeatureScatter(Pal_3, feature1 = "nCount_RNA", feature2 = "percent.mt")
graphics.off()
# Filtering step. Exclude cells with min & max nFeatures, percent.mito and nCount_RNA
Pal_R3 <- subset(Pal_3, subset = nFeature_RNA > 2000 & nFeature_RNA < 5000 & nCount_RNA > 5000 & nCount_RNA < 50000 & percent.mt>3 & percent.mt<10)

# Identification of highly variable features (feature selection)
#-----------------------Wide variable genes to interesting ones-----------------------------
Pal_R1 <- NormalizeData(Pal_R1, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
Pal_R1 <- FindVariableFeatures(Pal_R1, selection.method = "vst", nfeatures = 3000)
x <- VariableFeatures(Pal_R1)
setdiff(markers_short,x)
feat=unique(c(markers_short,x))
feat=feat[feat!=''] #Remove empty elements
VariableFeatures(Pal_R1)<-feat
length(feat)
#---------------------Run the standard workflow for visualization and clustering-----------
Pal_R1 <- ScaleData(Pal_R1, verbose = FALSE)
Pal_R1 <- RunPCA(Pal_R1, npcs = 50, verbose = FALSE)
saveRDS(Pal_R1, "Pal_R1.rds")
Pal_R1<-readRDS("Pal_R1.rds")

# Identification of highly variable features (feature selection)
#-----------------------Wide variable genes to interesting ones-----------------------------
Pal_R4 <- NormalizeData(Pal_R4, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
Pal_R4 <- FindVariableFeatures(Pal_R4, selection.method = "vst", nfeatures = 3000)
x <- VariableFeatures(Pal_R4)
setdiff(markers_short,x)
feat=unique(c(markers_short,x))
feat=feat[feat!=''] #Remove empty elements
VariableFeatures(Pal_R4)<-feat
length(feat)
#---------------------Run the standard workflow for visualization and clustering-----------
Pal_R4 <- ScaleData(Pal_R4, verbose = FALSE)
Pal_R4 <- RunPCA(Pal_R4, npcs = 50, verbose = FALSE)
saveRDS(Pal_R4, "Pal_R4.rds")

# Identification of highly variable features (feature selection)
#-----------------------Wide variable genes to interesting ones-----------------------------
Pal_R3 <- NormalizeData(Pal_R3, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
Pal_R3 <- FindVariableFeatures(Pal_R3, selection.method = "vst", nfeatures = 3000)
x <- VariableFeatures(Pal_R3)
setdiff(markers_short,x)
feat=unique(c(markers_short,x))
feat=feat[feat!=''] #Remove empty elements
VariableFeatures(Pal_R3)<-feat
length(feat)
#---------------------Run the standard workflow for visualization and clustering-----------
Pal_R3 <- ScaleData(Pal_R3, verbose = FALSE)
Pal_R3 <- RunPCA(Pal_R3, npcs = 50, verbose = FALSE)
saveRDS(Pal_R3, "Pal_R3.rds")
Pal_R3<-readRDS("Pal_R3.rds")
#-------------------------------Integration------------------------------------------------
Palmitate<- list("Pal_R1"= Pal_R1, "Pal_R4"= Pal_R4, "Pal_R3"=Pal_R3)

Palmitate.anchors <- FindIntegrationAnchors(object.list = Palmitate, k.anchor = 5, dims = 1:20)

Palmitate.integrated <- IntegrateData(anchorset = Palmitate.anchors, dims = 1:20)
Palmitate.integrated <- ScaleData(Palmitate.integrated, verbose = FALSE)
Palmitate.integrated <- RunPCA(Palmitate.integrated, npcs = 20, verbose = FALSE)
#
#Embeding
Palmitate.integrated <- RunUMAP(Palmitate.integrated, reduction = "pca", dims = 1:15, verbose = FALSE)

Palmitate.integrated=FindNeighbors(Palmitate.integrated, reduction = "pca", dims = 1:15, k.param=20) 
Palmitate.integrated=FindClusters(Palmitate.integrated, resolution=0.4)

saveRDS(Palmitate.integrated, file="Palmitate.integrated.rds")
#Nuevo analisis
Palmitate.integrated_new<-Palmitate.integrated
#Cargar analisis con nombres
#Palmitate.integrated<-readRDS("/media/grupo_srs/Datos/Ana/Proyecto_Procr/1_Palmitates/5_Integrated/1_Seurat/Palmitate.integrated_dim1-30.rds")
table(Idents(Palmitate.integrated))



#------------Relevant genes to be included in the analysis, to name populations--
#------------------- Adult Markers no endothelial, no immune, no neuronal-----------
markers_islet=read.csv('/home/ach/Doc/Proyecto_Procr/1_Mouse_Islets/5_Integrated/1_Seurat/marcadores_RevFinal.csv')
markers<-markers_islet
markers_islet$endothelial<-markers_islet$B.cell.T.cell<-markers_islet$macrophage<-NULL  
markers_mouse=(unique(unlist(markers, use.names = FALSE)))
markers_short<-(unique(unlist(markers_islet, use.names = FALSE)))
markers_mouse<-c(markers_mouse, 'Mki67','Cdk1','Pcna','MSTRG.16201')



#---------------------------------------------------------------------------------------------
# UMAP and Clustering
PPal_R4almitate.integrated <- RunUMAP(Palmitate.integrated, reduction = "pca", dims = 1:20)
Palmitate.integrated <- FindNeighbors(Palmitate.integrated, reduction = "pca", dims = 1:20)
Palmitate.integrated <- FindClusters(Palmitate.integrated, resolution = 0.4)
DefaultAssay(Palmitate.integrated)<- "RNA" 
saveRDS(Palmitate.integrated,file='Pal_Palmitate.rds')
Pal_1<-readRDS('Pal_Palmitate.rds')
#-----------------------Exploring graph-----------------------------------------------------
pdf(paste0("Vivoli_Palmitate_IsletMarkers.pdf"), height=6, width=16)
DimPlot(Palmitate.integrated, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 4)
DotPlot(Palmitate.integrated, features = markers_mouse, cols = c("white", "red"),col.min = -2.5, col.max = 2.5, dot.min = 0, cluster.idents=TRUE, dot.scale = 8) + RotatedAxis() 
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

Pal_1_r<-RenameIdents(Pal_1,`0` = Cl0, `1` = Cl1, `2` = Cl2, `3` = Cl3, `4` = Cl4, `5` = Cl5, `6` = Cl6, `7` = Cl7, `8` = Cl8, `9` = Cl9)
DefaultAssay(Pal_1_r)<- "RNA" 
pdf(paste0("Vivoli_Palmitate_IsletMarkers.pdf"), height=6, width=16)
DimPlot(Pal_1_r, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 4)
DotPlot(Pal_1_r, features = c(markers_mouse, 'MSTRG.7597','MSTRG.13389','MSTRG.4070','MSTRG.23603','MSTRG.9663','MSTRG.12772','MSTRG.7595') , cols = c("white", "red"),col.min = -2.5, col.max = 2.5, dot.min = 0, cluster.idents=TRUE, dot.scale = 8) + RotatedAxis() 
VlnPlot(Pal_1_r, features = c('Mki67','Cdk1','Pcna','MSTRG.16201','MSTRG.7597','MSTRG.13389','MSTRG.4070','MSTRG.23603','MSTRG.9663','MSTRG.12772','MSTRG.7595'), pt.size = 0, combine = TRUE)  
graphics.off()



