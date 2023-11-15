rm(list = ls())
library('Seurat')
library('patchwork') 
library('dplyr')
library('scran') 
library('cowplot')
PATH <- '/media/grupo_srs/Datos/Ana/Proyecto_INGAP/sc-Vivoli/SeuratAnalysis_Vivoli'
setwd(paste0(PATH,'filtered_feature_bc_matrix_V_Rep1'))
Viv_0<-ReadMtx(mtx = "matrix.mtx.gz", features = "features.tsv.gz", cells = "barcodes.tsv.gz")
Viv_1 <- CreateSeuratObject(counts = Viv_0, project ="Veh_Rep1")
Viv_1[["percent.mt"]] <- PercentageFeatureSet(Viv_1, pattern = "^Mt-")

setwd(paste0(PATH,'SeuratAnalysis_Vehicle'))
# Visualize QC metrics as a violin plot
pdf("01_nGene_percent.mito_nUMI_VlnPlot_Rep1.pdf", height=10, width=20)
VlnPlot(Viv_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
graphics.off()
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
pdf("02_nGene_percent.mito_nUMI_GenePlot_Rep1.pdf", height=10, width=20)
FeatureScatter(Viv_1, feature1 = "nCount_RNA", feature2 = "percent.mt")
graphics.off()
# Filtering step. Exclude cells with min & max nFeatures, percent.mito and nCount_RNA
Viv_R1 <- subset(Viv_1, subset = nFeature_RNA > 2000 & nFeature_RNA < 5000 & nCount_RNA > 5000 & nCount_RNA < 50000 & percent.mt>1 & percent.mt<10)

PATH <- '/media/grupo_srs/Datos/Ana/Proyecto_INGAP/sc-Vivoli/SeuratAnalysis_Vivoli'
setwd(paste0(PATH,'filtered_feature_bc_matrix_V_Rep2'))
Viv_0<-ReadMtx(mtx = "matrix.mtx.gz", features = "features.tsv.gz", cells = "barcodes.tsv.gz")
Viv_2 <- CreateSeuratObject(counts = Viv_0, project ="Veh_Rep2")
Viv_2[["percent.mt"]] <- PercentageFeatureSet(Viv_2, pattern = "^Mt-")

setwd(paste0(PATH,'SeuratAnalysis_Vehicle'))
# Visualize QC metrics as a violin plot
pdf("01_nGene_percent.mito_nUMI_VlnPlot_Rep2.pdf", height=10, width=20)
VlnPlot(Viv_2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
graphics.off()
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
pdf("02_nGene_percent.mito_nUMI_GenePlot_Rep2.pdf", height=10, width=20)
FeatureScatter(Viv_2, feature1 = "nCount_RNA", feature2 = "percent.mt")
graphics.off()
# Filtering step. Exclude cells with min & max nFeatures, percent.mito and nCount_RNA
Viv_R2 <- subset(Viv_2, subset = nFeature_RNA > 2000 & nFeature_RNA < 5000 & nCount_RNA > 5000 & nCount_RNA < 40000 & percent.mt>1 & percent.mt<10)

PATH <- '/media/grupo_srs/Datos/Ana/Proyecto_INGAP/sc-Vivoli/SeuratAnalysis_Vivoli'
setwd(paste0(PATH,'filtered_feature_bc_matrix_V_Rep3'))
Viv_0<-ReadMtx(mtx = "matrix.mtx.gz", features = "features.tsv.gz", cells = "barcodes.tsv.gz")
Viv_3 <- CreateSeuratObject(counts = Viv_0, project ="Veh_Rep3")
Viv_3[["percent.mt"]] <- PercentageFeatureSet(Viv_3, pattern = "^Mt-")

setwd(paste0(PATH,'SeuratAnalysis_Vehicle'))
# Visualize QC metrics as a violin plot
pdf("01_nGene_percent.mito_nUMI_VlnPlot_Rep3.pdf", height=10, width=20)
VlnPlot(Viv_3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
graphics.off()
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
pdf("02_nGene_percent.mito_nUMI_GenePlot_Rep3.pdf", height=10, width=20)
FeatureScatter(Viv_3, feature1 = "nCount_RNA", feature2 = "percent.mt")
graphics.off()
# Filtering step. Exclude cells with min & max nFeatures, percent.mito and nCount_RNA
Viv_R3 <- subset(Viv_3, subset = nFeature_RNA > 2000 & nFeature_RNA < 5000 & nCount_RNA > 5000 & nCount_RNA < 50000 & percent.mt >2 & percent.mt<8)

#------------------- Adult Markers ------------------------------------------------------
markers_islet=read.csv('/media/grupo_srs/Datos/Ana/Proyecto_Procr/1_Mouse_Islets/5_Integrated/1_Seurat/marcadores_RevFinal.csv')
markers<-markers_islet
markers_islet$endothelial<-markers_islet$B.cell.T.cell<-markers_islet$macrophage<-NULL  
markers_mouse=(unique(unlist(markers, use.names = FALSE)))
markers_short<-(unique(unlist(markers_islet, use.names = FALSE)))
markers_mouse<-unique(c(markers_mouse, 'Mki67','Cdk1','Pcna','MSTRG.16201','MSTRG.7595'))


# Identification of highly variable features (feature selection)
#-----------------------Wide variable genes to interesting ones-----------------------------
Viv_R1 <- NormalizeData(Viv_R1, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
Viv_R1 <- FindVariableFeatures(Viv_R1, selection.method = "vst", nfeatures = 3000)
x <- VariableFeatures(Viv_R1)
setdiff(markers_short,x)
feat=unique(c(markers_short,x))
feat=feat[feat!=''] #Remove empty elements
VariableFeatures(Viv_R1)<-feat
length(feat)
#---------------------Run the standard workflow for visualization and clustering-----------
Viv_R1 <- ScaleData(Viv_R1, verbose = FALSE)
Viv_R1 <- RunPCA(Viv_R1, npcs = 50, verbose = FALSE)
saveRDS(Viv_R1, "Veh_R1.rds")

# Identification of highly variable features (feature selection)
#-----------------------Wide variable genes to interesting ones-----------------------------
Viv_R2 <- NormalizeData(Viv_R2, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
Viv_R2 <- FindVariableFeatures(Viv_R2, selection.method = "vst", nfeatures = 3000)
x <- VariableFeatures(Viv_R2)
setdiff(markers_short,x)
feat=unique(c(markers_short,x))
feat=feat[feat!=''] #Remove empty elements
VariableFeatures(Viv_R2)<-feat
length(feat)
#---------------------Run the standard workflow for visualization and clustering-----------
Viv_R2 <- ScaleData(Viv_R2, verbose = FALSE)
Viv_R2 <- RunPCA(Viv_R2, npcs = 50, verbose = FALSE)
saveRDS(Viv_R2, "Veh_R2.rds")

# Identification of highly variable features (feature selection)
#-----------------------Wide variable genes to interesting ones-----------------------------
Viv_R3 <- NormalizeData(Viv_R3, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
Viv_R3 <- FindVariableFeatures(Viv_R3, selection.method = "vst", nfeatures = 3000)
x <- VariableFeatures(Viv_R3)
setdiff(markers_short,x)
feat=unique(c(markers_short,x))
feat=feat[feat!=''] #Remove empty elements
VariableFeatures(Viv_R3)<-feat
length(feat)
#---------------------Run the standard workflow for visualization and clustering-----------
Viv_R3 <- ScaleData(Viv_R3, verbose = FALSE)
Viv_R3 <- RunPCA(Viv_R3, npcs = 50, verbose = FALSE)
saveRDS(Viv_R3, "Veh_R3.rds")

Vehicle<- list("Veh_R1"= Viv_R1, "Veh_R2"= Viv_R2, "Veh_R3"=Viv_R3)
Vehicle <- lapply(X = Vehicle, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})


Vehicle.anchors <- FindIntegrationAnchors(object.list = Vehicle, k.anchor = 5, dims = 1:20)

Vehicle.integrated <- IntegrateData(anchorset = Vehicle.anchors, dims = 1:20)
DefaultAssay(Vehicle.integrated) <- "integrated"
Vehicle.integrated <- ScaleData(Vehicle.integrated, verbose = FALSE)
Vehicle.integrated <- RunPCA(Vehicle.integrated, npcs = 20, verbose = FALSE)
#
#Embeding
Vehicle.integrated <- RunUMAP(Vehicle.integrated, reduction = "pca", dims = 1:15, verbose = FALSE)
Vehicle.integrated=FindNeighbors(Vehicle.integrated, reduction = "pca", dims = 1:15, k.param=20) 
Vehicle.integrated=FindClusters(Vehicle.integrated, resolution=0.4)

saveRDS(Vehicle.integrated, file="Vehicle.integrated.rds")
table(Idents(Vehicle.integrated))

#-------------------------------Integration------------------------------------------------
Viv_R1=readRDS("/media/grupo_srs/Datos/Ana/Proyecto_INGAP/sc-Vivoli/SeuratAnalysis_Vivoli/Veh_R1.rds")
Viv_R2=readRDS("/media/grupo_srs/Datos/Ana/Proyecto_INGAP/sc-Vivoli/SeuratAnalysis_Vivoli/Veh_R2.rds")
Viv_R3=readRDS("/media/grupo_srs/Datos/Ana/Proyecto_INGAP/sc-Vivoli/SeuratAnalysis_Vivoli/Veh_R3.rds")
Pal_R1=readRDS("/media/grupo_srs/Datos/Ana/Proyecto_INGAP/sc-Vivoli/SeuratAnalysis_Vivoli/Pal_R1.rds")
Pal_R3=readRDS("/media/grupo_srs/Datos/Ana/Proyecto_INGAP/sc-Vivoli/SeuratAnalysis_Vivoli/Pal_R3.rds")
Pal_R4=readRDS("/media/grupo_srs/Datos/Ana/Proyecto_INGAP/sc-Vivoli/Pal_R4.rds")

Oleate_R1=readRDS("/media/grupo_srs/Datos/Ana/Proyecto_INGAP/sc-Vivoli/SeuratAnalysis_Vivoli/Oleate_R1.rds")
Oleate_R2=readRDS("/media/grupo_srs/Datos/Ana/Proyecto_INGAP/sc-Vivoli/SeuratAnalysis_Vivoli/Oleate_R2.rds")
Oleate_R3=readRDS("/media/grupo_srs/Datos/Ana/Proyecto_INGAP/sc-Vivoli/Oleate_R3.rds")
Oleate_R4=readRDS("/media/grupo_srs/Datos/Ana/Proyecto_INGAP/sc-Vivoli/Oleate_R4.rds")
All<- list("Veh_R1"= Viv_R1, "Veh_R2"= Viv_R2, "Veh_R3"=Viv_R3, "Pal_R1"= Pal_R1, "Pal_R4"= Pal_R4, "Pal_R3"=Pal_R3, "Oleate_R1"= Oleate_R1, "Oleate_R2"= Oleate_R2, "Oleate_R3"=Oleate_R3, "Oleate_R4" = Oleate_R4)
#-------------------------------------------------------------------------------------------


#--------------------------------------------------------
Vehicle.integrated<-readRDS(file="Vehicle.integrated.rds")
Palmitate.integrated<-readRDS(file="Palmitate.integrated.rds")
Oleate.integrated<-readRDS(file="Oleate.integrated.rds")
#--------------------------------------------------------
All<- list("All"= Vehicle.integrated, "Palmitate"= Palmitate.integrated, "Oleate"=Oleate.integrated)
All <- lapply(X = All, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})


All.anchors <- FindIntegrationAnchors(object.list = All, k.anchor = 5, dims = 1:20)

All.integrated <- IntegrateData(anchorset = All.anchors, dims = 1:20)
DefaultAssay(All.integrated) <- "integrated"
All.integrated <- ScaleData(All.integrated, verbose = FALSE)
All.integrated <- RunPCA(All.integrated, npcs = 20, verbose = FALSE)
#
#Embeding
All.integrated <- RunUMAP(All.integrated, reduction = "pca", dims = 1:20, verbose = FALSE)
All.integrated=FindNeighbors(All.integrated, reduction = "pca", dims = 1:20, k.param=30) 
DefaultAssay(All.integrated)<- 'integrated'
All.integrated=FindClusters(All.integrated, resolution=0.2)
saveRDS(All.integrated,'All.integrated.rds')

#-----------------------Exploring graph-----------------------------------------------------
setwd('/media/grupo_srs/Datos/Ana/Proyecto_INGAP/sc-Vivoli/SeuratAnalysis_VivoliSeuratAnalysis_Integrated/')
DefaultAssay(All.integrated)<- 'RNA'
pdf(paste0("Vivoli_AllIntegrated_IsletMarkers_1.pdf"), height=6, width=16)
DimPlot(All.integrated, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 4)
DotPlot(All.integrated, features = markers_mouse, cols = c("white", "red"),col.min = -2.5, col.max = 2.5, dot.min = 0, cluster.idents=TRUE, dot.scale = 8) + RotatedAxis() 
graphics.off()


Cl0<-'Beta1'
Cl1<-'Beta2'
Cl2<-'Beta3'
Cl3<-'Beta4'
Cl4<-'Alpha'
Cl5<-'Beta5'
Cl6<-'Beta6'
Cl7<-'Mesenchymal1'
Cl8<-'Delta'
Cl9<-'Beta7'
Cl10<-'Mesenchymal2'
Cl11<-'Immune'
Cl12<-'Beta Prog'
Cl13<-'Beta-Alfa'

Viv_1_r<-RenameIdents(All.integrated,`0` = Cl0, `1` = Cl1, `2` = Cl2, `3` = Cl3, `4` = Cl4, `5` = Cl5, `6` = Cl6, `7` = Cl7, `8` = Cl8, `9` = Cl9, `10` = Cl10, `11` = Cl11, `12` = Cl12, `13` = Cl13)
DefaultAssay(Viv_1_r)<- "RNA" 
integration.markers <- FindAllMarkers(Viv_1_r, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5_x_clus<-integration.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)


Noanno<-unique(c("MSTRG.12772", "MSTRG.13913", "MSTRG.15234", "MSTRG.16201", "MSTRG.20759", "MSTRG.23603", "MSTRG.23864", "MSTRG.26754", "MSTRG.8937", "MSTRG.9074", "MSTRG.9663", "MSTRG.23864", "MSTRG.26754", "MSTRG.29762"))

pdf(paste0("Vivoli_Integrated_IsletMarkers_top5Mark.pdf"), height=8, width=16)
DimPlot(Viv_1_r, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 4)
DotPlot(Viv_1_r, features = markers_mouse , cols = c("white", "red"),col.min = -2.5, col.max = 2.5, dot.min = 0, cluster.idents=TRUE, dot.scale = 8) + RotatedAxis() 
VlnPlot(Viv_1_r, features = c('Mki67','Cdk1','Pcna', 'MSTRG.16201','MSTRG.23603','MSTRG.9663','MSTRG.12772','MSTRG.7595'), pt.size = 0, combine = TRUE) 
DotPlot(Viv_1_r, features = unique(top5_x_clus$gene))+ RotatedAxis() 
graphics.off()

pdf(paste0("Vivoli_Integrated_NoAnno.pdf"), height=8, width=16)
DotPlot(Viv_1_r, features = Noanno , cols = c("white", "red"),col.min = -2.5, col.max = 2.5, dot.min = 0, cluster.idents=TRUE, dot.scale = 8) + RotatedAxis() 
VlnPlot(Viv_1_r, features = Noanno, pt.size = 0, combine = TRUE) 
graphics.off()


#Análisis células beta

#Beta_cells<-subset(Viv_1_r, idents = c('Beta1','Beta2','Beta3','Beta4','Beta5','Beta6', 'Beta7', 'Beta Prog','Beta-Alfa'))
Beta_cells<-subset(All_i.integrated, idents = c('0','1','2','3','5','8', '11'))
Beta_cells <- NormalizeData(Beta_cells, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
Beta_cells <- FindVariableFeatures(Beta_cells, selection.method = "vst", nfeatures = 2000)
#---------------------Run the standard workflow for visualization and clustering-----------
Beta_cells <- ScaleData(Beta_cells, verbose = FALSE)
Beta_cells <- RunPCA(Beta_cells, npcs = 30, verbose = FALSE)
Beta_cells <- RunUMAP(Beta_cells, reduction = "pca", dims = 1:15, verbose = FALSE)
Beta_cells=FindNeighbors(Beta_cells, reduction = "pca", dims = 1:15, k.param=20) 
Beta_cells=FindClusters(Beta_cells, resolution=0.2)
pdf(paste0("Vivoli_BetaCells_1.pdf"), height=6, width=16)
DimPlot(Beta_cells, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 4)
DotPlot(Beta_cells, features = markers_mouse , cols = c("white", "red"),col.min = -2.5, col.max = 2.5, dot.min = 0, cluster.idents=TRUE, dot.scale = 8) + RotatedAxis() 
VlnPlot(Beta_cells, features = c('Mki67','Cdk1','Pcna','MSTRG.16201','MSTRG.23603','MSTRG.9663','MSTRG.12772','MSTRG.7595'), pt.size = 0, combine = TRUE)  
graphics.off()

pdf(paste0("Vivoli_BetaCellsd_NoAnno.pdf"), height=8, width=16)
DotPlot(Beta_cells, features = Noanno , cols = c("white", "red"),col.min = -2.5, col.max = 2.5, dot.min = 0, cluster.idents=TRUE, dot.scale = 8) + RotatedAxis() 
VlnPlot(Beta_cells, features = Noanno, pt.size = 0, combine = TRUE) 
graphics.off()


Beta_cells$group='Vehicle'
Beta_cells$group[Beta_cells$orig.ident=='Pal_Rep1']='Palmitate'
Beta_cells$group[Beta_cells$orig.ident=='Pal_Rep3']='Palmitate'
Beta_cells$group[Beta_cells$orig.ident=='Pal_Rep4']='Palmitate'
Beta_cells$group[Beta_cells$orig.ident=='Oleate_Rep1']='Oleate'
Beta_cells$group[Beta_cells$orig.ident=='Oleate_Rep2']='Oleate'
Beta_cells$group[Beta_cells$orig.ident=='Oleate_Rep3']='Oleate'
saveRDS(Beta_cells,'Beta_cells.rds')
Beta_cells<-readRDS('Beta_cells.rds')
Beta_cells$celltype.condition <- paste(Idents(Beta_cells), Beta_cells$group, sep="_")
Beta_cells$celltype <- Idents(Beta_cells)
Idents(Beta_cells) <- "celltype.condition"

for (i in 0:4){ #or however many clusters you have
try({
ident1 <- paste0(i,"_Oleate")
ident2 <- paste0(i,"_Vehicle")
condition.diffgenes <- FindMarkers(Beta_cells, ident.1 = ident1, ident.2=ident2, min.pct=0.25, logfc.threshold=0.25)
write.csv(condition.diffgenes, file=paste0("Oleate_Beta",i,".csv"))
})
}
for (i in 0:4){ #or however many clusters you have
try({
ident1 <- paste0(i,"_Palmitate")
ident2 <- paste0(i,"_Vehicle")
condition.diffgenes <- FindMarkers(Beta_cells, ident.1 = ident1, ident.2=ident2, min.pct=0.25, logfc.threshold=0.25)
write.csv(condition.diffgenes, file=paste0("Palmitate_Beta",i,".csv"))
})
}
Idents(Beta_cells)<-Beta_cells$celltype
integration.markers <- FindAllMarkers(Beta_cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

library('xlsx')
length(unique(Beta_cells$celltype))
for (i in 1:length(unique(integration.markers$cluster))){
	a<-integration.markers[integration.markers$cluster == i,]
	#write.xlsx(a, file = "Cluster_Markers.xlsx", sheetName = paste0(i), append = TRUE)	
	write.csv(x = a, file = paste0(i,"_clusterMarkers_BetaCells.csv"))
	}
wb = createWorkbook()
for (i in 0:4) {
  a<-integration.markers[integration.markers$cluster == i,]
  sheet_name = paste('cl', i)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, a)
}

saveWorkbook(wb, 'ClusterMarkers_Betacells.xlsx')


Idents(Beta_cells)<- Beta_cells$group 
condition.diffgenes <- FindMarkers(Beta_cells, ident.1 = "Palmitate", ident.2="Vehicle", min.pct=0.25, logfc.threshold=0.25)

condition.diffgenes_Ol <- FindMarkers(Beta_cells, ident.1 = "Oleate", ident.2="Vehicle", min.pct=0.25, logfc.threshold=0.25)


Idents(Beta_cells)<-Beta_cells$celltype
pdf("ViolinPlot_MSTRG.16201_Ucn3_Mafa_7datasets.pdf", height=3, width=10)
VlnPlot(Beta_cells, 
	features = c('MSTRG.16201','Ucn3','Mafa'),
	split.by = "group",
	pt.size = 0, combine = FALSE) 
graphics.off()
library(openxlsx)
for (i in 0:4)){
	a<-betacells.markers[betacells.markers$cluster == i,]
	#write.xlsx(a, file = "Cluster_Markers.xlsx", sheetName = paste0(i), append = TRUE)	
	write.csv(x = a, file = paste0(i,"_Beta_Markers.csv"))
	}
saveWorkbook(wb, 'my_workbook.xlsx')
top5_x_clus_Betacells<-betacells.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)



All_i<- list("Veh_R1"= Viv_R1, "Veh_R2"= Viv_R2, "Veh_R3"=Viv_R3, "Pal_R1"= Pal_R1, "Pal_R3"=Pal_R3, "Oleate_R1"= Oleate_R1, "Oleate_R2"= Oleate_R2)
#--------------------------------------------------------

All_i <- lapply(X = All_i, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})


All_i.anchors <- FindIntegrationAnchors(object.list = All_i, k.anchor = 5, dims = 1:20)

All_i.integrated <- IntegrateData(anchorset = All_i.anchors, dims = 1:20)
DefaultAssay(All_i.integrated) <- "integrated"
All_i.integrated <- ScaleData(All_i.integrated, verbose = FALSE)
All_i.integrated <- RunPCA(All_i.integrated, npcs = 20, verbose = FALSE)
#
#Embeding
All_i.integrated <- RunUMAP(All_i.integrated, reduction = "pca", dims = 1:20, verbose = FALSE)
All_i.integrated=FindNeighbors(All_i.integrated, reduction = "pca", dims = 1:20, k.param=30) 
DefaultAssay(All_i.integrated)<- 'integrated'
All_i.integrated=FindClusters(All_i.integrated, resolution=0.4)
saveRDS(All_i.integrated,'All.integrated_7datasets.rds')

DefaultAssay(All_i.integrated)<- 'RNA'
pdf(paste0("Vivoli_AllIntegrated_IsletMarkers_7datasets.pdf"), height=6, width=16)
DimPlot(All_i.integrated, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 4)
DotPlot(All_i.integrated, features = markers_mouse, cols = c("white", "red"),col.min = -2.5, col.max = 2.5, dot.min = 0, cluster.idents=TRUE, dot.scale = 8) + RotatedAxis() 
graphics.off()



