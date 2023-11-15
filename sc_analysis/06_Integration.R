rm(ls())
library('Seurat')
library('patchwork') 
library('dplyr')
library('scran') 
library('cowplot')
setwd('/media/grupo_srs/Datos/Ana/Proyecto_INGAP/sc-Vivoli/SeuratAnalysis_Vivoli/Integration/SANTI')
#------------------- Adult Markers ------------------------------------------------------
markers_islet=read.csv('/media/grupo_srs/Datos/Ana/Proyecto_Procr/1_Mouse_Islets/5_Integrated/1_Seurat/marcadores_RevFinal.csv')
markers<-markers_islet
markers_islet$endothelial<-markers_islet$B.cell.T.cell<-markers_islet$macrophage<-NULL  
markers_mouse=(unique(unlist(markers, use.names = FALSE)))
markers_short<-(unique(unlist(markers_islet, use.names = FALSE)))
markers_mouse<-unique(c('Ins2',markers_mouse, 'Mki67','Cdk1','Pcna','MSTRG.16201','MSTRG.7595'))

#-------------------------------Viv------------------------------------------------
Viv_R1=readRDS("/media/grupo_srs/Datos/Ana/Proyecto_INGAP/sc-Vivoli/SeuratAnalysis_Vivoli/data/Vehicle/Veh_R1.rds")
Viv_R2=readRDS("/media/grupo_srs/Datos/Ana/Proyecto_INGAP/sc-Vivoli/SeuratAnalysis_Vivoli/data/Vehicle/Veh_R2.rds")
Viv_R3=readRDS("/media/grupo_srs/Datos/Ana/Proyecto_INGAP/sc-Vivoli/SeuratAnalysis_Vivoli/data/Vehicle/Veh_R3.rds")

# para las muestras de palmitato filtrar nfeatures>2000 igual que en las otras samples
Pal_R1=readRDS("/media/grupo_srs/Datos/Ana/Proyecto_INGAP/sc-Vivoli/SeuratAnalysis_Vivoli/data/Palmitate/Pal_R1.rds")
Pal_R4=readRDS("/media/grupo_srs/Datos/Ana/Proyecto_INGAP/sc-Vivoli/SeuratAnalysis_Vivoli/data/Palmitate/Pal_R4.rds")
Pal_R3=readRDS("/media/grupo_srs/Datos/Ana/Proyecto_INGAP/sc-Vivoli/SeuratAnalysis_Vivoli/data/Palmitate/Pal_R3.rds")


Oleate_R1=readRDS("/media/grupo_srs/Datos/Ana/Proyecto_INGAP/sc-Vivoli/SeuratAnalysis_Vivoli/data/Oleate/Oleate_R1.rds")
Oleate_R2=readRDS("/media/grupo_srs/Datos/Ana/Proyecto_INGAP/sc-Vivoli/SeuratAnalysis_Vivoli/data/Oleate/Oleate_R2.rds")
Oleate_R3=readRDS("/media/grupo_srs/Datos/Ana/Proyecto_INGAP/sc-Vivoli/SeuratAnalysis_Vivoli/data/Oleate/Oleate_R3.rds")
Oleate_R4=readRDS("/media/grupo_srs/Datos/Ana/Proyecto_INGAP/sc-Vivoli/SeuratAnalysis_Vivoli/data/Oleate/Oleate_R4.rds")

Viv_R1_r<-RenameCells(Viv_R1, new.names=paste0(colnames(Viv_R1),'_1'))
Viv_R2_r<-RenameCells(Viv_R2, new.names=paste0(colnames(Viv_R2),'_2'))
Viv_R3_r<-RenameCells(Viv_R3, new.names=paste0(colnames(Viv_R3),'_3'))
Pal_R1_r<-RenameCells(Pal_R1, new.names=paste0(colnames(Pal_R1),'_4'))
Pal_R4_r<-RenameCells(Pal_R4, new.names=paste0(colnames(Pal_R4),'_5'))
Pal_R3_r<-RenameCells(Pal_R3, new.names=paste0(colnames(Pal_R3),'_6'))
Oleate_R1_r<-RenameCells(Oleate_R1, new.names=paste0(colnames(Oleate_R1),'_7'))
Oleate_R2_r<-RenameCells(Oleate_R2, new.names=paste0(colnames(Oleate_R2),'_8'))
Oleate_R3_r<-RenameCells(Oleate_R3, new.names=paste0(colnames(Oleate_R3),'_9'))
Oleate_R4_r<-RenameCells(Oleate_R4, new.names=paste0(colnames(Oleate_R4),'_10'))

All<- list("Veh_R1"= Viv_R1, "Veh_R2"= Viv_R2, "Veh_R3"=Viv_R3, "Pal_R1"= Pal_R1, "Pal_R4"= Pal_R4, "Pal_R3"=Pal_R3, "Oleate_R1"= Oleate_R1, "Oleate_R2"= Oleate_R2, "Oleate_R3"=Oleate_R3, "Oleate_R4" = Oleate_R4)
#-------------------------------------------------------------------------------------------

#All<- list("Veh_R1"= Viv_R1, "Veh_R2"= Viv_R2, "Veh_R3"=Viv_R3, "Pal_R1"= Pal_R1, "Pal_R3"=Pal_R3, "Oleate_R1"= Oleate_R1, "Oleate_R2"= Oleate_R2)
All <- lapply(X = All, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})

# features <- SelectIntegrationFeatures(object.list = ifnb.list)
# immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)

Vehicle.anchors <- FindIntegrationAnchors(object.list = All, k.anchor = 5, dims = 1:30)

All.integrated <- IntegrateData(anchorset = Vehicle.anchors, dims = 1:30)
DefaultAssay(All.integrated) <- "integrated"
All.integrated <- ScaleData(All.integrated, verbose = FALSE)
All.integrated <- RunPCA(All.integrated, npcs = 30, verbose = FALSE)
#
#---- 2D Embeding----------------------------------------------------------------------------------
All.integrated <- RunUMAP(All.integrated, reduction = "pca", dims = 1:30, verbose = FALSE)
DefaultAssay(All.integrated)<- 'integrated'
All.integrated=FindNeighbors(All.integrated, reduction = "pca", dims = 1:30, k.param=20)
All.integrated=FindClusters(All.integrated, resolution=0.4)

saveRDS(All.integrated, file="Vivoli.integrated_3000features_30npcs_res04.rds")
#----------------------- READ RDS-----------------------------------------------------------
All.integrated<-readRDS(file="Vivoli.integrated_3000features_30npcs_res04.rds")
table(Idents(All.integrated))
#-----------------------Exploring graph-----------------------------------------------------
# DefaultAssay(All.integrated)<- 'RNA'
pdf(paste0("01_Vivoli_10dataset_IsletMarkers_UMAP.pdf"), height=6, width=16)
DimPlot(All.integrated, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 4)
# DotPlot(All.integrated, features = markers_mouse, cols = c("white", "red"),col.min = -2.5, col.max = 2.5, dot.min = 0, cluster.idents=TRUE, dot.scale = 8) + RotatedAxis() 
graphics.off()


DefaultAssay(All.integrated)<- 'RNA'
pdf(paste0("02_Vivoli_10dataset_IsletMarkers_DOTPLOT.pdf"), height=6, width=16)
# DimPlot(All.integrated, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 4)
DotPlot(All.integrated, features = markers_mouse, cols = c("white", "red"),col.min = -2.5, col.max = 2.5, dot.min = 0, cluster.idents=TRUE, dot.scale = 8) + RotatedAxis() 
graphics.off()


#-------Selection of beta cells---------------------------------------------------
Beta_cells<-subset(All.integrated, idents = c('6','5','4','0','10','3','1', '12'))
Beta_id<-WhichCells(Beta_cells)
Beta_Veh1<-subset(Viv_R1_r, cells=Beta_id)
Beta_Veh2<-subset(Viv_R2_r, cells=Beta_id)
Beta_Veh3<-subset(Viv_R3_r, cells=Beta_id)
Beta_Pal1<-subset(Pal_R1_r, cells=Beta_id)
Beta_Pal4<-subset(Pal_R4_r, cells=Beta_id)
Beta_Pal3<-subset(Pal_R3_r, cells=Beta_id)
Beta_Ol1<-subset(Oleate_R1_r, cells=Beta_id)
Beta_Ol2<-subset(Oleate_R2_r, cells=Beta_id)
Beta_Ol3<-subset(Oleate_R3_r, cells=Beta_id)
Beta_Ol4<-subset(Oleate_R4_r, cells=Beta_id)

list_beta<- list("Veh_R1"= Beta_Veh1, "Veh_R2"= Beta_Veh2, "Veh_R3"=Beta_Veh3, "Pal_R1"= Beta_Pal1, "Pal_R4"= Beta_Pal4, "Pal_R3"=Beta_Pal3, "Oleate_R1"= Beta_Ol1, "Oleate_R2"= Beta_Ol2, "Oleate_R3"=Beta_Ol3, "Oleate_R4" = Beta_Ol4)

#-------------------------------------------------------------------------------------------

Betas <- lapply(X = list_beta, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})


# features <- SelectIntegrationFeatures(object.list = ifnb.list)
# immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)

Betas.anchors <- FindIntegrationAnchors(object.list = Betas, k.anchor = 5, dims = 1:30)

Betas.integrated <- IntegrateData(anchorset = Betas.anchors, dims = 1:30)
DefaultAssay(Betas.integrated) <- "integrated"
Betas.integrated <- ScaleData(Betas.integrated, verbose = FALSE)
Betas.integrated <- RunPCA(Betas.integrated, npcs = 30, verbose = FALSE)
#
#Embeding
Betas.integrated <- RunUMAP(Betas.integrated, reduction = "pca", dims = 1:30, verbose = FALSE)
DefaultAssay(Betas.integrated)<- 'integrated'
Betas.integrated=FindNeighbors(Betas.integrated, reduction = "pca", dims = 1:30, k.param=20)
Betas.integrated=FindClusters(Betas.integrated, resolution=0.1)
saveRDS(Betas.integrated, 'Betas.integrated.rds')


#*********************
Betas.integrated<-readRDS(file="Betas.integrated.rds")
table(Idents(Betas.integrated))

Betas.integrated=FindClusters(Betas.integrated, resolution=1)
table(Idents(Betas.integrated))

DefaultAssay(Betas.integrated)<- 'RNA'
pdf(paste0("BetaCells_res01.pdf"), height=6, width=16)
DimPlot(Betas.integrated, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 4)
DotPlot(Betas.integrated, features = markers_mouse , cols = c("white", "red"),col.min = -2.5, col.max = 2.5, dot.min = 0, cluster.idents=TRUE, dot.scale = 8) + RotatedAxis() 
VlnPlot(Betas.integrated, features = c('Mki67','Cdk1','Pcna','MSTRG.16201','MSTRG.23603','MSTRG.9663','MSTRG.12772','MSTRG.7595'), pt.size = 0, combine = TRUE)  
#DotPlot(Beta_cells, features = unique(top5_x_clus_Betacells$gene))+ RotatedAxis()
graphics.off()

write.csv(table(Idents(Betas.integrated)),file='n_cels_betaCl.csv')

# Celulas beta con Gcg
cl4_cells<-WhichCells(Betas.integrated, idents='4')
pdf(paste0("Integrated_cl4HighLighted.pdf"), height=6, width=16)
DimPlot(All.integrated, reduction = "umap", label = TRUE, cells.highlight= cl4_cells, cols.highlight ="darkblue", pt.size = 0.5, label.size = 4)
graphics.off()

cells_clbeta4<-subset(All.integrated, cells=cl4_cells)
DefaultAssay(cells_clbeta4)<-'integrated'
pdf("Clbeta4_nGene_percent.mito_nUMI_VlnPlot.pdf", height=10, width=20)
VlnPlot(cells_clbeta4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
graphics.off()


#Proliferativas
cl3_cells_1<-WhichCells(Betas.integrated, idents='3')
pdf(paste0("Integrated_cl3_1HighLighted.pdf"), height=6, width=16)
DimPlot(All.integrated, reduction = "umap", label = TRUE, cells.highlight= cl3_cells_1, cols.highlight ="darkblue", pt.size = 0.5, label.size = 4)
graphics.off()

cells_clbeta3<-subset(All.integrated, cells=cl3_cells_1)
DefaultAssay(cells_clbeta3)<-'integrated'
pdf("Clbeta3_nGene_percent.mito_nUMI_VlnPlot.pdf", height=10, width=20)
VlnPlot(cells_clbeta3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
graphics.off()

DefaultAssay(Betas.integrated)<-'integrated'
pdf("AllBetas_nGene_percent.mito_nUMI_VlnPlot.pdf", height=10, width=20)
VlnPlot(Betas.integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
graphics.off()


#---------------Second beta cell filtering----------------------------------------
Beta_cells_2<-subset(Betas.integrated, idents = c('16','19'), invert=TRUE)
Beta_id_2<-WhichCells(Beta_cells_2)
Beta_Veh1_2<-subset(Viv_R1_r, cells=Beta_id_2)
Beta_Veh2_2<-subset(Viv_R2_r, cells=Beta_id_2)
Beta_Veh3_2<-subset(Viv_R3_r, cells=Beta_id_2)
Beta_Pal1_2<-subset(Pal_R1_r, cells=Beta_id_2)
Beta_Pal4_2<-subset(Pal_R4_r, cells=Beta_id_2)
Beta_Pal3_2<-subset(Pal_R3_r, cells=Beta_id_2)
Beta_Ol1_2<-subset(Oleate_R1_r, cells=Beta_id_2)
Beta_Ol2_2<-subset(Oleate_R2_r, cells=Beta_id_2)
Beta_Ol3_2<-subset(Oleate_R3_r, cells=Beta_id_2)
Beta_Ol4_2<-subset(Oleate_R4_r, cells=Beta_id_2)

list_beta_2<- list("Veh_R1"= Beta_Veh1_2, "Veh_R2"= Beta_Veh2_2, "Veh_R3"=Beta_Veh3_2, "Pal_R1"= Beta_Pal1_2, "Pal_R4"= Beta_Pal4_2, "Pal_R3"=Beta_Pal3_2, "Oleate_R1"= Beta_Ol1_2, "Oleate_R2"= Beta_Ol2_2, "Oleate_R3"=Beta_Ol3_2, "Oleate_R4" = Beta_Ol4_2)

#-------------------------------------------------------------------------------------------

Betas_2 <- lapply(X = list_beta_2, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})

features <- SelectIntegrationFeatures(object.list = Betas_2)
Betas.anchors_2 <- FindIntegrationAnchors(object.list = Betas_2, k.anchor = 5, dims = 1:30, anchor.features = features)

Betas.integrated_2 <- IntegrateData(anchorset = Betas.anchors_2, dims = 1:30)


DefaultAssay(Betas.integrated_2) <- "integrated"
Betas.integrated_2 <- ScaleData(Betas.integrated_2, verbose = FALSE)
Betas.integrated_2 <- RunPCA(Betas.integrated_2, npcs = 30, verbose = FALSE)
#
#Embeding
Betas.integrated_2 <- RunUMAP(Betas.integrated_2, reduction = "pca", dims = 1:30, verbose = FALSE)
DefaultAssay(Betas.integrated_2)<- 'integrated'
Betas.integrated_2=FindNeighbors(Betas.integrated_2, reduction = "pca", dims = 1:30, k.param=20)
Betas.integrated_2=FindClusters(Betas.integrated_2, resolution=0.9)

write.csv(table(Idents(Betas.integrated_2)),file='n_cels_betaCl_2.csv')
saveRDS(Betas.integrated_2, 'Betas.integrated_2.rds')

top5_x_clus_Betacells<-betacells.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)
  
head(Betas.integrated_2@meta.data)  
DefaultAssay(Betas.integrated_2)<- 'integrated'
Idents(Betas.integrated_2)<-Betas.integrated_2$integrated_snn_res.0.9
Betas.integrated_2=FindClusters(Betas.integrated_2, resolution=0.7)

DefaultAssay(Betas.integrated_2)<- 'RNA'

pdf(paste0("BetaCells_Ins2_res09_2.pdf"), height=6, width=16)
DimPlot(Betas.integrated_2, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 4)
DotPlot(Betas.integrated_2, features = markers_mouse , cols = c("white", "red"),col.min = -2.5, col.max = 2.5, dot.min = 0, cluster.idents=TRUE, dot.scale = 8) + RotatedAxis() 
VlnPlot(Betas.integrated_2, features = c('Mki67','Cdk1','Pcna','MSTRG.16201','MSTRG.23603','MSTRG.9663','MSTRG.12772','MSTRG.7595'), pt.size = 0, combine = TRUE)  
#DotPlot(Beta_cells, features = unique(top5_x_clus_Betacells$gene))+ RotatedAxis()
graphics.off()


# Proliferativas
cl4_cells_2<-WhichCells(Betas.integrated_2, idents='4')
pdf(paste0("Integrated_cl4_2_HighLighted.pdf"), height=6, width=16)
DimPlot(All.integrated, reduction = "umap", label = TRUE, cells.highlight= cl4_cells_2, cols.highlight ="darkblue", pt.size = 0.5, label.size = 4)
graphics.off()

cells_clbeta3<-subset(All.integrated, cells=cl3_cells)
table(Idents(cells_clbeta3))

DimPlot(Betas.integrated_2, reduction = "umap", label = TRUE, cells.highlight= cl4_cells, cols.highlight ="darkblue", pt.size = 0.5, label.size = 4)



#-------------------------------------------------------------------------------------------


top5_x_clus_Betacells<-betacells.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)
 

pdf(paste0("Vivoli_BetaCellsd_NoAnno.pdf"), height=8, width=16)
DotPlot(Betas.integrated, features = Noanno , cols = c("white", "red"),col.min = -2.5, col.max = 2.5, dot.min = 0, cluster.idents=TRUE, dot.scale = 8) + RotatedAxis() 
VlnPlot(Betas.integrated, features = Noanno, pt.size = 0, combine = TRUE) 
graphics.off()

#----Verify if MSTRG.16201 is a Marker of Oleate or Palmitate treatment

Betas.integrated$group='Vehicle'
Betas.integrated$group[Betas.integrated$orig.ident=='Pal_Rep1']='Palmitate'
Betas.integrated$group[Betas.integrated$orig.ident=='Pal_Rep3']='Palmitate'
Betas.integrated$group[Betas.integrated$orig.ident=='Oleate_Rep1']='Oleate'
Betas.integrated$group[Betas.integrated$orig.ident=='Oleate_Rep2']='Oleate'

Betas.integrated$cluster<-Idents(Betas.integrated)
Betas.integrated$condition<-paste0(Betas.integrated$cluster,'_', Betas.integrated$group)
head(Betas.integrated@meta.data)
Idents(Betas.integrated)<-Betas.integrated$condition
library('openxlsx')
wb_O = createWorkbook()
for (i in 0:4){ #or however many clusters you have
try({
ident1 <- paste0(i,"_Oleate")
ident2 <- paste0(i,"_Vehicle")
condition.diffgenes <- FindMarkers(Betas.integrated, ident.1 = ident1, ident.2=ident2, min.pct=0.25, logfc.threshold=0.25)
condition.diffgenes$gene<-row.names(condition.diffgenes)
sheet_name = paste('Cl', i)
addWorksheet(wb_O, sheet_name)
writeData(wb_O, sheet_name, condition.diffgenes)
})
}
saveWorkbook(wb_O, 'Markers_OleateTreat.xlsx')

wb_P = createWorkbook()
for (i in 0:4){ #or however many clusters you have
try({
ident1 <- paste0(i,"_Palmitate")
ident2 <- paste0(i,"_Vehicle")
condition.diffgenes <- FindMarkers(Betas.integrated, ident.1 = ident1, ident.2=ident2, min.pct=0.25, logfc.threshold=0.25)
condition.diffgenes$gene<-row.names(condition.diffgenes)
sheet_name = paste('Cl', i)
addWorksheet(wb_P, sheet_name)
writeData(wb_P, sheet_name, condition.diffgenes)
})
}
saveWorkbook(wb_P, 'Markers_PalmitateTreat.xlsx')

Idents(Betas.integrated)<-Betas.integrated$cluster
pdf(paste0("VlnPlot_Mstrg16201_Mafa_Ucn3.pdf"), height=4, width=8)
VlnPlot(Betas.integrated, features = c('MSTRG.16201','Mafa','Ucn3'), pt.size = 0,split.by ='group', combine = FALSE, log = TRUE) 
graphics.off()


