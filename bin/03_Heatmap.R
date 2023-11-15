library("FactoMineR")
library("factoextra")
library("Hmisc") # library for correlation function
require("gplots")
root="/media/grupo_srs/Datos/Ana/Proyecto_INGAP/De_novo_transcript_3/E_Exp_Tissue_PairedEnd_PooledReference/STRINGTIE_EB_NewReference/Tablas_resultado"
setwd(root)
expression_filename <- paste0(root,"/Liver_Ingap_Islet24h_Brain_mayor1TPMalmenos1tejido.csv")
table1 <- read.csv(expression_filename, header=TRUE)
# Asignacion de gene_id como indice
gene_names = table1[,c(1,2)]
row.names(table1) <- table1$gene_id
table1=table1[-1] # Borrar indices
table1=table1[-1] # Borrar codigos ID genes
names(table1)
# Asignacion nombres mas simples
colnames(table1)=c("C9","I9","C11","I11","C16","I16","I1","I2","I3","L1","L2","L3","B1","B2","B3")

# Replace 0 values and calculate LOG2!!!!
matrix <- data.matrix(table1)
matrix <- matrix+1
#matrix[matrix<0.01] <- 0.01
matrix.nonzero <- as.data.frame(matrix)
matrix.log <- log(matrix.nonzero,2)

matrix3.log.scaled <- t(scale(t(matrix.log), scale=TRUE))

# This avoids using NA / Inf values from expression data that is all zeros!! 
matrix3.log.scaled.2<-matrix3.log.scaled[complete.cases(matrix3.log.scaled), ]

my_palette = colorRampPalette(c("midnightblue","dodgerblue3","white","goldenrod1","darkorange2"))
library(RColorBrewer)
coul <- colorRampPalette(brewer.pal(8, "RdBu"))(25)

plot_filename=paste0(root,"/Heatmap_PruebaMaca.pdf") 
pdf(plot_filename)#, width=70, height=70, pointsize=16
pl<-heatmap.2(matrix3.log.scaled.2, col=coul, density.info=c("none"), trace=c("none"),srtCol=90)
dev.off()

#hclust coloured dendrogram
cols_branches <-rainbow(7) #c("chartreuse4","red","blue","aquamarine3")
#Customize hclust coloured dendrogram
pheno_dend <- color_branches(as.hclust(pl$rowDendrogram), k=7, col = cols_branches)
plot(pheno_dend)
pdf(paste0(root,"/Heatmap_colDend_hclust.pdf"))
heatmap.2(matrix3.log.scaled.2, col=coul, density.info=c("none"), trace=c("none"), Rowv=pheno_dend, main = "Genes")
dev.off()
write.table(matrix3.log.scaled.2, "/media/grupo_srs/Datos/Ana/Proyecto_INGAP/De_novo_transcript_3/ballgown_j2_a10_f01_Rn62020_2/Tablas_resultado/Input_data_colouredHeatmap.csv", row.names=TRUE,sep = "," ) 
# Reguperar genes del cluster 1 naranja y 4 rojo
cluster_data=table(matrixC, matrix3.log.scaled.2)
hc <- hclust(dist(matrix3.log.scaled.2), method = "complete")
plot(hc)
B=cutree(hc,4)
matrixB <- data.matrix(B)
plot(B, main="AAA")
head(B)

