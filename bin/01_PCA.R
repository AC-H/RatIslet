# PCA 3d
library("FactoMineR")
library("factoextra")
setwd("/media/grupo_srs/Datos/Ana/Proyecto_INGAP/De_novo_transcript_3/E_Exp_Tissue_PairedEnd_PooledReference/STRINGTIE_EB_NewReference/Tablas_resultado")
root="/media/grupo_srs/Datos/Ana/Proyecto_INGAP/De_novo_transcript_3/E_Exp_Tissue_PairedEnd_PooledReference/STRINGTIE_EB_NewReference/Tablas_resultado"
expression_filename <- paste0(root,"/ PCA_inputData_allTissues.csv")
table <- read.csv(expression_filename, header=TRUE)
# Pre-processing data
# Replace 0 values and calculate LOG2!!!!
matrix <- data.matrix(table)
matrix <- matrix+1
#matrix[matrix<0.01] <- 0.01
matrix.nonzero <- as.data.frame(matrix)
matrix.log <- log(matrix.nonzero,2)

# 1. Perform PCA analysis with "PCA" function
pca1 = PCA(t(matrix.log), graph=F,scale.unit = TRUE)
pca1_eig=get_eigenvalue(pca1)
summary(pca1)
pca1$eig
pca1$var$coord
# # correlations between variables and PCs
pca1$ind$coord
head(pca1$x)
# 2. Perform PCA analysis with "prcomp" function
pca2=prcomp(matrix.log, scale=TRUE)
PCA_rotation=pca2[2]
summary(pca2)
  PropOfVar=summary(pca2)$importance[2,1:3]
pca2$sdev
#PCs aka Scores
x_pca2=pca2$x
head(x_pca2)
dim(x_pca2)
max(x_pca2[,1])
max(x_pca2[,2])
max(x_pca2[,3])
min(x_pca2[,1])
min(x_pca2[,2])
min(x_pca2[,3])
PCA=as.data.frame(PCA_rotation)
head(x_pca2[,1])
dimension.pca1<-dimdesc(pca2, axes=c(1,2,3,4))


PCAcolors=c('#caa0ff', '#a55af4', '#8DC3E9', '#4C88BE','#e78ea5', '#c74767',
            '#ff990e','#ff990e','#ff990e',
            "#db291b","#db291b","#db291b",
            "#000000","#000000","#000000")
Legend=c("Brain",
         "Liver",
         "Islet_24h",
         "C9","I9","C11","I11", "C16", "I16")
#,"Islet_WST_8w"
LegendColors=c("#000000",
               "#db291b",
               "#ff990e",
               '#caa0ff', '#a55af4', '#8DC3E9', '#4C88BE','#e78ea5', '#c74767')

pdf(file="3DPCA_prcomp.pdf", height=5,width=15)
par(oma=c(0,0,0,0))
par(mar=c(10,10,4,2))
PCA_plot=matrix(1:6, nrow=1, ncol=3, byrow=T)
layout(PCA_plot)
require(lattice)
data(PCA)
par.set=list(axis.line=list(col='transparent'),clip=list(panel="off"))
barplot(PropOfVar, names.arg=c("PC1","PC2","PC3"), cex.axis=3, cex.lab=2.5, cex.names=2.2, ylim=c(0,0.8), xlab="Proportion of variance")#,main = ""
plot(PCA$rotation.PC1, PCA$rotation.PC2, xlab="PC1", ylab="PC2", col = PCAcolors, pch = 19,cex = 3.5, cex.axis=2.2, cex.lab=3, cex.names=2.8)
plot(PCA$rotation.PC1, PCA$rotation.PC3, xlab="PC1", ylab="PC3", col= PCAcolors,  pch = 19,cex = 3.5, cex.axis=2.2, cex.lab=3, cex.names=2.8)
plot(cloud(rotation.PC3 ~ rotation.PC2 * rotation.PC1, data=PCA, pch=20, col=PCAcolors, cex=2, lwd=10, par.settings=par.set ))
plot(1, col="white", xlim=c(0,1), ylim=c(-5,5), xaxt="n", yaxt="n", cex.axis=1, cex.lab=1, cex.names=1)
# ATTENTION! FOR PLOTING COLORED DOTS AT THE LEGEND IT IS NOT NECESSARY TO FOLLOW THE SAME ORDER AS SPECIFIED IN THE TABLE1 FILE, BECAUSE ORDER IS ONLY IMPORTANT THE USING THE "PLOT" COMMAND ABOVE!
op <- par(cex = 0.9)
legend("center", legend=Legend, col=LegendColors ,  cex=1, border="white", pch=19, ncol = 1,box.lty=0)
dev.off()

pdf(file="3DPCA_PCA.pdf", height=5,width=15)
par(oma=c(0,0,0,0))
par(mar=c(10,10,4,2))
PCA_plot=matrix(1:6, nrow=1, ncol=3, byrow=T)
layout(PCA_plot)
require(lattice)
data(PCA)
par.set=list(axis.line=list(col='transparent'),clip=list(panel="off"))
barplot(pca1$eig[1:3,2]/100, names.arg=c("PC1","PC2","PC3"), cex.axis=3, cex.lab=2.5, cex.names=2.2, ylim=c(0,0.5), xlab="Proportion of variance")#,main = ""
plot(pca1$ind$coord[,1], pca1$ind$coord[,2], xlab="PC1", ylab="PC2", col = PCAcolors, pch = 19,cex = 3.5, cex.axis=2.2, cex.lab=3, cex.names=2.8)
plot(pca1$ind$coord[,1], pca1$ind$coord[,3], xlab="PC1", ylab="PC3", col= PCAcolors,  pch = 19,cex = 3.5, cex.axis=2.2, cex.lab=3, cex.names=2.8)
plot(cloud(pca1$ind$coord[,1] ~ pca1$ind$coord[,2] * pca1$ind$coord[,3], data=PCA, pch=20, col=PCAcolors, cex=2, lwd=10, par.settings=par.set ))
plot(1, col="white", xlim=c(0,1), ylim=c(-5,5), xaxt="n", yaxt="n", cex.axis=1, cex.lab=1, cex.names=1)
# ATTENTION! FOR PLOTING COLORED DOTS AT THE LEGEND IT IS NOT NECESSARY TO FOLLOW THE SAME ORDER AS SPECIFIED IN THE TABLE1 FILE, BECAUSE ORDER IS ONLY IMPORTANT THE USING THE "PLOT" COMMAND ABOVE!
op <- par(cex = 0.9)
legend("center", legend=Legend, col=LegendColors ,  cex=1, border="white", pch=19, ncol = 1,box.lty=0)
dev.off()

