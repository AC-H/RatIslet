# PCA 3d
setwd("/media/grupo_srs/Datos/Ana/Proyecto_INGAP/De_novo_transcript_3/E_Exp_Tissue_PairedEnd_PooledReference/STRINGTIE_EB_NewReference/Tablas_resultado")
library(pca3d)
install.packages("pca3d")
install.packages("rgl")
install.packages("nat")
# Replace 0 values and calculate LOG2!!!!
matrix <- data.matrix(table1)
matrix <- matrix+1
#matrix[matrix<0.01] <- 0.01
matrix.nonzero <- as.data.frame(matrix)
matrix.log <- log(matrix.nonzero,2)
PCA_result=prcomp(matrix.log, scale=FALSE)
PCA_rotation=PCA_result[2]
PCA1=as.data.frame(PCA_rotation)
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

pdf(file="3DPCA_legend9.pdf", height=5,width=15)
par(oma=c(0,0,0,0))
par(mar=c(10,10,4,2))
PCA_plot=matrix(1:6, nrow=1, ncol=3, byrow=T)
layout(PCA_plot)
require(lattice)
data(PCA)
par.set=list(axis.line=list(col='transparent'),clip=list(panel="off"))
barplot(pca1$eig[1:3,2]/100, names.arg=c("PC1","PC2","PC3"), cex.axis=3, cex.lab=2.5, cex.names=2.2, ylim=c(0,0.5), xlab="Proportion of variance")#,main = ""
plot(pca1$ind$coord[,1], pca1$ind$coord[,2], xlab="PC1", ylab="PC2", col = PCAcolors, pch = 19,cex = 3.5, cex.axis=2.2, cex.lab=3, cex.names=2.2)
plot(pca1$ind$coord[,1], pca1$ind$coord[,3], xlab="PC1", ylab="PC3", col= PCAcolors,  pch = 19,cex = 3.5, cex.axis=2.2, cex.lab=3, cex.names=2.2)
plot(cloud(rotation.PC3 ~ rotation.PC2 * rotation.PC1, data=PCA, pch=20, col=PCAcolors, cex=2, lwd=10, par.settings=par.set ))
plot(1, col="white", xlim=c(0,1), ylim=c(-5,5), xaxt="n", yaxt="n", cex.axis=1, cex.lab=1, cex.names=1)
# ATTENTION! FOR PLOTING COLORED DOTS AT THE LEGEND IT IS NOT NECESSARY TO FOLLOW THE SAME ORDER AS SPECIFIED IN THE TABLE1 FILE, BECAUSE ORDER IS ONLY IMPORTANT THE USING THE "PLOT" COMMAND ABOVE!
op <- par(cex = 0.9)
legend("center", legend=Legend, col=LegendColors ,  cex=1, border="white", pch=19, ncol = 1,box.lty=0)
dev.off()

# Mayor1=table1[(table1$C9>1&table1$C11>1&table1$C16>1),]
