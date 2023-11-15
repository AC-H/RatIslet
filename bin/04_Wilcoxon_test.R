
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
quotient=data.frame(FC9=matrix.nonzero$I9/matrix.nonzero$C9,FC11=matrix.nonzero$I11/matrix.nonzero$C11,FC16=matrix.nonzero$I16/matrix.nonzero$C16)
row.names(quotient)=row.names(matrix)
FC.log <- log(quotient,2)
FC.log=data.matrix(FC.log)

output <- matrix(ncol=1, nrow=dim(FC.log)[1])
for (i in 1:dim(FC.log)[1]){
  W_test=wilcox.test(FC.log[569,], mu = 0, alternative = "two.sided")
  output[i,]=W_test$p.value
  }
genes_Wilcox_test=data.frame(output,row.names = row.names(matrix))
min=genes_Wilcox_test[genes_Wilcox_test[,1]<0.05,]
boxplot(genes_Wilcox_test[,1])
summary(genes_Wilcox_test[,1])
#//////////////////////////////////////////////////////////////////////
control=data.frame(matrix.log$C9,matrix.log$C11,matrix.log$C16)
tratamiento=data.frame(matrix.log$I9,matrix.log$I11,matrix.log$I16)
control=data.matrix(control)
tratamiento=data.matrix(tratamiento)
res <- wilcox.test(control[569,], tratamiento[569,], paired = TRUE)
res

# DEseq
library( "DESeq2" )
library(ggplot2)
setwd("/media/grupo_srs/Datos/Ana/Proyecto_INGAP/De_novo_transcript_3/E_Exp_Tissue_PairedEnd_PooledReference/STRINGTIE_EB_NewReference/Tablas_resultado")
countData_backUp <- read.csv('Tissue_Ingap_TPM.csv', header = TRUE, sep = ",")
countData <- data.frame("gene_id"=countData_backUp$gene_id, 
                        "C9"=countData_backUp$C9,"I9"=countData_backUp$I9,
                        "C11"=countData_backUp$C11,"I11"=countData_backUp$I11,
                        "C16"=countData_backUp$C16,"I16"=countData_backUp$I16)
row.names(countData)=countData$gene_id
countData$gene_id=NULL
countData <- data.matrix(countData)
head(countData)

metaData <- read.csv('phenodata.csv', header = TRUE, sep = ",")
metaData <- data.matrix(metaData)
metaData
dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design=~treat) #, tidy = TRUE
