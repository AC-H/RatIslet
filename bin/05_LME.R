rm(list=ls())
setwd("/media/grupo_srs/Datos/Ana/Proyecto_INGAP/De_novo_transcript_3/E_Exp_Tissue_PairedEnd_PooledReference/STRINGTIE_EB_NewReference/Tablas_resultado")

library("readxl")
library(genefilter)
library(dplyr)
library(devtools)
library(limma)
library(edgeR)
library(DESeq2)
file="/media/grupo_srs/Datos/Ana/Proyecto_INGAP/De_novo_transcript_3/E_Exp_Tissue_PairedEnd_PooledReference/STRINGTIE_EB_NewReference/Tablas_resultado/TPMdata_ingap.csv"
TPMdata_ingap=read.csv(file, row.names = 1)
head(TPMdata_ingap)

#col_order <- c("C9", "I9", "C11",
#               "I11", "C16", "I16")
#rawdata_ingap <- rawdata_ingap[, col_order]
#head(rawdata_ingap)
# Filtrar genes que no tienen expresion
TPMdata=TPMdata_ingap[rowSums(TPMdata_ingap[])>0,]
dim(TPMdata) # 24467
matrix_ingap <- data.matrix(TPMdata)
#1) raw count data --> VST transformation with
#matrix_ingap_Stab=varianceStabilizingTransformation(matrix_ingap) #using DESeq  (this was done years ago).
#head(matrix_ingap_Stab)
#2) extract log2 transformed expression data from above
matrix_ingap <- matrix_ingap+1
#matrix.nonzero <- as.data.frame(matrix)
matrix_ingap.log <- log(matrix_ingap,2)
head(matrix_ingap.log)
df.expression<-matrix_ingap.log
rawtreat_ingap=read.csv("/media/grupo_srs/Datos/Ana/Proyecto_INGAP/De_novo_transcript_3/E_Exp_Tissue_PairedEnd_PooledReference/STRINGTIE_EB_NewReference/Ingap/phenodata.csv")
#Batch <- factor(rawtreat$population)
#Treat <- factor(rawtreat$tratamiento,levels=c("ctrl","glt"))
#tiempo<-factor(rawtreat$tiempo, levels=c("3","24"))

#3) lme model is built one gene at a time, you will need numerical
# values of gene expression data of the gene, experiment factors (treatment, time etc), and donorID.
# the model for the single gene looks like + time, time=rawtreat$tiempo
library(nlme)
p<-c();
FC_LME<-c();
#i=1
for(i in 1:dim(df.expression)[1]) {
  x<-t(df.expression[i,]);
  data0<-data.frame(gene=as.numeric(x), trt=rawtreat_ingap$treat, donor=rawtreat_ingap$population);
  err<-try(m2<-lme(gene ~ trt , random=~1|donor, data=data0))
  if (inherits(err, "try-error")) {
    p[i]<-1;
  }else {
    p[i]<-summary(m2)$tTable[2,5]
    FC_LME[i]<-summary(m2)$tTable[2,1]
  }
}
sum(p<0.05, na.rm = TRUE) #2008
# el valor de la pendiente de cada modelo lineal, es de alguna manera el nuevo FC?
gene_pvalue_ingap<-data.frame(gene=rownames(TPMdata), p_val_LME=p)
head(gene_pvalue_ingap)
library(tidyr)
#gene_pvalue_ingap=gene_pvalue_ingap %>% separate(gene, c("gene_id","gene_name"), sep = "([|])")
#head(gene_pvalue_ingap)
#data_stab_ingap=data.frame(matrix_ingap_Stab)
#data_stab_ingap$id=rownames(data_stab_ingap)
#head(data_stab_ingap)
#data_stab_ingap=data_stab_ingap %>% separate(id, c("gene_id","gene_name"), sep = "([|])")
# Reconstruccion tabla entera

# Comparacion con resultados del t-test. Agarrar el flag 1 o 0 que indica si todos los genes van en el mismo sentido.
A=read_xlsx("/media/grupo_srs/Datos/Ana/Proyecto_INGAP/De_novo_transcript_3/E_Exp_Tissue_PairedEnd_PooledReference/STRINGTIE_EB_NewReference/Ingap/Tabla_datos_Ingap_tejidosPublicos_ANALYSIS_2020-12-11.xlsx", sheet="replace_0x0.01")
names(A)[names(A) == "Avg FC"] <- "avg_FC"
names(A)[names(A) == "p_val"] <- "p_val_ttest"
names(A)[names(A) == "all same direction"] <- "FlagFC"

head(A)
names(A)

Ingap_modif_table=merge(x=A, y=gene_pvalue_ingap, by.x="gene_id", by.y="gene", all.x = TRUE )

head(Ingap_modif_table)
dim(Ingap_modif_table)
boxplot(p)

write.csv(Ingap_modif_table,file="pvalue_INGAP_Stabilized_TPM.csv")

# Comparacion con resultados del t-test
santi_ttest=Ingap_modif_table %>% select(gene_id, ref_gene_name_x_x, TPM_C9,TPM_I9,TPM_C11,TPM_I11,TPM_C16,TPM_I16,FC9, FC11, FC16,avg_FC, FlagFC, p_val_ttest,p_val_LME)
santi_ttest$p_val_LME[is.na(santi_ttest$p_val_LME)] <- 1
head(santi_ttest)
names(santi_ttest)
# Cuantos genes de los 2008 LME_signif tienen flag 1?
signifLME_FlagFC=santi_ttest[(santi_ttest$FlagFC=1 & santi_ttest$p_val_LME<0.05),]#
dim(signifLME_FlagFC) #2008
signifttest_FlagFC=santi_ttest[(santi_ttest$FlagFC=1 & santi_ttest$p_val_ttest<0.05),]#
dim(signifttest_FlagFC)#1521

#Generacion de tabla comparativa entre TPM y counts con varianza estabilizada
#santi_ttest_conPvalLME=merge(x=santi_ttest, y=gene_pvalue_ingap, by.x="gene_id", by.y="gene", all.x = TRUE)
#head(santi_ttest_conPvalLME)

#losmascapitos=santi_ttest_conPvalLME[(santi_ttest_conPvalLME$p_val_ttest<0.05|santi_ttest_conPvalLME$p_val_LME<0.05),]#2135
#losmascapitos_intersec=santi_ttest_conPvalLME[(santi_ttest_conPvalLME$p_val_ttest<0.05&santi_ttest_conPvalLME$p_val_LME<0.05),]#1394
#sum(santi_ttest_conPvalLME$p_val_ttest<0.05, na.rm = TRUE) #1521
#plot(losmascapitos_intersec$p_val_ttest,losmascapitos_intersec$p_val_LME)  
# Filtrar match por minimo valor de modulacion en FC 20%
capitos20porciento=santi_ttest_conPvalLME[(santi_ttest_conPvalLME$avg_FC>1.2|santi_ttest_conPvalLME$avg_FC<0.8),]
Capitos20porciento=capitos20porciento[(capitos20porciento$p_val_LME<0.05),]
dim(Capitos20porciento)
head(Capitos20porciento)
Capitos20porciento_minExpr=Capitos20porciento[(Capitos20porciento$TPM_C9>0.5&Capitos20porciento$TPM_C11>0.5&Capitos20porciento$TPM_C16>0.5)|
                                              (Capitos20porciento$TPM_I9>0.5&Capitos20porciento$TPM_I11>0.5&Capitos20porciento$TPM_I16>0.5),]
dim(Capitos20porciento_minExpr)
write.csv(Capitos20porciento_minExpr, file="Capitos20porciento_minExpr_LME_Mayo7.csv")

Capitos20porciento_tt=capitos20porciento[(capitos20porciento$p_val_ttest<0.05),]
dim(Capitos20porciento_tt)
head(Capitos20porciento_tt)


write.csv(santi_ttest_conPvalLME, file="Ingap_TPM_conPvalueVogel.csv")
dim(santi_ttest_conPvalLME)


B=santi_ttest_conPvalLME[santi_ttest_conPvalLME$p_val_ttest<0.05,]
head(B)#1521 genes signif
dim(B)
# Merge genes con p value signif en ttest y LME

match_significativo=santi_ttest_conPvalLME[(santi_ttest_conPvalLME$p_val_LME<0.05 & santi_ttest_conPvalLME$p_val_ttest<0.05),]
dim(match_significativo)
head(match_significativo)

# Filtrar match por un minimo en TPM
match_minExpr=match[(match$TPM_C9>0.5&match$TPM_C11>0.5&match$TPM_C16>0.5)|(match$TPM_I9>0.5&match$TPM_I11>0.5&match$TPM_I16>0.5),]
dim(match_minExpr)
head(match_minExpr)
barplot(mx,main='Hours By Sprint',ylab='Hours', xlab='Sprint',beside = TRUE, 
        col=colours, ylim=c(0,max(mx)*1.3))

for(i in 1:dim(match_minExpr)[1]) {
  x<-data.frame(match_minExpr[i,]);
  gene_id=x$gene_id
  gene=x$ref_gene_name_x_x
  data0<-data.frame(TPM=c(x$TPM_C9, x$TPM_I9,x$TPM_C11,x$TPM_I11,x$TPM_C16,x$TPM_I16),
                    Counts=c(x$C9,x$I9,x$C11,x$I11,x$C16,x$I16))
  png(paste0(gene,"_",gene_id,".png"))
  barplot(t(as.matrix(data0)), col=c("darkblue","red"), main=paste0(gene,"_",gene_id),
          legend = names(data0), beside=TRUE, names.arg=c("C9", "I9", "C11", "I11","C16", "I16"))
  graphics.off()
}

# Grouped Bar Plot
counts <- table(mtcars$vs, mtcars$gear)
barplot(counts, main="Car Distribution by Gears and VS",
        xlab="Number of Gears", col=c("darkblue","red"),
        legend = rownames(counts), beside=TRUE)
data0<-data.frame(TPM_C9=x$TPM_C9,TPM_I9=x$TPM_I9,
                  TPM_C11=x$TPM_C11,TPM_I11=x$TPM_I11,
                  TPM_C16=x$TPM_C16,TPM_I16=x$TPM_I16,
                  C9=x$C9,I9=x$I9,
                  C11=x$C11,I11=x$I11,
                  C16=x$C16,I16=x$I16)
                  
# Analisis 29-04-22
TPMdata_ingap$gene<-row.names(TPMdata_ingap)
Ingap_LMEpval=merge(x=TPMdata_ingap, y=gene_pvalue_ingap, by.x='gene', by.y="gene", all.y = TRUE)
Ingap_LMEpval_signif=Ingap_LMEpval[Ingap_LMEpval$p_val_LME<0.05,]
dim(Ingap_LMEpval_signif)

Ingap_LMEpval_signif_1=Ingap_LMEpval_signif[(Ingap_LMEpval_signif$C9>0.5&Ingap_LMEpval_signif$C11>0.5&Ingap_LMEpval_signif$C16>0.5) | (Ingap_LMEpval_signif$I9>0.5&Ingap_LMEpval_signif$I11>0.5&Ingap_LMEpval_signif$I16>0.5),]
dim(Ingap_LMEpval_signif_1)

