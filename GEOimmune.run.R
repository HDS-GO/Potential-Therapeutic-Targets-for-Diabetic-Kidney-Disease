####

#install.packages('e1071')

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("preprocessCore")


inputFile="merge.normalize.txt"      
setwd("C:\\Users\\wode3\\Desktop\\7.4�ֶ�GEO���ݿ�+eQTL�ϵ¶�����\\24����������ϸ�����������")      
source("GEOimmune.CIBERSORT.R")     

outTab=CIBERSORT("ref.txt", inputFile, perm=1000, QN=TRUE)

outTab=outTab[outTab[,"P-value"]<1,]
outTab=as.matrix(outTab[,1:(ncol(outTab)-3)])
outTab=rbind(id=colnames(outTab),outTab)
write.table(outTab, file="CIBERSORT-Results.txt", sep="\t", quote=F, col.names=F)


##
