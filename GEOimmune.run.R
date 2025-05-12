####

#install.packages('e1071')

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("preprocessCore")


inputFile="merge.normalize.txt"      
setwd("C:\\Users\\wode3\\Desktop\\7.4分多GEO数据库+eQTL孟德尔分析\\24疾病的免疫细胞浸润的评估")      
source("GEOimmune.CIBERSORT.R")     

outTab=CIBERSORT("ref.txt", inputFile, perm=1000, QN=TRUE)

outTab=outTab[outTab[,"P-value"]<1,]
outTab=as.matrix(outTab[,1:(ncol(outTab)-3)])
outTab=rbind(id=colnames(outTab),outTab)
write.table(outTab, file="CIBERSORT-Results.txt", sep="\t", quote=F, col.names=F)


##
