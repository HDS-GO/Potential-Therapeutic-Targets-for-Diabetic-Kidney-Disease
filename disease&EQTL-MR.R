######

#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")  # 如果没有安装BiocManager包，则安装
#BiocManager::install("VariantAnnotation")  # 通过BiocManager安装VariantAnnotation包

#BiocManager::install("rtracklayer")
#install.packages("devtools")  # 安装devtools包，用于安装GitHub上的包
#devtools::install_github("mrcieu/gwasglue", force = TRUE)  # 通过devtools安装gwasglue包，force = TRUE强制重新安装

#install.packages("remotes")  # 安装remotes包，用于安装GitHub上的包
#remotes::install_github("MRCIEU/TwoSampleMR")  # 通过remotes安装TwoSampleMR包

#引用包
library(VariantAnnotation)  # 引用VariantAnnotation包，用于处理变异注释数据
library(gwasglue)  # 引用gwasglue包，用于GWAS数据的处理
library(TwoSampleMR)  # 引用TwoSampleMR包，用于进行双样本Mendelian Randomization分析

exposureFile="exposure.F.csv"         #暴露数据文件
outcomeID="ukb-b-16309"        #结局数据id(需修改)
outcomeName="Epilepsy"       #图形中展示疾病的名称
setwd("C:\\Users\\hp\\Desktop\\111")     #设置工作目录

#读取暴露数据
exposure_dat=read_exposure_data(filename=exposureFile,  # 从文件中读取暴露数据
                                sep = ",",  # 文件中的分隔符为逗号
                                snp_col = "SNP",  # SNP列名
                                beta_col = "beta.exposure",  # beta列名
                                se_col = "se.exposure",  # 标准误列名
     pval_col = "pval.exposure",  # P值列名
      effect_allele_col="effect_allele.exposure",  # 效应等位基因列名
  other_allele_col = "other_allele.exposure",  # 另一等位基因列名
      eaf_col = "eaf.exposure",  # 频率列名
     phenotype_col = "exposure",  # 表型列名
   id_col = "id.exposure",  # ID列名
     samplesize_col = "samplesize.exposure",  # 样本量列名
                                chr_col="chr.exposure", pos_col = "pos.exposure",  # 染色体和位置列名
                                clump=FALSE)  # 不进行聚类

#读取结局数据
outcomeData=extract_outcome_data(snps=exposure_dat$SNP, outcomes=outcomeID)  # 提取结局数据
write.csv(outcomeData, file="outcome.csv", row.names=F)  # 将结局数据写入CSV文件

#将暴露数据和结局数据合并
outcomeData$outcome=outcomeName  # 添加一个新列，表示疾病名称
dat=harmonise_data(exposure_dat, outcomeData)  # 将暴露数据和结局数据合并

#输出用于孟德尔随机化的工具变量
outTab=dat[dat$mr_keep=="TRUE",]  # 选择用于MR的工具变量
write.csv(outTab, file="table.SNP.csv", row.names=F)  # 将用于MR的工具变量写入CSV文件

#孟德尔随机化分析
mrResult=mr(dat)  # 进行MR分析

#对结果进行OR值的计算
mrTab=generate_odds_ratios(mrResult)  # 计算结果的奥氏比值
write.csv(mrTab, file="table.MRresult.csv", row.names=F)  # 将奥氏比值写入CSV文件

#异质性分析
heterTab=mr_heterogeneity(dat)  # 进行异质性分析
write.csv(heterTab, file="table.heterogeneity.csv", row.names=F)  # 将异质性结果写入CSV文件

#多效性检验
pleioTab=mr_pleiotropy_test(dat)  # 进行多效性检验
write.csv(pleioTab, file="table.pleiotropy.csv", row.names=F)  # 将多效性检验结果写入CSV文件

#绘制散点图
pdf(file="pic.scatter_plot.pdf", width=7.5, height=7)  # 打开PDF文件，设置宽度和高度
mr_scatter_plot(mrResult, dat)  # 绘制散点图
dev.off()  # 关闭PDF文件

#森林图
res_single=mr_singlesnp(dat)      #得到每个工具变量对结局的影响
pdf(file="pic.forest.pdf", width=7, height=5.5)  # 打开PDF文件，设置宽度和高度
mr_forest_plot(res_single)  # 绘制森林图
dev.off()  # 关闭PDF文件

#漏斗图
pdf(file="pic.funnel_plot.pdf", width=7, height=6.5)  # 打开PDF文件，设置宽度和高度
mr_funnel_plot(singlesnp_results = res_single)  # 绘制漏斗图
dev.off()  # 关闭PDF文件

#留一法敏感性分析
pdf(file="pic.leaveoneout.pdf", width=7, height=5.5)  # 打开PDF文件，设置宽度和高度
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))  # 进行留一法敏感性分析并绘制图形
dev.off()  # 关闭PDF文件


###```
