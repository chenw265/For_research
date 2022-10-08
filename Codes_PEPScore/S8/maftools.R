
#if (!require("BiocManager"))
#    install.packages("BiocManager")
#BiocManager::install("maftools")


library(maftools)       #引用包
setwd("E:\\My_projects\\LUSC\\11Maftools")      #设置工作目录

#读取模型文件，获取基因列表
geneRT=read.table("multiCox.txt", header=T, sep="\t", check.names=F, row.names=1)
gene=row.names(geneRT)

#读取临床数据文件
clinical=read.table("clinical.txt", header=T, sep="\t", check.names=F, row.names=1)
outTab=rbind(Tumor_Sample_Barcode=colnames(clinical), clinical)
write.table(outTab, file="ann.txt", sep="\t", quote=F, col.names=F)

#绘制瀑布图
pdf(file="oncoplot.pdf", width=8, height=7)
maf=read.maf(maf="new.maf", clinicalData="ann.txt")
oncoplot(maf=maf, clinicalFeatures=c("Age","Gender","Stage"), genes=gene)
dev.off()
