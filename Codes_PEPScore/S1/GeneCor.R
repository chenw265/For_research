
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("corrplot")


#引用包
library(limma)
library(corrplot)

expFile="tcga.pyroptosisExp.txt"             #表达数据文件
geneFile="gene.txt"              #基因列表文件
#riskFile="risk.TCGAtrain.txt"      #风险文件
surFile="tcgatime.txt"
setwd("E:\\My_projects\\LUSC\\01 PyroptosisExp")     #设置工作目录

#读取基因表达文件,并对数据进行处理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

#删掉正常样品
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]

#提取目标基因表达量
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),]
data=t(data)
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
data=avereps(data)
data=log2(data+1)

#读取风险文件
#risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#数据合并
#sameSample=intersect(row.names(data), row.names(risk))
#data=cbind(risk[sameSample,"riskScore",drop=F], data[sameSample,,drop=F])

#读取生存文件
#sur=read.table(surFile, header = T, sep="\t", check.names = F, row.names = 1)
#数据合并
#sameSample=intersect(row.names(data), row.names(sur))
#data=cbind(sur[sameSample, "futime", drop=F], data[sameSample,,drop=F])


#相关性矩阵
M=cor(data)
res1=cor.mtest(data, conf.level=0.95)

#绘制相关性图形
pdf(file="geneCor.pdf", width=7, height=7)
corrplot(M,
         order="original",
         method = "color",
         type = "full",
         tl.cex=0.6, pch=T,
         p.mat = res1$p,
         insig = "label_sig",
         pch.cex = 1.6,
         sig.level=0.05,
         number.cex = 1,
         col=colorRampPalette(c("#0099B4FF", "white", "#ED0000FF"))(50),
         tl.col="black")
dev.off()

