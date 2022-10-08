#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


library(limma)       #引用包
setwd("E:\\My_projects\\LUSC\\01 PyroptosisExp")     #设置工作目录

#读取输入文件，并对数据进行处理
rt=read.table("symbol.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#提取细胞焦亡基因的表达量
gene=read.table("gene.txt", header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), rownames(data))
geneExp=data[sameGene,]

#输出结果
out=rbind(ID=colnames(geneExp),geneExp)
write.table(out,file="tcga.pyroptosisExp.txt",sep="\t",quote=F,col.names=F)
