#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


library(limma)       #���ð�
setwd("E:\\My_projects\\LUSC\\01 PyroptosisExp")     #���ù���Ŀ¼

#��ȡ�����ļ����������ݽ��д���
rt=read.table("symbol.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#��ȡϸ����������ı�����
gene=read.table("gene.txt", header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), rownames(data))
geneExp=data[sameGene,]

#������
out=rbind(ID=colnames(geneExp),geneExp)
write.table(out,file="tcga.pyroptosisExp.txt",sep="\t",quote=F,col.names=F)