
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("corrplot")


#���ð�
library(limma)
library(corrplot)

expFile="tcga.pyroptosisExp.txt"             #���������ļ�
geneFile="gene.txt"              #�����б��ļ�
#riskFile="risk.TCGAtrain.txt"      #�����ļ�
surFile="tcgatime.txt"
setwd("E:\\My_projects\\LUSC\\01 PyroptosisExp")     #���ù���Ŀ¼

#��ȡ��������ļ�,�������ݽ��д���
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

#ɾ��������Ʒ
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]

#��ȡĿ����������
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),]
data=t(data)
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
data=avereps(data)
data=log2(data+1)

#��ȡ�����ļ�
#risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#���ݺϲ�
#sameSample=intersect(row.names(data), row.names(risk))
#data=cbind(risk[sameSample,"riskScore",drop=F], data[sameSample,,drop=F])

#��ȡ�����ļ�
#sur=read.table(surFile, header = T, sep="\t", check.names = F, row.names = 1)
#���ݺϲ�
#sameSample=intersect(row.names(data), row.names(sur))
#data=cbind(sur[sameSample, "futime", drop=F], data[sameSample,,drop=F])


#����Ծ���
M=cor(data)
res1=cor.mtest(data, conf.level=0.95)

#���������ͼ��
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
