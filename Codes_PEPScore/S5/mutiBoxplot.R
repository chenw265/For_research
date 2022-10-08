
library(limma)
library(reshape2)
library(ggpubr)

setwd("E:\\My_projects\\LUSC\\20RiskPyroptosis\\boxplot")
expFile="symbol.txt"
riskFile="risk.TCGA.txt"
geneFile="gene.txt"


#��ȡ��������ļ�
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#ɾ��������Ʒ
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]

#��ȡĿ�����
gene=read.table(geneFile, header = F, sep="\t", check.names = F)
gene=as.vector(gene[,1])
#��ȡĿ����������
sameGene=intersect(as.vector(gene), rownames(data))
geneExp=data[sameGene,]
exp=t(geneExp)
row.names(exp)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\",  row.names(exp))
exp=avereps(exp)
data=exp
data=log2(data+1)

#��ȡ�����ļ�
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#���ݺϲ�
sameSample=intersect(row.names(data), row.names(risk))
data=cbind(risk[sameSample,"Risk",drop=F], data[sameSample,,drop=F])

##################��������ͼ##################
#������ת����ggplot2�����ļ�
data=melt(data, id.vars=c("Risk"))
colnames(data)=c("Risk", "Gene", "Expression")
#��������ͼ
group=levels(factor(data$Risk))
data$Risk=factor(data$Risk, levels=c("low","high"))
bioCol=c("#00468BFF", "#ED0000FF")
bioCol=bioCol[1:length(group)]
boxplot=ggboxplot(data, x="Gene", y="Expression", fill="Risk",
                  xlab="",
                  ylab="Expression",
                  legend.title="Risk",
                  width=0.8,
                  palette=bioCol,
                  #size = 0.2,
                  #outlier.shape = NA,
                  #notch = T,
                  #add = "jitter",
                  #add.params = list(size=0.2, color="Risk")
                  )+
  rotate_x_text(90)+
  stat_compare_means(aes(group=Risk),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "")), label="p.signif")
  
#���ͼƬ
pdf(file="pyroptosisExp.diff.pdf", width=15, height=6)
print(boxplot)
dev.off()