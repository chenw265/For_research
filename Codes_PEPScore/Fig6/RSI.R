

modelFile="RSImodel.txt"
expFile="symbol.txt"
riskFile="risk.TCGA.txt"
setwd("E:\\My_projects\\LUSC\\17RSI")
library(limma)
library(ggpubr)
library(ggplot2)
#��ȡ�ļ�,���������ļ���������
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.5,]

#ɾ��������Ʒ
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]
data=t(data)
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
data=t(avereps(data))
#data=log2(data+1)

#��ȡģ���ļ�
model=read.table(modelFile, header = T, sep = "\t", check.names = F)
modelgene=model$gene
samegene=intersect(modelgene,rownames(data))
data=data[samegene,]
data=t(data)
data=as.data.frame(data)
RSIscore=0
for (i in 1:ncol(data)) {
  RSIscore=data[,i]*model$coef[i]+RSIscore
}
data=cbind(data, RSIscore)

#��ȡ�����ļ�
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#�����ļ���ģ�ͷ�ֵ�ļ��ϲ�
sameSample=intersect(row.names(risk), row.names(data))
risk=risk[sameSample, "Risk", drop=F]
Score=data[sameSample,]
rt=cbind(risk, Score)

#���ñȽ���
rt$Risk=factor(rt$Risk, levels = c("low", "high"))
type=levels(factor(rt[,"Risk"]))
comp=combn(type,2)
my_comparisons=list()
for (i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#��������ͼ
boxplot=ggboxplot(rt, x="Risk", y="RSIscore", fill ="Risk",
                  xlab = "Risk",
                  ylab = "RSI",
                  legend.title="Risk",
                  palette = c("#00468BFF", "#ED0000FF"),
                  size = 0.5,
                  outlier.shape = NA,
                  notch = T,
                  add = "jitter",
                  add.params = list(size=1.5, color="Risk"
                                    ))+
  stat_compare_means(comparisons = my_comparisons) #Ĭ��ΪWilxon
pdf(file = "RSI.pdf", width = 3, height = 3.5)
print(boxplot)
dev.off()
