

modelFile="RSImodel.txt"
expFile="symbol.txt"
riskFile="risk.TCGA.txt"
setwd("E:\\My_projects\\LUSC\\17RSI")
library(limma)
library(ggpubr)
library(ggplot2)
#读取文件,并对输入文件进行整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.5,]

#删掉正常样品
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]
data=t(data)
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
data=t(avereps(data))
#data=log2(data+1)

#读取模型文件
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

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#风险文件与模型分值文件合并
sameSample=intersect(row.names(risk), row.names(data))
risk=risk[sameSample, "Risk", drop=F]
Score=data[sameSample,]
rt=cbind(risk, Score)

#设置比较组
rt$Risk=factor(rt$Risk, levels = c("low", "high"))
type=levels(factor(rt[,"Risk"]))
comp=combn(type,2)
my_comparisons=list()
for (i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#绘制箱线图
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
  stat_compare_means(comparisons = my_comparisons) #默认为Wilxon
pdf(file = "RSI.pdf", width = 3, height = 3.5)
print(boxplot)
dev.off()

