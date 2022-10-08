
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("corrplot")


#引用包
library(limma)
library(corrplot)

#expFile="tcga.pyroptosisExp.txt"             #表达数据文件
#geneFile="gene.txt"              #基因列表文件
cluster="cluster.txt"             #分型文件
riskFile="risk.TCGA.txt"      #风险文件
#surFile="tcgatime.txt"
setwd("E:\\My_projects\\LUSC\\21ClusterandRisk")     #设置工作目录

#读取基因表达文件,并对数据进行处理
#rt=read.table(expFile, header=T, sep="\t", check.names=F)
#rt=as.matrix(rt)
#rownames(rt)=rt[,1]
#exp=rt[,2:ncol(rt)]
#dimnames=list(rownames(exp), colnames(exp))
#data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
#data=avereps(data)

#删掉正常样品
#group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
#group=sapply(strsplit(group,""), "[", 1)
#group=gsub("2", "1", group)
#data=data[,group==0]

#提取目标基因表达量
#geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
#data=data[as.vector(geneRT[,1]),]
#data=t(data)
#rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
#data=avereps(data)
#data=log2(data+1)

#读取文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
data=read.table(cluster, header=T, sep="\t", check.names=F, row.names=1)


#数据合并
sameSample=intersect(row.names(data), row.names(risk))
data=cbind(risk[sameSample,"Risk",drop=F], data[sameSample,,drop=F])
data$Risk[which(data$Risk=="high")] <- 1
data$Risk[which(data$Risk=="low")] <- 2
data$Cluster[which(data$Cluster=="C1")] <- 1
data$Cluster[which(data$Cluster=="C2")] <- 2

#读取生存文件
#sur=read.table(surFile, header = T, sep="\t", check.names = F, row.names = 1)
#数据合并
#sameSample=intersect(row.names(data), row.names(sur))
#data=cbind(sur[sameSample, "futime", drop=F], data[sameSample,,drop=F])

attach(data)
mytable  <-  xtabs(~Risk + Cluster)
library(gmodels)
CrossTable(Risk, Cluster)
chisq.test(mytable)
detach (data)

data=cbind(data[sameSample,,drop=F], risk[sameSample,"riskScore",drop=F])
library(pROC)
rocobj=roc(data$Cluster, data$riskScore)
pdf(file="ROC.pdf",width=6.5,height=6.5)
plot(rocobj,
     legacy.axes = TRUE,
     lwd=2,
     col="#00468BFF",
     identity.col="black",
     identity.lty=2,
     identity.lwd=2,
     print.auc=T,
     thresholds="best", # 基于youden指数选择roc曲线最佳阈值点
     print.thres="best",
     print.thres.col="#ED0000FF",
     print.thres.adj=1.05,
     print.auc.col="black",
     print.auc.x=0.3,
     print.auc.y=0.1,
     auc.polygon=T,
     grid=T,
     grid.lwd=2) # 在roc曲线上显示最佳阈值点
dev.off()
