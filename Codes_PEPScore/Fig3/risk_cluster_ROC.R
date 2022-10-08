
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("corrplot")


#���ð�
library(limma)
library(corrplot)

#expFile="tcga.pyroptosisExp.txt"             #���������ļ�
#geneFile="gene.txt"              #�����б��ļ�
cluster="cluster.txt"             #�����ļ�
riskFile="risk.TCGA.txt"      #�����ļ�
#surFile="tcgatime.txt"
setwd("E:\\My_projects\\LUSC\\21ClusterandRisk")     #���ù���Ŀ¼

#��ȡ��������ļ�,�������ݽ��д���
#rt=read.table(expFile, header=T, sep="\t", check.names=F)
#rt=as.matrix(rt)
#rownames(rt)=rt[,1]
#exp=rt[,2:ncol(rt)]
#dimnames=list(rownames(exp), colnames(exp))
#data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
#data=avereps(data)

#ɾ��������Ʒ
#group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
#group=sapply(strsplit(group,""), "[", 1)
#group=gsub("2", "1", group)
#data=data[,group==0]

#��ȡĿ����������
#geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
#data=data[as.vector(geneRT[,1]),]
#data=t(data)
#rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
#data=avereps(data)
#data=log2(data+1)

#��ȡ�ļ�
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
data=read.table(cluster, header=T, sep="\t", check.names=F, row.names=1)


#���ݺϲ�
sameSample=intersect(row.names(data), row.names(risk))
data=cbind(risk[sameSample,"Risk",drop=F], data[sameSample,,drop=F])
data$Risk[which(data$Risk=="high")] <- 1
data$Risk[which(data$Risk=="low")] <- 2
data$Cluster[which(data$Cluster=="C1")] <- 1
data$Cluster[which(data$Cluster=="C2")] <- 2

#��ȡ�����ļ�
#sur=read.table(surFile, header = T, sep="\t", check.names = F, row.names = 1)
#���ݺϲ�
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
     thresholds="best", # ����youdenָ��ѡ��roc���������ֵ��
     print.thres="best",
     print.thres.col="#ED0000FF",
     print.thres.adj=1.05,
     print.auc.col="black",
     print.auc.x=0.3,
     print.auc.y=0.1,
     auc.polygon=T,
     grid=T,
     grid.lwd=2) # ��roc��������ʾ�����ֵ��
dev.off()