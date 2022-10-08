#install.packages("survival")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("ConsensusClusterPlus")


#引用包
library(limma)
library(survival)
library(ConsensusClusterPlus)
expFile="tcga.pyroptosisExp.txt"     #表达数据文件
#cliFile="time.txt"                   #生存数据文件
workDir="E:\\My_projects\\LUSC\\02Cluster"      #工作目录
setwd(workDir)       #设置工作目录

#读取表达输入文件
data=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)

#删掉正常样品
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]
data=t(data)
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
data=avereps(data)
data=log2(data+1)

#读取生存数据
#cli=read.table(cliFile,sep="\t",check.names=F,header=T,row.names=1)     #读取临床文件

#数据合并并输出结果
#sameSample=intersect(row.names(data),row.names(cli))
#data=data[sameSample,]
#cli=cli[sameSample,]
#rt=cbind(cli,data)

#单因素COX分析
#sigGenes=c()
#for(i in colnames(rt)[3:ncol(rt)]){
#	cox=coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
#	coxSummary=summary(cox)
#	coxP=coxSummary$coefficients[,"Pr(>|z|)"]
#	if(coxP<0.05){ sigGenes=c(sigGenes,i) }
#}

#聚类
maxK=9
data=t(data)
#data=t(data[,sigGenes])

results=ConsensusClusterPlus(data,
              maxK=maxK,
              reps=50,
              pItem=0.8,
              pFeature=1,
              title=workDir,
              clusterAlg="pam",
              distance="euclidean",
              seed=123456,
              plot="pdf")

#输出分型结果
clusterNum=2        #分几类
Cluster=results[[clusterNum]][["consensusClass"]]
Cluster=as.data.frame(Cluster)
Cluster[,1]=paste0("C", Cluster[,1])
ClusterOut=rbind(ID=colnames(Cluster), Cluster)
write.table(ClusterOut, file="cluster.txt", sep="\t", quote=F, col.names=F)
