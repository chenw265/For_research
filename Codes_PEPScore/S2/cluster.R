#install.packages("survival")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("ConsensusClusterPlus")


#���ð�
library(limma)
library(survival)
library(ConsensusClusterPlus)
expFile="tcga.pyroptosisExp.txt"     #���������ļ�
#cliFile="time.txt"                   #���������ļ�
workDir="E:\\My_projects\\LUSC\\02Cluster"      #����Ŀ¼
setwd(workDir)       #���ù���Ŀ¼

#��ȡ���������ļ�
data=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)

#ɾ��������Ʒ
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]
data=t(data)
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
data=avereps(data)
data=log2(data+1)

#��ȡ��������
#cli=read.table(cliFile,sep="\t",check.names=F,header=T,row.names=1)     #��ȡ�ٴ��ļ�

#���ݺϲ���������
#sameSample=intersect(row.names(data),row.names(cli))
#data=data[sameSample,]
#cli=cli[sameSample,]
#rt=cbind(cli,data)

#������COX����
#sigGenes=c()
#for(i in colnames(rt)[3:ncol(rt)]){
#	cox=coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
#	coxSummary=summary(cox)
#	coxP=coxSummary$coefficients[,"Pr(>|z|)"]
#	if(coxP<0.05){ sigGenes=c(sigGenes,i) }
#}

#����
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

#������ͽ��
clusterNum=2        #�ּ���
Cluster=results[[clusterNum]][["consensusClass"]]
Cluster=as.data.frame(Cluster)
Cluster[,1]=paste0("C", Cluster[,1])
ClusterOut=rbind(ID=colnames(Cluster), Cluster)
write.table(ClusterOut, file="cluster.txt", sep="\t", quote=F, col.names=F)