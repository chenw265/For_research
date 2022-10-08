
#install.packages("Rtsne")
#install.packages("ggplot2")


#���ð�
library(Rtsne)
library(ggplot2)
library(limma)
library(survival)
library(ConsensusClusterPlus)
expFile="tcga.pyroptosisExp.txt"     #���������ļ�
clusterFile="cluster.txt"                   #���������ļ�
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

cl=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)

#���ݺϲ�
sameSample=intersect(row.names(data),row.names(cl))
data=data[sameSample,]
cl=cl[sameSample,]
rt=cbind(cl,data)

write.table(rt, file="clusterExpression.txt", sep="\t", quote=F, col.names=T)

bioPCA=function(inputFile=null, pcaFile=null, tsneFile=null){
	#��ȡ�����ļ�,��ȡ����
	rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
	data=rt[c(2:(ncol(rt)))]
	risk=rt[,"cl"]

  #PCA����
	data.pca=prcomp(data, scale. = TRUE)
	pcaPredict=predict(data.pca)
	PCA = data.frame(PC1 = pcaPredict[,1], PC2 = pcaPredict[,2],risk=risk)	
	#����PCAͼ
	pdf(file=pcaFile, height=4.5, width=5.5)       #����������ļ�
	p=ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = risk)) +
	  stat_ellipse(aes(fill=risk), alpha=1/5)+
		scale_colour_manual(name="Cluster",  values =c("#00468BFF", "#ED0000FF"))+
	    theme_bw()+
	    theme(plot.margin=unit(rep(1.5,4),'lines'))+
	    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	print(p)
	dev.off()
	

	#t-SNE����
	tsneOut=Rtsne(data, dims=2, perplexity=10, verbose=F, max_iter=500,check_duplicates=F)
	tsne=data.frame(tSNE1 = tsneOut$Y[,1], tSNE2 = tsneOut$Y[,2],risk=risk)	
	#����tSNEͼ
	pdf(file=tsneFile, height=4.5, width=5.5)       #����������ļ�
	p=ggplot(data = tsne, aes(tSNE1, tSNE2)) + geom_point(aes(color = risk),size=2.5) + 
	  #stat_ellipse(aes(fill=risk), alpha=1/5)+
		scale_colour_manual(name="Cluster",  values =c("#00468BFF", "#ED0000FF"))+
	    theme_bw()+
	    theme(plot.margin=unit(rep(1.5,4),'lines'))+
	    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	print(p)
	dev.off()
}
bioPCA(inputFile="clusterExpression.txt", pcaFile="cluster.PCA.pdf", tsneFile="cluster.t-SNE.pdf")