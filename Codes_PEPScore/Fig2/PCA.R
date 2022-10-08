
#install.packages("Rtsne")
#install.packages("ggplot2")


#引用包
library(Rtsne)
library(ggplot2)
library(limma)
library(survival)
library(ConsensusClusterPlus)
expFile="tcga.pyroptosisExp.txt"     #表达数据文件
clusterFile="cluster.txt"                   #生存数据文件
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

cl=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)

#数据合并
sameSample=intersect(row.names(data),row.names(cl))
data=data[sameSample,]
cl=cl[sameSample,]
rt=cbind(cl,data)

write.table(rt, file="clusterExpression.txt", sep="\t", quote=F, col.names=T)

bioPCA=function(inputFile=null, pcaFile=null, tsneFile=null){
	#读取输入文件,提取数据
	rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
	data=rt[c(2:(ncol(rt)))]
	risk=rt[,"cl"]

  #PCA分析
	data.pca=prcomp(data, scale. = TRUE)
	pcaPredict=predict(data.pca)
	PCA = data.frame(PC1 = pcaPredict[,1], PC2 = pcaPredict[,2],risk=risk)	
	#绘制PCA图
	pdf(file=pcaFile, height=4.5, width=5.5)       #保存输入出文件
	p=ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = risk)) +
	  stat_ellipse(aes(fill=risk), alpha=1/5)+
		scale_colour_manual(name="Cluster",  values =c("#00468BFF", "#ED0000FF"))+
	    theme_bw()+
	    theme(plot.margin=unit(rep(1.5,4),'lines'))+
	    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	print(p)
	dev.off()
	

	#t-SNE分析
	tsneOut=Rtsne(data, dims=2, perplexity=10, verbose=F, max_iter=500,check_duplicates=F)
	tsne=data.frame(tSNE1 = tsneOut$Y[,1], tSNE2 = tsneOut$Y[,2],risk=risk)	
	#绘制tSNE图
	pdf(file=tsneFile, height=4.5, width=5.5)       #保存输入出文件
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
