#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggplot2")


#引用包
library(limma)
library(ggplot2)
expFile="symbol.txt"        #表达输入文件
cluFile="cluster.txt"       #分型结果文件
logFCfilter=1           #logFC过滤条件(logFC=1: Fold change=2     logFC=0.585: Fold change=1.5)
fdrFilter=0.05              #fdr过滤条件
setwd("E:\\My_projects\\LUSC\\04ClusterdiffExp")     #设置工作目录

#读取表达数据文件,并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

##去除正常样品
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]
colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", colnames(data))

#读取分型结果文件
Type=read.table(cluFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(colnames(data), row.names(Type))
data=data[,sameSample,drop=F]
Type=Type[sameSample,,drop=F]

#提取不同分型的样品
c1=Type[Type$Cluster=="C1",,drop=F]
c2=Type[Type$Cluster=="C2",,drop=F]
dataC1=data[,row.names(c1)]
dataC2=data[,row.names(c2)]
data=cbind(dataC1, dataC2)
data=data[rowMeans(data)>0.5,]
conNum=ncol(dataC1)
treatNum=ncol(dataC2)
Type=c(rep(1,conNum), rep(2,treatNum))

#差异分析
outTab=data.frame()
for(i in row.names(data)){
	rt=data.frame(expression=data[i,], Type=Type)
	wilcoxTest=wilcox.test(expression ~ Type, data=rt)
	pvalue=wilcoxTest$p.value
	conGeneMeans=mean(data[i,1:conNum])
	treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
	logFC=log2(treatGeneMeans)-log2(conGeneMeans)
	conMed=median(data[i,1:conNum])
	treatMed=median(data[i,(conNum+1):ncol(data)])
	diffMed=treatMed-conMed
	if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
		outTab=rbind(outTab,cbind(gene=i,Mean1=conGeneMeans,Mean2=treatGeneMeans,logFC=logFC,pValue=pvalue))
	}
}
pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)), method="fdr")
outTab=cbind(outTab, fdr=fdr)

#输出差异表格
outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter),]
write.table(outDiff, file="diff.txt", sep="\t", row.names=F, quote=F)

#输出差异基因的表达文件
diffExp=rbind(ID=colnames(data[as.vector(outDiff[,1]),]), data[as.vector(outDiff[,1]),])
write.table(diffExp, file="diffGeneExp.txt", sep="\t", col.names=F, quote=F)
