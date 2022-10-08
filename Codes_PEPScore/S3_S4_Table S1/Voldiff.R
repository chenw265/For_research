#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("pheatmap")
#install.packages("ggplot2")


#引用包
library(limma)
library(ggplot2)
library(pheatmap)

logFCfilter=1         #logFC过滤条件
fdrFilter=0.05        #fdr过滤条件
expFile="symbol.txt"        #表达数据文件
cluFile="cluster.txt" 
geneFile="gene.txt"
setwd("E:\\My_projects\\LUSC\\04ClusterdiffExp")     #设置工作目录

#读取表达数据文件,并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.5,]

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
		outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
	}
}
pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)), method="fdr")
outTab=cbind(outTab, fdr=fdr)

#输出所有基因的差异情况
write.table(outTab,file="tcga.all.txt",sep="\t",row.names=F,quote=F)

#输出差异表格
outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter),]
write.table(outDiff,file="tcga.diff.txt",sep="\t",row.names=F,quote=F)

#输出差异基因的表达文件
diffExp=rbind(ID=colnames(data[as.vector(outDiff[,1]),]),data[as.vector(outDiff[,1]),])
write.table(diffExp,file="tcga.diffExp.txt",sep="\t",col.names=F,quote=F)

#绘制差异基因热图
hmExp=log2(data[as.vector(outDiff[,1]),]+0.01)
Type=c(rep("C1",conNum),rep("C2",treatNum))
names(Type)=colnames(data)
Type=as.data.frame(Type)
pdf(file="heatmap.pdf",width=10,height=8)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c(rep("#00468BFF",5), "white", rep("#ED0000FF",5)))(50),
         cluster_cols =F,
         scale="row",
         show_colnames = F,
         show_rownames = T,
         fontsize = 8,
         fontsize_row=3,
         fontsize_col=8)
dev.off()



#定义显著性
outTab$fdr=as.numeric(outTab$fdr)
outTab$logFC=as.numeric(outTab$logFC)
Significant=ifelse((outTab$fdr<fdrFilter & abs(outTab$logFC)>logFCfilter), ifelse(outTab$logFC>logFCfilter,"Up","Down"), "Not")
#绘制火山图
p = ggplot(outTab, aes(logFC, -log10(fdr)))+
    geom_point(aes(col=Significant))+
    scale_color_manual(values=c("#42B540FF", "black", "#ED0000FF"))+
    labs(title = " ")+
    theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
p=p+theme_bw()

#为火山图添加带有感兴趣基因名称的图层
gene=read.table(geneFile, header = F, sep="\t", check.names = F)
samegene=intersect(gene$V1, outDiff$gene)
outTab1=outTab
rownames(outTab1)=outTab1$gene
outTab1=outTab1[samegene,]
for_label=data.frame(outTab1,row.names = NULL)

p=p +
  geom_point(size = 3, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = gene),
    data = for_label,
    color="black"
  )

#保存为图片
pdf("vol.pdf", width=6.2, height=5.5)
print(p)
dev.off()
