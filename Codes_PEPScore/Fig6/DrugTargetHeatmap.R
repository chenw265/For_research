#引用包
library(limma)
library(pheatmap)
set.seed(17980122)
expFile="symbol.txt"   #表达文件的输入
riskFile="risk.TCGA.txt"   #风险输入文件
geneFile="DrugTarget.txt"
setwd("E:\\My_projects\\LUSC\\16DrugsSensitive\\heatmap")  #工作路径设定

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

#读取基因文件
geneTable=read.table(geneFile, header=F, sep="\t", check.names=F)
gene=as.vector(geneTable[,1])

#提取Target的表达值
#samegene=intersect(gene, rownames(data))
geneExp=data[gene,,drop=F]
data=geneExp

#读取分型结果文件
Cli=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(colnames(data), row.names(Cli))
data=data[,sameSample,drop=F]
Cli=Cli[sameSample,,drop=F]

#提取不同分型的样品
highrisk=Cli[Cli$Risk=="high",,drop=F]
lowrisk=Cli[Cli$Risk=="low",,drop=F]
dataH=data[,row.names(highrisk)]
dataL=data[,row.names(lowrisk)]
data=cbind(dataL,dataH)
#data=data[rowMeans(data)>0.5,]
conNum=ncol(dataL)
treatNum=ncol(dataH)
sampleType=c(rep(1,conNum), rep(2,treatNum))

#差异分析
sigVec=c()
outTab=data.frame()
for(i in rownames(data)){
  if(sd(data[i,])<0.001){next}
  wilcoxTest=wilcox.test(data[i,] ~ sampleType)
  pvalue=wilcoxTest$p.value
  if(pvalue<0.99){
    Sig=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
    sigVec=c(sigVec, paste0(i, Sig))
    conGeneMeans=mean(data[i,1:conNum])
    treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
    logFC=log2(treatGeneMeans)-log2(conGeneMeans)
    outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
  }
}
#输出差异分析结果
write.table(outTab, file="diff.xls", sep="\t", row.names=F, quote=F)
write.table(outTab, file="diff.txt", sep="\t", row.names=F, quote=F)

#输出差异基因的表达文件
exp=data[as.vector(outTab[,1]),]
expOut=rbind(ID=colnames(exp), exp)
write.table(expOut, file="diffGeneExp.txt", sep="\t", col.names=F, quote=F)

#绘制差异基因热图
exp=log2(exp+0.1)
row.names(exp)=sigVec
Group=c(rep("Low-risk",conNum),rep("High-risk",treatNum))
names(Group)=colnames(data)
Group=as.data.frame(Group)

#rownames(geneTable)=geneTable[,1]
#geneTable=geneTable[samegene,,drop=F]
row.names(exp)=paste0(geneTable$V4,"-",sigVec)
Therapy_type=data.frame(Therapy_type=geneTable[,2],row.names = rownames(exp))

ann_colors=list()
clusterCol=c("#E64B35FF","#4DBBD5FF")
names(clusterCol)=levels(factor(Group$Group))
ann_colors[["Group"]]=clusterCol

Therapy_col=c("#9ECAE1FF","#A1D99BFF","#BCBDDCFF")
names(Therapy_col)=levels(factor(Therapy_type$Therapy_type))
ann_colors[["Therapy_type"]]=Therapy_col

pdf(file="heatmap.pdf", width=6, height=6)
pheatmap(exp, 
         annotation=Group,
         annotation_row = Therapy_type,
         annotation_colors = ann_colors,
         color = colorRampPalette(c(rep("#0099B4FF",5), "white", rep("#F05C3BFF",5)))(250),
         cluster_cols =F,
         cluster_rows =F,
         gaps_col = conNum,
         gaps_row = c(14,20),
         scale="row",
         show_colnames = F,
         show_rownames = T,
         fontsize = 8,
         fontsize_row=6,
         fontsize_col=8)
dev.off()



