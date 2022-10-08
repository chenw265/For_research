#install.packages("corrplot")
#install.packages("circlize")


#引用包
library(corrplot)
library(circlize)
library(ggsci)
library(scales)
library(limma)       #引用包
checkpoint=c("BTLA","VTCN1","CD276","CD274","CTLA4","IDO1","LAG3","PDCD1","TIGIT","HAVCR2","VSIR")
CF=c("CCL5","CCR2","CCR5","CXCL9","CXCR3")
#tmbFile="TMB.txt"                       #肿瘤突变文件
riskFile="risk.TCGA.txt"             #风险文件
#immuneFile="MCPcounter.result.txt"      #免疫细胞浸润结果
expFile="symbol.txt"
setwd("E:\\My_projects\\LUSC\\12ImmuneRelated\\06ImmuneGeneCor")     #设置工作目录

#读取免疫细胞浸润结果文件
#data=read.table(immuneFile, header=T, sep="\t", check.names=F, row.names=1)
#data=t(data)

#读取基因表达文件
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt, check.names=F)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#删掉正常样品
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]

#提取目标基因表达量
gene=checkpoint
sameGene=intersect(as.vector(gene), rownames(data))
geneExp=data[sameGene,]
exp=t(geneExp)
row.names(exp)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\",  row.names(exp))
exp=avereps(exp)
data=exp
data=log2(data+1)

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#数据合并
sameSample=intersect(row.names(data), row.names(risk))
data=cbind(risk[sameSample,"riskScore",drop=F], data[sameSample,,drop=F])

#合并肿瘤突变符合文件
#TMB=read.table(tmbFile, header=T, sep="\t", check.names=F, row.names=1)
#sameSample=intersect(row.names(TMB), row.names(data))
#rt=cbind(TMB[sameSample,,drop=F], data[sameSample,,drop=F])
rt=data
#计算相关相关性矩阵
cor1=cor(rt)
res1=cor.mtest(rt, conf.level=0.95)
res1[["p"]]

#设置图形颜色
col = c(rgb(1,0.6,0.45,seq(1,0,length=32)),rgb(0.4,0.6,0.95,seq(0,1,length=32)))
cor1[cor1==1]=0
c1 = ifelse(c(cor1)>=0,rgb(1,0.6,0.45,abs(cor1)),rgb(0.4,0.6,0.95,abs(cor1)))
col1 = matrix(c1,nc=ncol(rt))

#c("#6BAED6FF","#F39B7FFF","#C7C7C7FF")蓝橙灰
pal=c("#6BAED6FF","#6BAED6FF","#6BAED6FF","#C7C7C7FF","#6BAED6FF","#6BAED6FF","#6BAED6FF","#6BAED6FF","#6BAED6FF","#6BAED6FF","#6BAED6FF","#F05C3BFF")
#library(RColorBrewer)
#colCount=ncol(rt)
#getPalette=colorRampPalette(brewer.pal(12,"Set3"),alpha=0.5)

#绘制圈图
pdf(file="CHECKPOINTcircos.pdf", width=7, height=7)
par(mar=c(2,2,2,4))
circos.par(gap.degree=c(3,rep(2, nrow(cor1)-1)), start.degree = 180)
chordDiagram(cor1, grid.col=pal, col=col1, transparency = 0.5, symmetric = T)
#chordDiagram(cor1, grid.col=rainbow(ncol(rt),alpha = 0.5), col=col1, transparency = 0.5, symmetric = T)
par(xpd=T)
colorlegend(col, vertical = T,labels=c(1,0,-1),xlim=c(1.1,1.3),ylim=c(-0.4,0.4)) #绘制图例
dev.off()
circos.clear()



################################## 趋化因子圈图#########################################
#读取基因表达文件
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt, check.names=F)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#删掉正常样品
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]
#提取目标基因表达量
gene=CF
sameGene=intersect(as.vector(gene), rownames(data))
geneExp=data[sameGene,]
exp=t(geneExp)
row.names(exp)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\",  row.names(exp))
exp=avereps(exp)
data=exp
data=log2(data+1)
#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#数据合并
sameSample=intersect(row.names(data), row.names(risk))
data=cbind(risk[sameSample,"riskScore",drop=F], data[sameSample,,drop=F])

#合并肿瘤突变符合文件
#TMB=read.table(tmbFile, header=T, sep="\t", check.names=F, row.names=1)
#sameSample=intersect(row.names(TMB), row.names(data))
#rt=cbind(TMB[sameSample,,drop=F], data[sameSample,,drop=F])
rt=data
#计算相关相关性矩阵
cor1=cor(rt)
res1=cor.mtest(rt, conf.level=0.95)
res1[["p"]]
#设置图形颜色
col = c(rgb(1,0.6,0.45,seq(1,0,length=32)),rgb(0.4,0.6,0.95,seq(0,1,length=32)))
cor1[cor1==1]=0
c1 = ifelse(c(cor1)>=0,rgb(1,0.6,0.45,abs(cor1)),rgb(0.4,0.6,0.95,abs(cor1)))
col1 = matrix(c1,nc=ncol(rt))

#c("#6BAED6FF","#F39B7FFF","#C7C7C7FF")蓝橙灰
pal=c("#6BAED6FF","#6BAED6FF","#6BAED6FF","#6BAED6FF","#6BAED6FF","#F05C3BFF")
#library(RColorBrewer)
#colCount=ncol(rt)
#getPalette=colorRampPalette(brewer.pal(12,"Set3"),alpha=0.5)

#绘制圈图
pdf(file="CFcircos.pdf", width=7, height=7)
par(mar=c(2,2,2,4))
circos.par(gap.degree=c(3,rep(2, nrow(cor1)-1)), start.degree = 180)
chordDiagram(cor1, grid.col=pal, col=col1, transparency = 0.5, symmetric = T)
#chordDiagram(cor1, grid.col=rainbow(ncol(rt),alpha = 0.5), col=col1, transparency = 0.5, symmetric = T)
par(xpd=T)
colorlegend(col, vertical = T,labels=c(1,0,-1),xlim=c(1.1,1.3),ylim=c(-0.4,0.4)) #绘制图例
dev.off()
circos.clear()
