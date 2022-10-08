#install.packages("corrplot")
#install.packages("circlize")


#���ð�
library(corrplot)
library(circlize)
library(ggsci)
library(scales)
library(limma)       #���ð�
checkpoint=c("BTLA","VTCN1","CD276","CD274","CTLA4","IDO1","LAG3","PDCD1","TIGIT","HAVCR2","VSIR")
CF=c("CCL5","CCR2","CCR5","CXCL9","CXCR3")
#tmbFile="TMB.txt"                       #����ͻ���ļ�
riskFile="risk.TCGA.txt"             #�����ļ�
#immuneFile="MCPcounter.result.txt"      #����ϸ��������
expFile="symbol.txt"
setwd("E:\\My_projects\\LUSC\\12ImmuneRelated\\06ImmuneGeneCor")     #���ù���Ŀ¼

#��ȡ����ϸ���������ļ�
#data=read.table(immuneFile, header=T, sep="\t", check.names=F, row.names=1)
#data=t(data)

#��ȡ��������ļ�
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt, check.names=F)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#ɾ��������Ʒ
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]

#��ȡĿ����������
gene=checkpoint
sameGene=intersect(as.vector(gene), rownames(data))
geneExp=data[sameGene,]
exp=t(geneExp)
row.names(exp)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\",  row.names(exp))
exp=avereps(exp)
data=exp
data=log2(data+1)

#��ȡ�����ļ�
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#���ݺϲ�
sameSample=intersect(row.names(data), row.names(risk))
data=cbind(risk[sameSample,"riskScore",drop=F], data[sameSample,,drop=F])

#�ϲ�����ͻ������ļ�
#TMB=read.table(tmbFile, header=T, sep="\t", check.names=F, row.names=1)
#sameSample=intersect(row.names(TMB), row.names(data))
#rt=cbind(TMB[sameSample,,drop=F], data[sameSample,,drop=F])
rt=data
#�����������Ծ���
cor1=cor(rt)
res1=cor.mtest(rt, conf.level=0.95)
res1[["p"]]

#����ͼ����ɫ
col = c(rgb(1,0.6,0.45,seq(1,0,length=32)),rgb(0.4,0.6,0.95,seq(0,1,length=32)))
cor1[cor1==1]=0
c1 = ifelse(c(cor1)>=0,rgb(1,0.6,0.45,abs(cor1)),rgb(0.4,0.6,0.95,abs(cor1)))
col1 = matrix(c1,nc=ncol(rt))

#c("#6BAED6FF","#F39B7FFF","#C7C7C7FF")���Ȼ�
pal=c("#6BAED6FF","#6BAED6FF","#6BAED6FF","#C7C7C7FF","#6BAED6FF","#6BAED6FF","#6BAED6FF","#6BAED6FF","#6BAED6FF","#6BAED6FF","#6BAED6FF","#F05C3BFF")
#library(RColorBrewer)
#colCount=ncol(rt)
#getPalette=colorRampPalette(brewer.pal(12,"Set3"),alpha=0.5)

#����Ȧͼ
pdf(file="CHECKPOINTcircos.pdf", width=7, height=7)
par(mar=c(2,2,2,4))
circos.par(gap.degree=c(3,rep(2, nrow(cor1)-1)), start.degree = 180)
chordDiagram(cor1, grid.col=pal, col=col1, transparency = 0.5, symmetric = T)
#chordDiagram(cor1, grid.col=rainbow(ncol(rt),alpha = 0.5), col=col1, transparency = 0.5, symmetric = T)
par(xpd=T)
colorlegend(col, vertical = T,labels=c(1,0,-1),xlim=c(1.1,1.3),ylim=c(-0.4,0.4)) #����ͼ��
dev.off()
circos.clear()



################################## ��������Ȧͼ#########################################
#��ȡ��������ļ�
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt, check.names=F)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#ɾ��������Ʒ
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]
#��ȡĿ����������
gene=CF
sameGene=intersect(as.vector(gene), rownames(data))
geneExp=data[sameGene,]
exp=t(geneExp)
row.names(exp)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\",  row.names(exp))
exp=avereps(exp)
data=exp
data=log2(data+1)
#��ȡ�����ļ�
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#���ݺϲ�
sameSample=intersect(row.names(data), row.names(risk))
data=cbind(risk[sameSample,"riskScore",drop=F], data[sameSample,,drop=F])

#�ϲ�����ͻ������ļ�
#TMB=read.table(tmbFile, header=T, sep="\t", check.names=F, row.names=1)
#sameSample=intersect(row.names(TMB), row.names(data))
#rt=cbind(TMB[sameSample,,drop=F], data[sameSample,,drop=F])
rt=data
#�����������Ծ���
cor1=cor(rt)
res1=cor.mtest(rt, conf.level=0.95)
res1[["p"]]
#����ͼ����ɫ
col = c(rgb(1,0.6,0.45,seq(1,0,length=32)),rgb(0.4,0.6,0.95,seq(0,1,length=32)))
cor1[cor1==1]=0
c1 = ifelse(c(cor1)>=0,rgb(1,0.6,0.45,abs(cor1)),rgb(0.4,0.6,0.95,abs(cor1)))
col1 = matrix(c1,nc=ncol(rt))

#c("#6BAED6FF","#F39B7FFF","#C7C7C7FF")���Ȼ�
pal=c("#6BAED6FF","#6BAED6FF","#6BAED6FF","#6BAED6FF","#6BAED6FF","#F05C3BFF")
#library(RColorBrewer)
#colCount=ncol(rt)
#getPalette=colorRampPalette(brewer.pal(12,"Set3"),alpha=0.5)

#����Ȧͼ
pdf(file="CFcircos.pdf", width=7, height=7)
par(mar=c(2,2,2,4))
circos.par(gap.degree=c(3,rep(2, nrow(cor1)-1)), start.degree = 180)
chordDiagram(cor1, grid.col=pal, col=col1, transparency = 0.5, symmetric = T)
#chordDiagram(cor1, grid.col=rainbow(ncol(rt),alpha = 0.5), col=col1, transparency = 0.5, symmetric = T)
par(xpd=T)
colorlegend(col, vertical = T,labels=c(1,0,-1),xlim=c(1.1,1.3),ylim=c(-0.4,0.4)) #����ͼ��
dev.off()
circos.clear()