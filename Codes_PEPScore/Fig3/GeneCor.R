
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("corrplot")


#���ð�
library(limma)
library(corrplot)

expFile="tcga.pyroptosisExp.txt"             #���������ļ�
geneFile="gene.txt"              #�����б��ļ�
riskFile="risk.TCGA.txt"      #�����ļ�
surFile="tcgatime.txt"
setwd("E:\\My_projects\\LUSC\\20RiskPyroptosis")     #���ù���Ŀ¼

#��ȡ��������ļ�,�������ݽ��д���
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

#ɾ��������Ʒ
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]

#��ȡĿ����������
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),]
data=t(data)
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
data=avereps(data)
data=log2(data+1)

#��ȡ�����ļ�
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#���ݺϲ�
sameSample=intersect(row.names(data), row.names(risk))
data=cbind(risk[sameSample,"Risk",drop=F], data[sameSample,,drop=F])

#��ȡ�����ļ�
#sur=read.table(surFile, header = T, sep="\t", check.names = F, row.names = 1)
#���ݺϲ�
#sameSample=intersect(row.names(data), row.names(sur))
#data=cbind(sur[sameSample, "futime", drop=F], data[sameSample,,drop=F])

#�ֱ��ȡ�ߵͷ��յĻ�������ļ�
dataH=data[which(data$Risk=="high"),2:ncol(data)]
dataL=data[which(data$Risk=="low"),2:ncol(data)]

#����Ծ���
Mh=cor(dataH)
resH=cor.mtest(dataH, conf.level=0.95)
Mh[which(Mh==1)]=0

Ml=cor(dataL)
resL=cor.mtest(dataL, conf.level=0.95)

#���������ͼ��
pdf(file="HighRiskgeneCor.pdf", width=7, height=7)
corrplot(Mh,
         order="original",
         method = "color",
         type = "upper",
         tl.cex=0.6, pch=T,
         p.mat = resH$p,
         insig = "label_sig",
         pch.cex = 1.6,
         sig.level=0.05,
         number.cex = 1,
         col=colorRampPalette(c("#0099B4FF", "white", "#ED0000FF"))(50),
         tl.col="black",
         add = F)
dev.off()

pdf(file="LowRiskgeneCor.pdf", width=7, height=7)
corrplot(Ml,
         order="original",
         method = "color",
         type = "lower",
         tl.cex=0.6, pch=T,
         p.mat = resL$p,
         insig = "label_sig",
         pch.cex = 1.6,
         sig.level=0.05,
         number.cex = 1,
         col=colorRampPalette(c("#0099B4FF", "white", "#ED0000FF"))(30),
         tl.col="black",
         add = F)
dev.off()



library(igraph)
library(reshape2)
cutoff=0.2                      #�������ֵ

out="Hcorrelation.pdf"           #�������ļ�
#��������Ծ����һ��
mydata = Mh
upper = upper.tri(mydata)
mydata[upper] = NA

#������Ծ���ת��Ϊ���ݿ�
df = data.frame(gene=rownames(mydata),mydata)
dfmeltdata = melt(df,id="gene")
dfmeltdata = dfmeltdata[!is.na(dfmeltdata$value),]
dfmeltdata = dfmeltdata[dfmeltdata$gene!=dfmeltdata$variable,]
dfmeltdata = dfmeltdata[abs(dfmeltdata$value)>cutoff,]

#��������ͼ�Ľڵ�ͱ�
corweight = dfmeltdata$value
weight = corweight+abs(min(corweight))+5
d = data.frame(p1=dfmeltdata$gene,p2=dfmeltdata$variable,weight=dfmeltdata$value)
g = graph.data.frame(dfmeltdata,directed = FALSE)

#������ɫ���ڵ��С�������С
E(g)$weight = weight
E(g)$color = ifelse(corweight>0,rgb(254/255,67/255,101/255,abs(corweight)),rgb(0/255,0/255,255/255,abs(corweight)))
V(g)$size = 8
V(g)$shape = "circle"
V(g)$lable.cex = 1.2
V(g)$color = "white"

#���ӻ�
pdf(out, width=7, height=6)
layout(matrix(c(1,1,1,0,2,0),byrow=T,nc=3),height=c(6,1),width=c(3,4,3))
par(mar=c(1.5,2,2,2))
vertex.frame.color = NA
plot(g,layout=layout_with_kk,
     vertex.label.cex=V(g)$lable.cex,
     edge.width = E(g)$weight,edge.arrow.size=0,
     vertex.label.color="black",
     vertex.frame.color=vertex.frame.color,
     edge.color=E(g)$color,
     vertex.label.cex=V(g)$lable.cex,
     vertex.label.font=2,
     vertex.size=V(g)$size,
     edge.curved=0.4)

#����ͼ��
color_legend = c(rgb(254/255,67/255,101/255,seq(1,0,by=-0.01)),rgb(0/255,0/255,255/255,seq(0,1,by=0.01)))
par(mar=c(2,2,1,2),xpd = T,cex.axis=1.6,las=1)
barplot(rep(1,length(color_legend)),border = NA, space = 0,ylab="",xlab="",xlim=c(1,length(color_legend)),horiz=FALSE,
        axes = F, col=color_legend,main="")
axis(3,at=seq(1,length(color_legend),length=5),c(1,0.5,0,-0.5,-1),tick=FALSE)
dev.off()




out="Lcorrelation.pdf"           #�������ļ�
#��������Ծ����һ��
mydata = Ml
upper = upper.tri(mydata)
mydata[upper] = NA

#������Ծ���ת��Ϊ���ݿ�
df = data.frame(gene=rownames(mydata),mydata)
dfmeltdata = melt(df,id="gene")
dfmeltdata = dfmeltdata[!is.na(dfmeltdata$value),]
dfmeltdata = dfmeltdata[dfmeltdata$gene!=dfmeltdata$variable,]
dfmeltdata = dfmeltdata[abs(dfmeltdata$value)>cutoff,]

#��������ͼ�Ľڵ�ͱ�
corweight = dfmeltdata$value
weight = corweight+abs(min(corweight))+5
d = data.frame(p1=dfmeltdata$gene,p2=dfmeltdata$variable,weight=dfmeltdata$value)
g = graph.data.frame(dfmeltdata,directed = FALSE)

#������ɫ���ڵ��С�������С
E(g)$weight = weight
E(g)$color = ifelse(corweight>0,rgb(254/255,67/255,101/255,abs(corweight)),rgb(0/255,0/255,255/255,abs(corweight)))
V(g)$size = 8
V(g)$shape = "circle"
V(g)$lable.cex = 1.2
V(g)$color = "white"

#���ӻ�
pdf(out, width=7, height=6)
layout(matrix(c(1,1,1,0,2,0),byrow=T,nc=3),height=c(6,1),width=c(3,4,3))
par(mar=c(1.5,2,2,2))
vertex.frame.color = NA
plot(g,layout=layout_with_kk,
     vertex.label.cex=V(g)$lable.cex,
     edge.width = E(g)$weight,edge.arrow.size=0,
     vertex.label.color="black",
     vertex.frame.color=vertex.frame.color,
     edge.color=E(g)$color,
     vertex.label.cex=V(g)$lable.cex,
     vertex.label.font=2,
     vertex.size=V(g)$size,
     edge.curved=0.4)

#����ͼ��
color_legend = c(rgb(254/255,67/255,101/255,seq(1,0,by=-0.01)),rgb(0/255,0/255,255/255,seq(0,1,by=0.01)))
par(mar=c(2,2,1,2),xpd = T,cex.axis=1.6,las=1)
barplot(rep(1,length(color_legend)),border = NA, space = 0,ylab="",xlab="",xlim=c(1,length(color_legend)),horiz=FALSE,
        axes = F, col=color_legend,main="")
axis(3,at=seq(1,length(color_legend),length=5),c(1,0.5,0,-0.5,-1),tick=FALSE)
dev.off()