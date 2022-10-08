#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")


#���ð�
library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

expFile="symbol.txt"          #���������ļ�
riskFile="risk.TCGA.txt"      #�����ļ�
gmtFile="h.all.v7.5.1.symbols.gmt"     #�����ļ�
setwd("E:\\My_projects\\LUSC\\10GSEA")      #���ù���Ŀ¼

#��ȡ�ļ�,���������ļ���������
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.5,]

#ɾ��������Ʒ
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]
data=t(data)
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
data=t(avereps(data))

#��ȡ�����ļ�
Risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
data=data[,row.names(Risk)]

#�ߵͷ��ձȽϣ��õ�logFC
dataL=data[,row.names(Risk[Risk[,"Risk"]=="low",])]
dataH=data[,row.names(Risk[Risk[,"Risk"]=="high",])]
meanL=rowMeans(dataL)
meanH=rowMeans(dataH)
meanL[meanL<0.00001]=0.00001
meanH[meanH<0.00001]=0.00001
logFC=log2(meanH)-log2(meanL)
logFC=sort(logFC,decreasing=T)
genes=names(logFC)

#��������ļ�
gmt=read.gmt(gmtFile)

#��������
kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$pvalue<0.05,]
write.table(kkTab,file="GSEA.result.txt",sep="\t",quote=F,row.names = F)
	
#����߷��ո�����ͼ��
termNum=5     #չʾǰ5��ͨ·
kkUp=kkTab[kkTab$NES>0,]
if(nrow(kkUp)>=termNum){
	showTerm=row.names(kkUp)[1:termNum]
	gseaplot=gseaplot2(kk, showTerm, base_size=8, title="Enriched in high risk group")
	pdf(file="GSEA.highRisk.pdf", width=7.5, height=5.5)
	print(gseaplot)
	dev.off()
}

#����ͷ��ո�����ͼ��
termNum=5     #չʾǰ5��ͨ·
kkDown=kkTab[kkTab$NES<0,]
if(nrow(kkDown)>=termNum){
	showTerm=row.names(kkDown)[1:termNum]
	gseaplot=gseaplot2(kk, showTerm, base_size=8, title="Enriched in low risk group")
	pdf(file="GSEA.lowRisk.pdf", width=7.5, height=5.5)
	print(gseaplot)
	dev.off()
}
