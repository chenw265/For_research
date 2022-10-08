
#���ð�
library(limma)
library(ggplot2)
library(dplyr)
corFilter=0.4                  #���ϵ����������
fdrFilter=0.05                #FDR���˱�׼
gene1File="gene1.txt"      #�����б��ļ�1
gene2File="gene2.txt"                #�����б��ļ�2
expFile="symbol.txt"      #Ԥ����ػ����б��ļ�
setwd("E:\\My_projects\\LUSC\\24modelGenesANDpyroptosisGenes")     #���ù���Ŀ¼

#��ȡ�����ļ����������ݽ��д���
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
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
data=t(data)
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
data=avereps(data)
data=log2(data+1)
data=t(data)

#��ȡϸ����������ı�����
gene1=read.table(gene1File, header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene1[,1]), rownames(data))
geneExp1=data[sameGene,]

gene2=read.table(gene2File, header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene2[,1]), rownames(data))
geneExp2=data[sameGene,]

#����Լ���
outTab=data.frame()
for(i in row.names(geneExp2)){
	if(sd(geneExp2[i,])>0.01){
		for(j in row.names(geneExp1)){
			if(i==j){next}
			x=as.numeric(geneExp2[i,])
			y=as.numeric(geneExp1[j,])
			corT=cor.test(x, y)
			cor=corT$estimate
			pvalue=corT$p.value
			text=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
			outTab=rbind(outTab,cbind(gene1s=i, gene2s=j, cor, text, pvalue))
		}
	}
}

#�����������ͼ
outTab$cor=as.numeric(outTab$cor)
#outTab$cor=round(outTab$cor,2)
pdf(file="cor.pdf", width=5, height=15)
ggplot(outTab, aes(gene1s, gene2s)) + 
  geom_tile(aes(fill = cor), colour = "black", size = NA)+ 
  scale_fill_gradient2(low = "#466983FF", high = "#C75127FF") + 
  geom_text(aes(label=text),col ="black",size = 3) +
  #theme_minimal() +    #ȥ������
  theme_classic(base_size = 12)+
  theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(),axis.ticks.y=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 10, face = "bold"),   #x������
        axis.text.y = element_text(size = 10, face = "bold")) +       #y������
  labs(fill =paste0("***  p<0.001","\n", "**  p<0.01","\n", " *  p<0.05","\n", "\n","Correlation")) +   #����ͼ��
  #labs(fill =paste0("***  p<0.001","\n", "**  p<0.01","\n", " *  p<0.05","\n", "\n","Correlation")) +   #����ͼ��
  scale_x_discrete(position = "bottom")      #X��������ʾλ��
dev.off()
