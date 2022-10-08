#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggplot2")
#install.packages("ggpubr")
#install.packages("ggExtra")


#���ð�
library(limma)
library(ggplot2)
library(ggpubr)
library(ggExtra)
setwd("E:\\My_projects\\LUSC\\14TMB")      #���ù���Ŀ¼

#��ȡ����ͻ�为���ļ�
tmb=read.table("TMB.txt", header=T, sep="\t", check.names=F, row.names=1)
	
#��ȡ���������ļ�
risk=read.table("risk.TCGA.txt", header=T, sep="\t", check.names=F, row.names=1)
	
#�ϲ�����
sameSample=intersect(row.names(tmb), row.names(risk))
tmb=tmb[sameSample,,drop=F]
risk=risk[sameSample,,drop=F]
data=cbind(tmb, risk)
data$TMB[data$TMB>quantile(data$TMB,0.99)]=quantile(data$TMB,0.99)
	
#���ñȽ���
data$Risk=ifelse(data$Risk=="high", "High-risk", "Low-risk")
group=levels(factor(data$Risk))
data$Risk=factor(data$Risk, levels=c("Low-risk", "High-risk"))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
	
#����boxplot
boxplot=ggboxplot(data, x="Risk", y="TMB", color="Risk",
			      xlab="",
			      ylab="Tumor tmbation burden",
			      legend.title="",
			      palette = c("#00468BFF", "#ED0000FF"),
			      size = 1,
			      notch = T,
			      add = "jitter")+ 
	stat_compare_means(comparisons = my_comparisons)
	
#���ͼƬ
pdf(file="riskTMB.pdf", width=3.5, height=4)
print(boxplot)
dev.off()

#����Է���
xlab="riskScore"
ylab="TMB"
x=as.numeric(data[,xlab])
y=as.numeric(data[,ylab])
df1=as.data.frame(cbind(x,y))
p1=ggplot(df1, aes(x, y)) + 
		xlab("Risk score") + ylab("Tumor tmbation burden")+
		geom_point() + geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
		stat_cor(method = 'spearman', aes(x =x, y =y))
p2=ggMarginal(p1, type="density", xparams=list(fill = "orange"), yparams=list(fill = "blue"))
#��������ͼ��
pdf(file="cor.pdf", width=5.2, height=5)
print(p2)
dev.off()