#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggplot2")
#install.packages("ggpubr")
#install.packages("ggExtra")


#引用包
library(limma)
library(ggplot2)
library(ggpubr)
library(ggExtra)
setwd("E:\\My_projects\\LUSC\\14TMB")      #设置工作目录

#读取肿瘤突变负荷文件
tmb=read.table("TMB.txt", header=T, sep="\t", check.names=F, row.names=1)
	
#读取风险数据文件
risk=read.table("risk.TCGA.txt", header=T, sep="\t", check.names=F, row.names=1)
	
#合并数据
sameSample=intersect(row.names(tmb), row.names(risk))
tmb=tmb[sameSample,,drop=F]
risk=risk[sameSample,,drop=F]
data=cbind(tmb, risk)
data$TMB[data$TMB>quantile(data$TMB,0.99)]=quantile(data$TMB,0.99)
	
#设置比较组
data$Risk=ifelse(data$Risk=="high", "High-risk", "Low-risk")
group=levels(factor(data$Risk))
data$Risk=factor(data$Risk, levels=c("Low-risk", "High-risk"))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
	
#绘制boxplot
boxplot=ggboxplot(data, x="Risk", y="TMB", color="Risk",
			      xlab="",
			      ylab="Tumor tmbation burden",
			      legend.title="",
			      palette = c("#00468BFF", "#ED0000FF"),
			      size = 1,
			      notch = T,
			      add = "jitter")+ 
	stat_compare_means(comparisons = my_comparisons)
	
#输出图片
pdf(file="riskTMB.pdf", width=3.5, height=4)
print(boxplot)
dev.off()

#相关性分析
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
#输出相关性图形
pdf(file="cor.pdf", width=5.2, height=5)
print(p2)
dev.off()
