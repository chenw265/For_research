#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("reshape2")
#install.packages("ggpubr")


#引用包
library(limma)
library(reshape2)
library(ggpubr)

riskFile="risk.TCGA.txt"            #风险文件
immFile="CIBERSORT-Results.txt"     #免疫细胞浸润结果文件
pFilter=0.05                        #免疫细胞浸润结果的过滤条件
setwd("E:\\My_projects\\LUSC\\12ImmuneRelated\\02ImmuneCor")      #设置工作目录

#读取免疫细胞结果文件，并对数据进行整理
immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)
immune=immune[immune[,"P-value"]<pFilter,]
data=as.matrix(immune[,1:(ncol(immune)-3)])

#删除正常样品
group=sapply(strsplit(row.names(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[group==0,]
row.names(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(data))
data=avereps(data)

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(data), row.names(risk))
rt=cbind(data[sameSample,,drop=F], risk[sameSample,"Risk",drop=F])
rt=rt[order(rt$Risk, decreasing=T),]
conNum=nrow(rt[rt$Risk=="low",])
treatNum=nrow(rt[rt$Risk=="high",])

##########绘制柱状图##########
data=t(rt[,-ncol(rt)])
pdf("barplot.pdf", height=10, width=18)
col=rainbow(nrow(data), s=0.5, v=0.5)
par(las=1,mar=c(8,5,4,20),mgp=c(3,0.1,0),cex.axis=1.5)
a1=barplot(data,col=col,yaxt="n",ylab="Relative Percent",xaxt="n",cex.lab=1.8)
a2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
par(srt=0,xpd=T)
rect(xleft = a1[1], ybottom = -0.01, xright = a1[conNum], ytop = -0.06,col="#00468BFF")
text(a1[conNum]/2,-0.035,"Low risk",cex=2)
rect(xleft = a1[conNum], ybottom = -0.01, xright =a1[length(a1)] , ytop = -0.06,col="#ED0000FF")
text((a1[length(a1)]+a1[conNum])/2,-0.035,"High risk",cex=2)
ytick2 = cumsum(data[,ncol(data)])
ytick1 = c(0,ytick2[-length(ytick2)])
legend(par('usr')[2]*0.96,par('usr')[4],legend=rownames(data),col=col,pch=15,bty="n",cex=1.65)
dev.off()


##################绘制箱线图##################
#把数据转换成ggplot2输入文件
data=rt
data=melt(data, id.vars=c("Risk"))
colnames(data)=c("Risk", "Immune", "Expression")
#绘制箱线图
group=levels(factor(data$Risk))
data$Risk=factor(data$Risk, levels=c("low","high"))
bioCol=c("#00468BFF", "#ED0000FF")
bioCol=bioCol[1:length(group)]
boxplot=ggboxplot(data, x="Immune", y="Expression", fill="Risk",
				  xlab="",
				  ylab="Fraction",
				  legend.title="Risk",
				  width=0.8,
				  palette=bioCol)+
				  rotate_x_text(50)+
	stat_compare_means(aes(group=Risk),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "")), label="p.signif")
#输出图片
pdf(file="immune.diff.pdf", width=8, height=6)
print(boxplot)
dev.off()
