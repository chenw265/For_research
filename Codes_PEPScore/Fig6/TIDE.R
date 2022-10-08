
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggpubr")
#install.packages("Rmisc")

#���ð�
library(limma)
library(ggpubr)
library(ggplot2)
library(tidyr)
library(Rmisc)
library(reshape2)
library(survival)
library(survminer)
library(timeROC)
library(ggsci)
library(scales)
expFile="symbol.txt"
tideFile="TIDE.txt"          #TIDE�ļ�
riskFile="risk.TCGA.txt"     #�����ļ�
tisFile="TISgene.txt"            #TIS�����ļ�
setwd("E:\\My_projects\\LUSC\\13TIDE")     #���ù���Ŀ¼

#��ȡ�����ļ�
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

#��ȡTIS���������
geneRT=read.table(tisFile, header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(geneRT[,1]), row.names(data))
geneExp=data[sameGene,]
logData=log2(geneExp+1)
TISscore=colMeans(logData)

#��ȡTIDE����ļ�
TIDEscore=read.table(tideFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(TIDEscore), names(TISscore))
score=cbind(TIDEscore[sameSample,,drop=F], TIS=TISscore[sameSample])
score=score[,c("TIDE","TIS")]

#ȥ��������Ʒ
group=sapply(strsplit(row.names(score),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
score=score[group==0,]
score=avereps(score)
row.names(score)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(score))
score=avereps(score)

#��ȡ���������ļ�
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(score), row.names(risk))
score=score[sameSample, , drop=F]
risk=risk[sameSample, c("futime","fustat","riskScore"), drop=F]
rt=cbind(risk, score)

######����ģ�ͱȽϵ�ROC����######
predictTime=3     #����Ԥ������
aucText=c()
#bioCol=rainbow(ncol(rt)-2, s=0.9, v=0.9)
#bioCol=pal_lancet()(ncol(rt-2))
bioCol=c("#ED0000FF","#42B540FF","#00468BFF")

pdf(file="multiROC.pdf", width=5.5, height=5.5)
#���Ʒ��յ÷ֵ�ROC����
i=3
ROC_rt=timeROC(T=risk$futime,
               delta=risk$fustat,
               marker=risk$riskScore, cause=1,
               weighting='aalen',
               times=c(predictTime),ROC=TRUE)
plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=3)
aucText=c(paste0("Risk", ", AUC=", sprintf("%.3f",ROC_rt$AUC[2])))
abline(0,1)
#��TIS��ֺ�TIDE��ֽ���ROC���ߵĻ���
for(i in 4:ncol(rt)){
  ROC_rt=timeROC(T=rt$futime,
                 delta=rt$fustat,
                 marker=rt[,i], cause=1,
                 weighting='aalen',
                 times=c(predictTime),ROC=TRUE)
  plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=3, add=TRUE)
  aucText=c(aucText, paste0(colnames(rt)[i],", AUC=",sprintf("%.3f",ROC_rt$AUC[2])))
}
#����ͼ�����õ�ROC�����µ����
legend("bottomright", aucText,lwd=2,bty="n",col=bioCol[1:(ncol(rt)-1)])
dev.off()

#################�ְ�С����ͼ############################
	
#�ϲ���������
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(score), row.names(risk))
score=score[sameSample, , drop=F]
risk=risk[sameSample, "Risk", drop=F]
data=cbind(score, risk)
	
#���ñȽ���
data$Risk=ifelse(data$Risk=="high", "High-risk", "Low-risk")
group=levels(factor(data$Risk))
data$Risk=factor(data$Risk, levels=c("Low-risk", "High-risk"))
group=levels(factor(data$Risk))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#��tide��ֽ���ѭ��,�ֱ����С����ͼ
for(i in colnames(data)[1:(ncol(data)-1)]){
	gg1=ggviolin(data, x="Risk", y=i, fill = "Risk", 
	         xlab="", ylab=i,
	         palette=c("#00468BFF", "#ED0000FF"),
	         legend.title="Risk",
	         add = "boxplot")+ 
	         stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
	
	pdf(file=paste0("violin.", i, ".pdf"), width=6, height=5)
	print(gg1)
	dev.off()
}

#�Ը���ֽ��л��Ʒְ�С����ͼ
data1=gather(data,"category","Score",-Risk)
Data_summary <- summarySE(data1, measurevar="Score", groupvars=c("Risk","category"))

#�ְ�С����ͼ��ͼ��������
#Source:https://gist.github.com/Karel-Kroeze/746685f5613e01ba820a31e57f87ec87
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, draw_group = function(self, data, ..., draw_quantiles = NULL){
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1,'group']
  newdata <- plyr::arrange(transform(data, x = if(grp%%2==1) xminv else xmaxv), if(grp%%2==1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 
                                              1))
    quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
  }
})
geom_split_violin <- function (mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position, show.legend = show.legend, inherit.aes = inherit.aes, params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

#mytheme
if(T){
  mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                   axis.title = element_text(size = 12,color ="black"), 
                   axis.text = element_text(size= 12,color = "black"),
                   panel.grid.minor.y = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   axis.text.x = element_text(angle = 45, hjust = 1 ),
                   panel.grid=element_blank(),
                   legend.position = "top",
                   legend.text = element_text(size= 12),
                   legend.title= element_text(size= 12)
  ) 
}


gene_split_violin <- ggplot(data1,aes(x=category ,y=Score,fill=Risk))+
  geom_split_violin(trim= F,color="white",scale = "area") + #���Ʒְ��С����ͼ
  geom_point(data = Data_summary,aes(x=category ,y=Score),pch=19,
             position=position_dodge(0.5),size= 1)+ #���ƾ�ֵΪ��ͼ
  geom_errorbar(data = Data_summary,aes(ymin = Score-ci, ymax= Score+ci), 
                width= 0.05, 
                position= position_dodge(0.5), 
                color="black",
                alpha = 0.8,
                size= 0.5) +
  scale_fill_manual(values = c("#56B4E9", "#E69F00"))+ 
  labs(y=("Score"),x=NULL,title = "") + 
  theme_bw()+ mytheme +
  stat_compare_means(aes(group = Risk),
                     label = "p.signif",
                     method = "anova",
                     label.y = max(data1$Score),
                     hide.ns = T)
gene_split_violin;ggsave(gene_split_violin,
                         filename = "./gene_split_violin.pdf",
                         height = 10,width = 16,units = "cm")