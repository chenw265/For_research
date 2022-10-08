
#install.packages("survival")
#install.packages("survminer")


#引用包
library(survival)
library(dplyr)
library(survminer)
library(gridExtra)
tmbFile="TMB.txt"              #肿瘤突变负荷文件
riskFile="risk.TCGA.txt"       #风险文件
setwd("E:\\My_projects\\LUSC\\14TMB")       #修改工作目录

#读取输入文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)      #读取风险文件
tmb=read.table(tmbFile, header=T, sep="\t", check.names=F, row.names=1)        #读取TMB数据文件

#合并数据
sameSample=intersect(row.names(tmb), row.names(risk))
tmb=tmb[sameSample,,drop=F]
risk=risk[sameSample,,drop=F]
data=cbind(risk, tmb)

#获取最优cutoff
res.cut=surv_cutpoint(data, time = "futime", event = "fustat", variables =c("TMB"))
cutoff=as.numeric(res.cut$cutpoint[1])
tmbType=ifelse(data[,"TMB"]<=cutoff, "L-TMB", "H-TMB")
scoreType=ifelse(data$Risk=="low", "low risk", "high risk")
mergeType=paste0(tmbType, "+", scoreType)

#定义生存分析函数
bioSurvival=function(surData=null, outFile=null){
	diff=survdiff(Surv(futime, fustat) ~ group, data=surData)
	length=length(levels(factor(surData[,"group"])))
	pValue=1-pchisq(diff$chisq, df=length-1)
	if(pValue<0.001){
		pValue="p<0.001"
	}else{
		pValue=paste0("p=",sprintf("%.03f",pValue))
	}
	fit <- survfit(Surv(futime, fustat) ~ group, data = surData)
	#print(surv_median(fit))
	
	#绘制生存曲线
	bioCol=c("#00468BFF", "#ED0000FF","#42B540FF","#0099B4FF","#925E9FFF","#DF8F44FF")
	bioCol=bioCol[1:length]
	surPlot=ggsurvplot(fit, 
			           data=surData,
			           conf.int=F,
			           #pval=pValue,
			           pval.size=6,
			           legend.title="",
			           legend.labs=levels(factor(surData[,"group"])),
			           font.legend=10,
			           legend = c(0.8, 0.8),
			           xlab="Time(years)",
			           break.time.by = 1,
			           palette = bioCol,
			           surv.median.line = "hv",
			           risk.table=T,
			           cumevents=F,
			           risk.table.height=.25)
	#输出图形
	pdf(file=outFile, onefile = FALSE, width=5.5, height=4.8)
	print(surPlot)
	dev.off()
}

#绘制肿瘤突变负荷的生存曲线
data$group=tmbType
bioSurvival(surData=data, outFile="TMB.survival.pdf")

#绘制肿瘤突变负荷联合病人风险的生存曲线
library(dplyr)
library(gridExtra)
canc.plots<-list()
data$group=mergeType
write.csv(data, "data.csv")
rt=data
survos<-Surv(rt$futime,rt$fustat=='1')
os<-survfit(survos~rt$group)
oslogrank<-survdiff(survos~rt$group)
bioCol=c("#00468BFF", "#ED0000FF","#42B540FF","#0099B4FF","#925E9FFF","#DF8F44FF")
canc.plots[[1]]<-ggsurvplot(
  os, 
  data=data,
  conf.int=F,
  #pval=pValue,
  pval.size=6,
  legend.title="",
  legend.labs=levels(factor(data[,"group"])),
  font.legend=10,
  legend = c(0.8, 0.8),
  xlab="Time(years)",
  break.time.by = 1,
  palette = bioCol,
  surv.median.line = "hv",
  risk.table=T,
  cumevents=F,
  risk.table.height=.25
)

OS.compared <- pairwise_survdiff(Surv(futime, fustat) ~ group,
                                 p.adjust.method = "none",
                                 data = rt)
OS.pvalue<-round(as.data.frame(OS.compared$p.value),4)
OS.pvalue<-as.data.frame(t(OS.pvalue))
#colnames(OS.pvalue)<-c("C2","C3")
#rownames(OS.pvalue)<-c("C1","C2")
OS.pvalue$" "<-row.names(OS.pvalue)
OS.pvalue1<-cbind(OS.pvalue$" ",OS.pvalue[1:3])
OS.pvalue1 <- OS.pvalue1 %>% 
  rename(
    '  ' = 'OS.pvalue$" "'
  )

ggplot.pureplots<-lapply(canc.plots, function(x){
  pureplot<-x$plot+
    theme(plot.title = element_text(margin=margin(0,0,0,0)))+
    theme(plot.margin = unit(c(0,2,0,0),"mm"),
          plot.title.position = "plot",
          plot.title = element_text(hjust=0.05)
    )+
    scale_x_continuous(expand = c(0.08,0),limits = c(-30,125),n.breaks = 10,  
                       breaks = seq(from = 0, to = 10, by = 1), 
                       labels = seq(from = 0, to = 10, by = 1))+
    scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.2),
                       #labels = seq(from = 0.00, to = 1.00, by = 0.20)
                       labels = c("0.00","0.20","0.40","0.60","0.80","1.00"))
  return(pureplot)
})

a<-tableGrob(OS.pvalue1,
             rows= NULL,
             theme = ttheme_default(core = list(fg_params=list(fontface=rep("bold"),cex = 0.6),
                                                bg_params = list(fill= "white")),
                                    colhead = list(fg_params=list(fontface=rep("bold"),cex = 0.6),
                                                   bg_params = list(fill= "white")),
                                    rowhead = list(fg_params=list(fontface=rep("bold")))))
a$heights <- unit(rep(0.38, 4), "cm")
a$widths <- unit(c(2,2,2,2), "cm")
ggplot.pureplots[[1]]<-ggplot.pureplots[[1]]+annotation_custom(a, 
                                                               xmin=0.8, xmax=4, ymin=-0.01, ymax=0.1)
pdf(file="TMB-risk.survival.pdf", onefile = FALSE, width=6.5, height=4.8)
ggplot.pureplots[[1]]
dev.off()



