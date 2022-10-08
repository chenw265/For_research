#install.packages("survival")
#install.packages("survminer")


#���ð�
library(survival)
library(survminer)
ClusterFile="cluster.txt"     #���ͽ���ļ�
cliFile="time.txt"            #���������ļ�
setwd("E:\\My_projects\\LUSC\\03ClusterSurv")      #���ù���Ŀ¼

#��ȡ�����ļ�
Cluster=read.table(ClusterFile, header=T, sep="\t", check.names=F, row.names=1)
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
colnames(cli)=c("futime", "fustat")
cli$futime=cli$futime/365

#���ݺϲ�
sameSample=intersect(row.names(Cluster), row.names(cli))
rt=cbind(cli[sameSample,,drop=F], Cluster[sameSample,,drop=F])

#����������
length=length(levels(factor(rt$Cluster)))
diff=survdiff(Surv(futime, fustat) ~ Cluster, data = rt)
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
	pValue="p<0.001"
}else{
	pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ Cluster, data = rt)
#print(surv_median(fit))

#������������
#bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=c("#0099B4FF", "#ED0000FF")
bioCol=bioCol[1:length]
surPlot=ggsurvplot(fit, 
		           data=rt,
		           conf.int=T,
		           pval=pValue,
		           pval.size=6,
		           legend.title="Cluster",
		           legend.labs=levels(factor(rt[,"Cluster"])),
		           legend = c(0.8, 0.8),
		           font.legend=10,
		           xlab="Time(years)",
		           break.time.by = 1,
		           palette = bioCol,
		           surv.median.line = "hv",
		           risk.table=T,
		           cumevents=F,
		           risk.table.height=.25)
pdf(file="survival.pdf",onefile = FALSE,width=7,height=5.5)
print(surPlot)
dev.off()
