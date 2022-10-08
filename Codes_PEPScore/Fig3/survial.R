
#���ð�
library(survival)
library(survminer)
setwd("E:\\My_projects\\LUSC\\08model_2\\01_05")       #���ù���Ŀ¼

#�����������ߺ���
bioSurvival=function(inputFile=null,outFile=null){
	#��ȡ�����ļ�
	rt=read.table(inputFile, header=T, sep="\t", check.names=F)
	#�Ƚϸߵͷ�����������죬�õ�������pֵ
	diff=survdiff(Surv(futime, fustat) ~ Risk,data = rt)
	pValue=1-pchisq(diff$chisq,df=1)
	if(pValue<0.001){
		pValue="p<0.001"
	}else{
		pValue=paste0("p=",sprintf("%.03f",pValue))
	}
	fit <- survfit(Surv(futime, fustat) ~ Risk, data = rt)
		
	#������������
	surPlot=ggsurvplot(fit, 
		           data=rt,
		           conf.int=F,
		           pval=pValue,
		           pval.size=6,
		           legend.title="Risk",
		           legend.labs=c("High risk", "Low risk"),
		           xlab="Time(years)",
		           break.time.by = 1,
		           palette=c("#ED0000FF","#00468BFF"),
		           surv.median.line = "hv",
		           risk.table=TRUE,
		           risk.table.title="",
		           risk.table.height=.25)
	pdf(file=outFile,onefile = FALSE,width = 6.5,height =5.5)
	print(surPlot)
	dev.off()
}

#���ú���,������������
bioSurvival(inputFile="risk.TCGA.txt", outFile="survival.TCGA.pdf")
bioSurvival(inputFile="risk.GEO5.txt", outFile="survival.GEO5.pdf")
bioSurvival(inputFile="risk.GEO1.txt", outFile="survival.GEO1.pdf")




bioSurvival(inputFile="risk.GEO3.txt", outFile="survival.GEO3.pdf")
bioSurvival(inputFile="risk.GEO2.txt", outFile="survival.GEO2.pdf")
bioSurvival(inputFile="risk.GEO4.txt", outFile="survival.GEO4.pdf")
