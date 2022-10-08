#引用包
library(limma)
library(survival)
library(survminer)
library(timeROC)
library(ggpubr)


expFile="symbol.txt"            #表达数据文件
riskFile="risk.TCGA.txt"     #风险文件
geneFiles="Li signature.txt"     #模型基因文件
setwd("E:\\My_projects\\LUSC\\22modelCompare")               #设置工作目录

#读取表达文件，并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#删掉正常样品
group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2", "1", group)
data=data[,group==0]

#读取风险文件
riskRT=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
riskRT=riskRT[,c("futime","fustat","riskScore")]
colnames(riskRT)=c("futime","fustat","ALL signature")


model=read.table(geneFiles, header = T, sep = "\t", check.names = F)
modelgene=model$gene
samegene=intersect(modelgene,rownames(data))
data1=data[samegene,]
colnames(data1)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", colnames(data1))
data1=t(data1)
data1=avereps(data1)
data1=as.data.frame(data1)

#合并生存数据文件
cli=riskRT[,c("futime", "fustat")]
sameSample=intersect(row.names(data1), row.names(cli))
data1=data1[sameSample,]
cli=cli[sameSample,]
data1=cbind(cli,data1)

riskScore=0
for (i in 3:ncol(data1)) {
  riskScore=data1[,i]*model$coef[i-2]+riskScore
}
data1=cbind(data1, riskScore)
data1=data1[row.names(riskRT),]
riskRT=cbind(riskRT, data1[,"riskScore"])
colnames(riskRT)[ncol(riskRT)]=geneFiles


#for(i in geneFiles){
	#读取基因列表
	header=unlist(strsplit(i, "\\."))
	model=read.table(i, header = T, sep = "\t", check.names = F)
	modelgene=model$gene
	samegene=intersect(modelgene,rownames(data))
	data1=data[samegene,]
	colnames(data1)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", colnames(data1))
	data1=t(data1)
	data1=avereps(data1)
	data1=as.data.frame(data1)
	
	#合并生存数据文件
	cli=riskRT[,c("futime", "fustat")]
	sameSample=intersect(row.names(data1), row.names(cli))
	data1=data1[sameSample,]
	cli=cli[sameSample,]
	data1=cbind(cli,data1)
	
	riskScore=0
	for (j in 1:ncol(data1)) {
	  riskScore=data[,j]*model$coef[j]+riskScore
	}
	data1=cbind(data1, riskScore)
	data1=data1[row.names(riskRT),]
	riskRT=cbind(riskRT, data1[,"riskScore"])
	colnames(riskRT)[ncol(riskRT)]=header[[1]]
}

#输出所有模型的风险打分
riskOut=rbind(ID=colnames(riskRT), riskRT)
write.table(riskOut, file="risk.models.txt", sep="\t", col.names=F, quote=F)


#########绘制生存曲线函数#########
bioSurvival=function(inputFile=null, outFile=null, varName=null){
	#读取输入文件
	rt=read.table(inputFile, header=T, sep="\t", check.names=F)
	#比较高低风险组生存差异，得到显著性p值
	rt$Type=ifelse(rt[,varName]>median(rt[,varName]), "high", "low")
	diff=survdiff(Surv(futime, fustat) ~ Type,data = rt)
	pValue=1-pchisq(diff$chisq,df=1)
	if(pValue<0.001){
		pValue="p<0.001"
	}else{
		pValue=paste0("p=",sprintf("%.03f",pValue))
	}
	fit <- survfit(Surv(futime, fustat) ~ Type, data = rt)
		
	#绘制生存曲线
	surPlot=ggsurvplot(fit, 
		           data=rt,
		           conf.int=F,
		           pval=pValue,
		           pval.size=6,
		           title=varName,
		           legend.title="Risk",
		           legend.labs=c("High risk", "Low risk"),
		           xlab="Time(years)",
		           break.time.by = 1,
		           palette=c("#ED0000FF","#00468BFF"),
		           risk.table=F,
		       	   risk.table.title="",
		           risk.table.col = "strata",
		           risk.table.height=.25)
	pdf(file=outFile, onefile = FALSE, width=6, height=5)
	print(surPlot)
	dev.off()
}

#########定义绘制ROC曲线函数#########
bioROC=function(inputFile=null, outFile=null, varName=null){
	#读取输入文件
	rt=read.table(inputFile, header=T, sep="\t", check.names=F)
	#ROC曲线
	ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
	               marker=rt[,varName], cause=1,
	               weighting='aalen',
	               times=c(1,3,5), ROC=TRUE)
	pdf(file=outFile, width=5, height=5)
	plot(ROC_rt,time=1,col="#42B540FF",title=FALSE,lwd=2)
	plot(ROC_rt,time=3,col="#ED0000FF",add=TRUE,title=FALSE,lwd=2)
	plot(ROC_rt,time=5,col="#00468BFF",add=TRUE,title=FALSE,lwd=2)
	text(0.75, 0.24, varName, cex=1.2)
	legend('bottomright',
	        c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
	          paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
	          paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
	        col=c("#42B540FF","#ED0000FF","#00468BFF"),lwd=2,bty = 'n')
	dev.off()
}

for(varName in colnames(riskRT)[3:ncol(riskRT)]){
	bioSurvival(inputFile="risk.models.txt", outFile=paste0("sur.",varName,".pdf"), varName=varName)
	bioROC(inputFile="risk.models.txt", outFile=paste0("ROC.",varName,".pdf"), varName=varName)
}
