
#install.packages("survival")
#install.packages("survminer")
#install.packages("timeROC")
#install.packages("ggDCA")


#引用包
library(survival)
library(survminer)
library(timeROC)
library(ggDCA)

predictTime=1       #预测时间
riskFile="nomoRisk.txt"        #列线图风险输入文件
cliFile="clinical.txt"         #临床数据文件
setwd("E:\\My_projects\\LUSC\\15nomo\\DCA")     #修改工作目录

#读取风险输入文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#读取临床数据文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

#合并数据
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli1=cli[samSample,,drop=F]
data=cbind(risk1, cli1)

#DCA分析
rt=cbind(risk1[,c("futime","fustat","Risk","Nomogram")], cli1)
rt[,"Age"]=ifelse(rt[,"Age"]>65, 1, 0)
rt[,"Nomogram"]=ifelse(rt[,"Nomogram"]>median(rt[,"Nomogram"]), 1, 0)
Nomogram<-coxph(Surv(futime,fustat)~Nomogram,rt)
Risk<-coxph(Surv(futime,fustat)~Risk,rt)
Age<-coxph(Surv(futime,fustat)~Age,rt)
#Gender<-coxph(Surv(futime,fustat)~Gender,rt)
#Grade<-coxph(Surv(futime,fustat)~Grade,rt)
#Stage<-coxph(Surv(futime,fustat)~Stage,rt)
T<-coxph(Surv(futime,fustat)~T,rt)
M<-coxph(Surv(futime,fustat)~M,rt)

#绘制决策曲线
pdf(file="DCA.pdf", width=6.5, height=5.2)
d_train=dca(Nomogram,Risk,Age,T,M, times=predictTime)
ggplot(d_train, linetype=1)
dev.off()


######绘制临床的ROC曲线######
rt=cbind(risk1[,c("futime","fustat","riskScore","Nomogram")], cli1)
aucText=c()
#bioCol=rainbow(ncol(rt)-1, s=0.9, v=0.9)
bioCol=c("#00468BFF", "#ED0000FF","#42B540FF","#0099B4FF","#925E9FFF","#DF8F44FF")
pdf(file="cliROC.pdf", width=6, height=6)
#绘制风险得分的ROC曲线
i=3
ROC_rt=timeROC(T=risk$futime,
               delta=risk$fustat,
               marker=risk$riskScore, cause=1,
               weighting='aalen',
               times=c(predictTime),ROC=TRUE)
plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=2)
aucText=c(paste0("Risk", ", AUC=", sprintf("%.3f",ROC_rt$AUC[2])))
abline(0,1)

#对临床数据进行循环，绘制临床数据的ROC曲线
for(i in 4:ncol(rt)){
	ROC_rt=timeROC(T=rt$futime,
				   delta=rt$fustat,
				   marker=rt[,i], cause=1,
				   weighting='aalen',
				   times=c(predictTime),ROC=TRUE)
	plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=2, add=TRUE)
	aucText=c(aucText, paste0(colnames(rt)[i],", AUC=",sprintf("%.3f",ROC_rt$AUC[2])))
}
#绘制图例，得到ROC曲线下的面积
legend("bottomright", aucText,lwd=2,bty="n",col=bioCol[1:(ncol(rt)-1)])
dev.off()
