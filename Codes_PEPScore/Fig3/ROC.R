#install.packages("survival")
#install.packages("survminer")
#install.packages("timeROC")


#引用包
library(survival)
library(survminer)
library(timeROC)
setwd("E:\\My_projects\\LUSC\\08model_2\\01_05")      #设置工作目录

#定义绘制ROC曲线函数
bioROC=function(inputFile=null, rocFile=null){
	#读取输入文件
	rt=read.table(inputFile, header=T, sep="\t", check.names=F)
	#ROC曲线
	ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
	               marker=rt$riskScore, cause=1,
	               weighting='aalen',
	               times=c(1,3,5), ROC=TRUE)
	pdf(file=rocFile,width=5,height=5)
	plot(ROC_rt,time=1,col="#42B540FF",title=FALSE,lwd=2)
	plot(ROC_rt,time=3,col='#00468BFF',add=TRUE,title=FALSE,lwd=2)
	plot(ROC_rt,time=5,col="#ED0000FF",add=TRUE,title=FALSE,lwd=2)
	legend('bottomright',
	        c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
	          paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
	          paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
	        col=c("#42B540FF","#00468BFF","#ED0000FF"),lwd=2,bty = 'n')
	dev.off()
}

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


bioROC(inputFile="risk.TCGA.txt", rocFile="TCGA.ROC.pdf")
bioROC(inputFile="risk.GEO1.txt", rocFile="GEO1.ROC.pdf")
bioROC(inputFile="risk.GEO5.txt", rocFile="GEO5.ROC.pdf")

dev.off()
