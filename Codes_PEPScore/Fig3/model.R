
#���ð�
library(survival)

trainFile="TCGA.uniSigExp.txt"      #train�������ļ�
testFile="GEO5.expTime.txt"          #test�������ļ�

setwd("E:\\My_projects\\LUSC\\08model_2\\01_05")                        #���ù���Ŀ¼
rt=read.table(trainFile, header=T, sep="\t", row.names=1,check.names=F)    #��ȡtrain�������ļ�

#COXģ�͹���
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)

#���ģ�͵Ĺ�ʽ
outMultiTab=data.frame()
outMultiTab=cbind(
				coef=multiCoxSum$coefficients[,"coef"],
				HR=multiCoxSum$conf.int[,"exp(coef)"],
				HR.95L=multiCoxSum$conf.int[,"lower .95"],
				HR.95H=multiCoxSum$conf.int[,"upper .95"],
				pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outMultiTab=cbind(id=row.names(outMultiTab),outMultiTab)
outMultiTab=outMultiTab[,1:2]
write.table(outMultiTab,file="multiCox.txt",sep="\t",row.names=F,quote=F)

#���train��ķ����ļ�(TCGA)
trainScore=predict(multiCox,type="risk",newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
Risk=as.vector(ifelse(trainScore>median(trainScore),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(trainScore), Risk)
write.table(cbind(id=rownames(outTab),outTab),file="risk.TCGA.txt",sep="\t",quote=F,row.names=F)

#���test��ķ����ļ�(GEO)
rt=read.table(testFile, header=T, sep="\t", row.names=1, check.names=F)
rt$futime=rt$futime/365
rt$fustat[which(rt$futime > 8)] <- 0
rt$futime[which(rt$futime > 8)] <- 8

testFinalGeneExp=rt[,coxGene]
testScore=predict(multiCox,type="risk",newdata=rt)
outCol=c("futime","fustat",coxGene)
Risk=as.vector(ifelse(testScore>median(trainScore),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(testScore), Risk)
write.table(cbind(id=rownames(outTab),outTab),file="risk.GEO5.txt",sep="\t",quote=F,row.names=F)


#���test��ķ����ļ�(GEO1)
rt=read.table("GEO1.expTime.txt", header=T, sep="\t", row.names=1, check.names=F)
rt$futime=rt$futime/365
rt$fustat[which(rt$futime > 8)] <- 0
rt$futime[which(rt$futime > 8)] <- 8

testFinalGeneExp=rt[,coxGene]
testScore=predict(multiCox,type="risk",newdata=rt)
outCol=c("futime","fustat",coxGene)
Risk=as.vector(ifelse(testScore>median(trainScore),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(testScore), Risk)
write.table(cbind(id=rownames(outTab),outTab),file="risk.GEO1.txt",sep="\t",quote=F,row.names=F)



#���test��ķ����ļ�(GEO3)
rt=read.table("GEO3.expTime.txt", header=T, sep="\t", row.names=1, check.names=F)
rt$futime=rt$futime/365
rt$fustat[which(rt$futime > 5)] <- 0
rt$futime[which(rt$futime > 5)] <- 5

testFinalGeneExp=rt[,coxGene]
testScore=predict(multiCox,type="risk",newdata=rt)
outCol=c("futime","fustat",coxGene)
Risk=as.vector(ifelse(testScore>median(trainScore),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(testScore), Risk)
write.table(cbind(id=rownames(outTab),outTab),file="risk.GEO3.txt",sep="\t",quote=F,row.names=F)



#���test��ķ����ļ�(GEO2)
rt=read.table("GEO2.expTime.txt", header=T, sep="\t", row.names=1, check.names=F)
rt$futime=rt$futime/365
rt$fustat[which(rt$futime > 5)] <- 0
rt$futime[which(rt$futime > 5)] <- 5

testFinalGeneExp=rt[,coxGene]
testScore=predict(multiCox,type="risk",newdata=rt)
outCol=c("futime","fustat",coxGene)
Risk=as.vector(ifelse(testScore>median(trainScore),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(testScore), Risk)
write.table(cbind(id=rownames(outTab),outTab),file="risk.GEO2.txt",sep="\t",quote=F,row.names=F)

#���test��ķ����ļ�(GEO4)
rt=read.table("GEO2.expTime.txt", header=T, sep="\t", row.names=1, check.names=F)
rt$futime=rt$futime/365
rt$fustat[which(rt$futime > 5)] <- 0
rt$futime[which(rt$futime > 5)] <- 5

testFinalGeneExp=rt[,coxGene]
testScore=predict(multiCox,type="risk",newdata=rt)
outCol=c("futime","fustat",coxGene)
Risk=as.vector(ifelse(testScore>median(trainScore),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(testScore), Risk)
write.table(cbind(id=rownames(outTab),outTab),file="risk.GEO4.txt",sep="\t",quote=F,row.names=F)
