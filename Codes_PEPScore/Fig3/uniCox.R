
#���ð�
library(survival)
library(survminer)

coxPfilter=0.01                   #�����Թ��˱�׼
inputFile="TCGA.expTime.txt"      #�����ļ�
setwd("E:\\My_projects\\LUSC\\08model_2\\01_05")      #���ù���Ŀ¼

#��ȡ�����ļ�
rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
rt$futime=rt$futime/365
rt$fustat[which(rt$futime > 10)] <- 0
rt$futime[which(rt$futime > 10)] <- 10


#�Ի������ѭ�����ҳ�Ԥ����صĻ���
outTab=data.frame()
sigGenes=c("futime","fustat")

for(i in colnames(rt[,3:ncol(rt)])){
  #cox����
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  #����Ԥ����صĻ���
  if(coxP<coxPfilter){
    sigGenes=c(sigGenes,i)
    outTab=rbind(outTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
    }
}
#��������صĽ��
write.table(outTab,file="TCGA.uniCox.txt",sep="\t",row.names=F,quote=F)

#�����������������ı�����
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="TCGA.uniSigExp.txt",sep="\t",row.names=F,quote=F)

library(forestplot)
data <- read.table("TCGA.uniCox.txt",header=T,sep="\t",row.names=1,check.names=F)
data=data.frame(id=rownames(data),data)
data$pvalue=as.numeric(data$pvalue)
data$HR=as.numeric(data$HR)
data$HR.95L=as.numeric(data$HR.95L)
data$HR.95H=as.numeric(data$HR.95H)
data$P.Value1 <- round(data$pvalue,3)
data$P.Value1[which(data$P.Value < 0.001)] <- "<0.001"
hr=round(data$HR,3)
hrLow=round(data$HR.95L,3)
hrHigh=round(data$HR.95H,3)
#hrLow  <- sprintf("%.3f",outTab$HR.95L)
#hrHigh <- sprintf("%.3f",data$"HR.95H")
data$aHR <- paste0(hr,"(",hrLow,"-",hrHigh,")")
tabletext <- cbind(c("Uni-COX",data$id),
                   c("HR (95%  CI)",data$aHR),
                   c("P Value",data$P.Value1))
data1 <- na.omit(data)
data1$Group <- ""
data1$Group[which(data1$HR> 1)] <- "red"
data1$Group[which(data1$HR < 1)] <- "green"

fn <- local({
  i = 0
  no_lines <- sum(!is.na(data$aHR))
  b_clrs = data1$Group
  l_clrs = colorRampPalette(colors=c("black", "black"))(no_lines)
  
  function(..., clr.line, clr.marker){
    i <<- i + 1
    fpDrawNormalCI(..., clr.line = l_clrs[i], clr.marker = b_clrs[i])
  }
})

nrow(data)+2

a<-forestplot(fn.ci_norm = fn,labeltext=tabletext,
              grid = structure(c(0,0.5,1,1.5), 
                               gp = gpar(lty = 2, col = "grey")),
              hrzl_lines = list("1" = gpar(lty=1,lwd=1.5),
                                "2" = gpar(lty=1,lwd=1.5),
                                #"12" = gpar(lty=1,lwd=1.5),
                                #"13" = gpar(lty=1,lwd=1.5),
                                "23" = gpar(lty=1,lwd=1.5)),
              graph.pos=2, #ΪPvalue����ͼ���ڵ�λ��
              mean=c(NA,data$HR),
              lower=c(NA,data$HR.95L),
              upper=c(NA,data$HR.95H),
              xticks = c(0,0.5,1,1.5),
              xlog = F,
              #�������
              # title="Hazard Ratio Plot",
              ##����x��
              #xlab="    ---Favor death",
              #fpTxtGp�����е�cex�������ø�������Ĵ�С
              txt_gp=fpTxtGp(cex = 0.8,
                             ticks=gpar(cex=0.7)),
              xlab  = gpar(fontfamily = "", cex = 1.1),
              ##fpColors����������ɫ
              #is.summary=c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,rep(FALSE,2),rep(TRUE,3)),
              col=fpColors(lines="black", zero = "gray50"),
              graphwidth = unit(2.5, "cm"),
              
              #����ͼ�л�׼�ߵ�λ��
              zero=1,
              cex=0.9, lineheight = "auto",
              colgap=unit(8,"mm"),
              #���Ӵ�С���ߵĿ���
              lwd.ci=1.3, boxsize=0.3,
              lwd.xaxis = 1.3,
              #����ͼ��������С���ߣ��߶�
              ci.vertices=TRUE, ci.vertices.height = 0.2,
              xticks.digits = 0.1)
pdf(file="forest4.pdf", height=6, width=5)
print(a)
dev.off()

















############����ɭ��ͼ����############
bioForest=function(coxFile=null,forestFile=null,forestCol=null){
	#��ȡ�����ļ�
	rt <- read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
	gene <- rownames(rt)
	hr <- sprintf("%.3f",rt$"HR")
	hrLow  <- sprintf("%.3f",rt$"HR.95L")
	hrHigh <- sprintf("%.3f",rt$"HR.95H")
	Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
	pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
		
	#���ͼ��
	height=nrow(rt)/12.5+5
	pdf(file=forestFile, width = 7,height = height)
	n <- nrow(rt)
	nRow <- n+1
	ylim <- c(1,nRow)
	layout(matrix(c(1,2),nc=2),width=c(3,2.5))
		
	#����ɭ��ͼ��ߵ��ٴ���Ϣ
	xlim = c(0,3)
	par(mar=c(4,2.5,2,1))
	plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
	text.cex=0.8
	text(0,n:1,gene,adj=0,cex=text.cex)
	text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
	text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
		
	#����ɭ��ͼ
	par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
	xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
	plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
	arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
	abline(v=1,col="black",lty=2,lwd=2)
	boxcolor = ifelse(as.numeric(hr) > 1, forestCol[1], forestCol[2])
	points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.6)
	axis(1)
	dev.off()
}

#���ú���������ɭ��ͼ
bioForest(coxFile="TCGA.uniCox.txt",forestFile="forest.pdf",forestCol=c("#ED0000FF","#42B540FF"))
