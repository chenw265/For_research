
#引用包
library(survival)
library(survminer)
library(limma)
coxPfilter=0.99                 #显著性过滤标准
cliFile="tcgatime.txt"
inputFile="tcga.pyroptosisExp.txt"      #输入文件
setwd("E:\\My_projects\\LUSC\\01 PyroptosisExp")      #设置工作目录

#读取TCGA基因表达文件,并对数据进行处理
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
tcga=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
tcga=avereps(tcga)
tcga=log2(tcga+1)

#删掉正常样品
group=sapply(strsplit(colnames(tcga),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
tcga=tcga[,group==0]
tcga=t(tcga)
rownames(tcga)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(tcga))
data=avereps(tcga)

#读取生存数据
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

#数据合并并输出结果
sameSample=intersect(row.names(data), row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]
out=cbind(cli, data)
out=cbind(id=row.names(out), out)
write.table(out, file="pytoptosisTCGA.expTime.txt", sep="\t", row.names=F, quote=F)

#读取输入文件
rt=out

#对基因进行循环，找出预后相关的基因
outTab=data.frame()
sigGenes=c("futime","fustat")

for(i in colnames(rt[,3:ncol(rt)])){
  #cox分析
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  #保留预后相关的基因
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
#输出单因素的结果
write.table(outTab,file="pyroptosisTCGA.uniCox.txt",sep="\t",row.names=F,quote=F)

#输出单因素显著基因的表达量
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="pytoptosisTCGA.uniSigExp.txt",sep="\t",row.names=F,quote=F)


library(forestplot)
data = outTab
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

a<-forestplot(fn.ci_norm = fn,labeltext=tabletext,
              grid = structure(c(0.25,0.5,1,2,4), 
                               gp = gpar(lty = 2, col = "grey")),
              hrzl_lines = list("1" = gpar(lty=1,lwd=1.5),
                                "2" = gpar(lty=1,lwd=1.5),
                                #"12" = gpar(lty=1,lwd=1.5),
                                #"13" = gpar(lty=1,lwd=1.5),
                                "53" = gpar(lty=1,lwd=1.5)),
              graph.pos=2, #为Pvalue箱线图所在的位置
              mean=c(NA,data$HR),
              lower=c(NA,data$HR.95L),
              upper=c(NA,data$HR.95H),
              xticks = c(0.25,0.5,1.000,2,4),
              xlog = T,
              #定义标题
              # title="Hazard Ratio Plot",
              ##定义x轴
              #xlab="    ---Favor death",
              #fpTxtGp函数中的cex参数设置各个组件的大小
              txt_gp=fpTxtGp(cex = 0.8,
                             ticks=gpar(cex=0.7)),
              xlab  = gpar(fontfamily = "", cex = 1.1),
              ##fpColors函数设置颜色
              #is.summary=c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,rep(FALSE,2),rep(TRUE,3)),
              col=fpColors(lines="black", zero = "gray50"),
              graphwidth = unit(2.5, "cm"),
              
              #箱线图中基准线的位置
              zero=1,
              cex=0.9, lineheight = "auto",
              colgap=unit(8,"mm"),
              #箱子大小，线的宽度
              lwd.ci=1.3, boxsize=0.3,
              lwd.xaxis = 1.3,
              #箱线图两端添加小竖线，高度
              ci.vertices=TRUE, ci.vertices.height = 0.2,
              xticks.digits = 0.1)
pdf(file="forest4.pdf", height=8, width=8)
print(a)
dev.off()

