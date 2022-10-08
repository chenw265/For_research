library(survival)       #引用包
library(forestplot)
setwd("E:\\My_projects\\LUSC\\09indep")     #设置工作目录
riskFile="risk.TCGA.txt"
cliFile="clinical.txt"
uniOutFile="uniCox.txt"
multiOutFile="multiCox.txt"
#uniForest="uniForest.pdf"
#multiForest="multiForest.pdf"


risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)    #读取风险文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)      #读取临床文件

#数据合并
sameSample=intersect(row.names(cli),row.names(risk))
risk=risk[sameSample,]
cli=cli[sameSample,]
rt=cbind(futime=risk[,1], fustat=risk[,2], cli, riskScore=risk[,(ncol(risk)-1)])

#单因素独立预后分析
uniTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  uniTab=rbind(uniTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
write.table(uniTab,file=uniOutFile,sep="\t",row.names=F,quote=F)
#bioForest(coxFile=uniOutFile, forestFile=uniForest, forestCol="#42B540FF")

#多因素独立预后分析
uniTab1=uniTab
uniTab=uniTab[as.numeric(uniTab[,"pvalue"])<0.05,]
rt1=rt[,c("futime", "fustat", as.vector(uniTab[,"id"]))]
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt1)
multiCoxSum=summary(multiCox)
multiTab=data.frame()
multiTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
multiTab=cbind(id=row.names(multiTab),multiTab)
write.table(multiTab,file=multiOutFile,sep="\t",row.names=F,quote=F)
#bioForest(coxFile=multiOutFile, forestFile=multiForest, forestCol="#ED0000FF")



a=vector(length = ncol(uniTab1))
a[1]="Muti-Cox"
a[-1]=""
data = rbind(uniTab1,a,multiTab)

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
#tabletext2 <- cbind(c("Muti-Cox", data$id[8:11]),
 #                   c("",data$aHR[8:11]),
  #                  c("",data$P.Value1[8:11]))
#tabletext=rbind(tabletext1, tabletext2)
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
              grid = structure(c(0.0,0.5,1.0,1.5,2.0,2.5,3.0), 
                               gp = gpar(lty = 2, col = "grey")),
              hrzl_lines = list("1" = gpar(lty=1,lwd=1.5),
                                "2" = gpar(lty=1,lwd=1.5),
                                "9" = gpar(lty=1,lwd=1.5),
                                "10" = gpar(lty=1,lwd=1.5),
                                "14" = gpar(lty=1,lwd=1.5)),
              graph.pos=2, #为Pvalue箱线图所在的位置
              mean=c(NA,data$HR),
              lower=c(NA,data$HR.95L),
              upper=c(NA,data$HR.95H),
              xticks = c(0.0,0.5,1.0,1.5,2.0,2.5,3.0),
              xlog = F,
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
pdf(file="indepforest4.pdf", height=4.5, width=5)
print(a)
dev.off()





