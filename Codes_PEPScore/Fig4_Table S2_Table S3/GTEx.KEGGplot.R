
#install.packages("digest")
#install.packages("GOplot")
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library("ggsci")
library("scales")
library(GOplot)
setwd("E:\\My_projects\\LUSC\\20RiskPyroptosis\\GOcir_KEGGcir")              #设置工作目录
#设置工作目录
inputFile="diff.txt"
pvalueFilter=0.05           #p值过滤条件
qvalueFilter=0.05           #矫正后的p值过滤条件
rt=read.table(inputFile,sep="\t",header=T,check.names=F)        #读取输入文件

#基因名字转换为基因id
genes=as.vector(rt[,1])
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=cbind(rt,entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        #去除基因id为NA的基因
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#kegg富集分析
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =1, qvalueCutoff =1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$gene[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
#保存富集结果
write.table(KEGG,file="KEGG.txt",sep="\t",quote=F,row.names = F)





ego=read.table("KEGG.txt", header = T,sep="\t",check.names=F)      #读取kegg富集结果文件
go=data.frame(Category = "All",ID = ego$ID,Term = ego$Description, Genes = gsub("/", ", ", ego$geneID), adj_pval = ego$p.adjust)

#读取基因的logFC文件
id.fc <- read.table(inputFile, header = T,sep="\t",check.names=F)
genelist <- data.frame(ID = id.fc$gene, logFC = id.fc$logFC)
row.names(genelist)=genelist[,1]

circ <- circle_dat(go, genelist)
termNum = 3                                     #限定term数目
geneNum = nrow(genelist)                        #限定基因数目

chord <- chord_dat(circ, genelist[1:geneNum,], go$Term[1:termNum])
pdf(file="circ.pdf",width = 9,height = 8)
GOChord(chord, 
        space = 0.021,           #基因之间的间距
        gene.order = 'logFC',    #按照logFC值对基因排序
        gene.space = 0.25,       #基因名跟圆圈的相对距离
        gene.size = 4,           #基因名字体大小 
        border.size =NA,       #线条粗细
        process.label = 9)     #term字体大小
dev.off()

termCol <- c("#223D6C","#D20A13","#FFD121","#088247","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")
pdf(file="KEGGcluster.pdf",width = 14,height = 10)
GOCluster(circ.gsym, 
          go$Term[1:termNum], 
          lfc.space = 0.2,                   #倍数跟树间的空隙大小
          lfc.width = 1,                     #变化倍数的圆圈宽度
          term.col = termCol[1:termNum],     #自定义term的颜色
          term.space = 0.2,                  #倍数跟term间的空隙大小
          term.width = 1)                    #富集term的圆圈宽度
dev.off()          
