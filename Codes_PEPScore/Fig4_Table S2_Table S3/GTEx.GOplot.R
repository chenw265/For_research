
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
gene=entrezIDs[entrezIDs!="NA"]        #去除基因id为NA的基因
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#GO富集分析
kk=enrichGO(gene = gene,OrgDb = org.Hs.eg.db, pvalueCutoff =1, qvalueCutoff = 1, ont="all", readable =T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
#保存富集结果
write.table(GO,file="GO.txt",sep="\t",quote=F,row.names = F)

ego=read.table("GO.txt", header = T,sep="\t",check.names=F)      #读取kegg富集结果文件
go=data.frame(Category = "All",ID = ego$ID,Term = ego$Description, Genes = gsub("/", ", ", ego$geneID), adj_pval = ego$p.adjust)

#读取基因的logFC文件
id.fc <- read.table(inputFile, header = T,sep="\t",check.names=F)
genelist <- data.frame(ID = id.fc$gene, logFC = id.fc$logFC)
row.names(genelist)=genelist[,1]

circ <- circle_dat(go, genelist)
termNum = 5                                     #限定term数目
geneNum = nrow(genelist)                        #限定基因数目

chord <- chord_dat(circ, genelist[1:geneNum,], go$Term[1:termNum])
pal= pal_lancet(alpha=0.5)(9)
pdf(file="GOcirc.pdf",width = 13,height = 11)
GOChord(chord,
        ribbon.col = palette(pal),
        space = 0.021,           #基因之间的间距
        gene.order = 'logFC',    #按照logFC值对基因排序
        gene.space = 0.25,       #基因名跟圆圈的相对距离
        gene.size = 5,           #基因名字体大小 
        border.size = NA,       #线条粗细
        process.label = 6)     #term字体大小
dev.off()

termCol <- c("#223D6C","#D20A13","#FFD121","#088247","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")
pdf(file="GOcluster.pdf",width = 14,height = 10)
GOCluster(circ.gsym, 
          go$Term[1:termNum], 
          lfc.space = 0.2,                   #倍数跟树间的空隙大小
          lfc.width = 1,                     #变化倍数的圆圈宽度
          term.col = termCol[1:termNum],     #自定义term的颜色
          term.space = 0.2,                  #倍数跟term间的空隙大小
          term.width = 1)                    #富集term的圆圈宽度
dev.off()          
