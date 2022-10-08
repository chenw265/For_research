#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")
#install.packages("digest")
#install.packages("GOplot")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")


#引用包
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library("GOplot")
library("ggsci")
library("scales")

pvalueFilter=0.05       #p值过滤条件
qvalueFilter=0.05       #矫正后的p值过滤条件
diffFile="diff.txt"
#定义颜色
colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}
	
setwd("E:\\My_projects\\LUSC\\05WGCNAinDiff\\GO")                            #设置工作目录

#读取模块中的基因
geneVec=c()
geneFiles=dir()       #获取目录下所有文件
geneFiles=grep("_genes.txt", geneFiles, value=T)     #提取_genes.txt结尾的文件
for(geneFile in geneFiles){
        gene=read.table(geneFile, header=F, sep="\t", check.names=F)
        geneVec=c(geneVec, as.vector(gene[,1]))
}

#rt=read.table("TCGA.immuneDiff.txt", header=T, sep="\t", check.names=F)     #读取输入文件
rt=as.data.frame(geneVec)     #读取输入文件
rownames(rt)=rt$geneVec
diff=read.table(diffFile, header=T, sep="\t", check.names=F)
rownames(diff)=diff$gene
gene=intersect(rt$geneVec, diff$gene)
rt=diff[gene,]

#基因名字转换为基因id
genes=as.vector(rt[,1])
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        #去除基因id为NA的基因
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#GO富集分析
kk=enrichGO(gene=gene,OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable =T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
#保存显著富集的结果
write.table(GO,file="GO.txt",sep="\t",quote=F,row.names = F)

#定义显示Term数目
showNum=10
if(nrow(GO)<30){
	showNum=nrow(GO)
}

#柱状图
pdf(file="barplot.pdf", width=9, height=7)
bar=barplot(kk, drop = TRUE, showCategory =showNum,split="ONTOLOGY",color = colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()
		
#气泡图
pdf(file="bubble.pdf", width=9, height=7)
bub=dotplot(kk,showCategory = showNum, orderBy = "GeneRatio",split="ONTOLOGY", color = colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bub)
dev.off()

#获取GO信息
go=data.frame(Category=GO$ONTOLOGY, ID=GO$ID, Term=GO$Description, Genes = gsub("/", ", ", GO$geneID), adj_pval = GO$p.adjust)
#读取基因的logFC
genelist <- data.frame(ID = rt$gene, logFC = rt$logFC)
row.names(genelist)=genelist[,1]
#设置圈图参数
circ <- circle_dat(go, genelist)
termNum =8         #限定GO数目
termNum=ifelse(nrow(go)<termNum,nrow(go),termNum)
geneNum=200        #限定基因数目
geneNum=ifelse(nrow(genelist)<geneNum, nrow(genelist), geneNum)
#绘制圈图
chord <- chord_dat(circ, genelist[1:geneNum,], go$Term[1:termNum])
pal= pal_lancet(alpha=0.5)(9)
pdf(file="GOcircos.pdf", width=11, height=11)
GOChord(chord, 
        ribbon.col = palette(pal),
        space = 0.021,           #基因之间的间距
        gene.order = 'logFC',    #按照logFC值对基因排序
        gene.space = 0.25,       #基因名跟圆圈的相对距离
        gene.size = 5,           #基因名字体大小 
        border.size = NA,       #线条粗细
        process.label = 6)       #GO字体大小
dev.off()



