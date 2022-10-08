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


#���ð�
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library("GOplot")
library("ggsci")
library("scales")

pvalueFilter=0.05       #pֵ��������
qvalueFilter=0.05       #�������pֵ��������
diffFile="diff.txt"
#������ɫ
colorSel="qvalue"
if(qvalueFilter>0.05){
        colorSel="pvalue"
}

setwd("E:\\My_projects\\LUSC\\05WGCNAinDiff\\GO")                            #���ù���Ŀ¼

#��ȡģ���еĻ���
geneVec=c()
geneFiles=dir()       #��ȡĿ¼�������ļ�
geneFiles=grep("_genes.txt", geneFiles, value=T)     #��ȡ_genes.txt��β���ļ�
for(geneFile in geneFiles){
        gene=read.table(geneFile, header=F, sep="\t", check.names=F)
        geneVec=c(geneVec, as.vector(gene[,1]))
}

#rt=read.table("TCGA.immuneDiff.txt", header=T, sep="\t", check.names=F)     #��ȡ�����ļ�
rt=as.data.frame(geneVec)     #��ȡ�����ļ�
rownames(rt)=rt$geneVec
diff=read.table(diffFile, header=T, sep="\t", check.names=F)
rownames(diff)=diff$gene
gene=intersect(rt$geneVec, diff$gene)
rt=diff[gene,]

#��������ת��Ϊ����id
genes=as.vector(rt[,1])
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=cbind(rt,entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        #ȥ������idΪNA�Ļ���
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#kegg��������
kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$gene[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
#�������������Ľ��
write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)

#������ʾͨ·����Ŀ
showNum=30
if(nrow(KEGG)<showNum){
	showNum=nrow(KEGG)
}

#��״ͼ
pdf(file="barplot.pdf", width=8, height=7)
barplot(kk, drop = TRUE, showCategory = showNum, color = colorSel)
dev.off()

#����ͼ
pdf(file="bubble.pdf", width=8, height=7)
dotplot(kk, showCategory = showNum, orderBy = "GeneRatio",color = colorSel)
dev.off()

#��ȡKEGG��Ϣ
kegg=data.frame(Category="ALL", ID = KEGG$ID, Term=KEGG$Description, Genes = gsub("/", ", ", KEGG$geneID), adj_pval = KEGG$p.adjust)
#��ȡ����Ĳ������
genelist <- data.frame(ID = rt$gene, logFC = rt$logFC)
row.names(genelist)=genelist[,1]
#����Ȧͼ����
circ <- circle_dat(kegg, genelist)
termNum =8         #�޶�ͨ·����Ŀ
termNum=ifelse(nrow(kegg)<termNum,nrow(kegg),termNum)
geneNum=200        #�޶��������Ŀ
geneNum=ifelse(nrow(genelist)<geneNum, nrow(genelist), geneNum)
#����Ȧͼ
pal=pal_lancet(alpha = 0.5)(8)
chord <- chord_dat(circ, genelist[1:geneNum,], kegg$Term[1:termNum])
pdf(file="KEGGcircos.pdf", width=11, height=11)
GOChord(chord, 
        ribbon.col = palette(pal),
        space = 0.021,           #����֮��ļ��
        gene.order = 'logFC',    #����logFCֵ�Ի�������
        gene.space = 0.25,       #��������ԲȦ����Ծ���
        gene.size = 5,           #�����������С 
        border.size = NA,       #������ϸ
        process.label = 6)       #term�����С
dev.off()