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
gene=entrezIDs[entrezIDs!="NA"]        #ȥ������idΪNA�Ļ���
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#GO��������
kk=enrichGO(gene=gene,OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable =T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
#�������������Ľ��
write.table(GO,file="GO.txt",sep="\t",quote=F,row.names = F)

#������ʾTerm��Ŀ
showNum=10
if(nrow(GO)<30){
	showNum=nrow(GO)
}

#��״ͼ
pdf(file="barplot.pdf", width=9, height=7)
bar=barplot(kk, drop = TRUE, showCategory =showNum,split="ONTOLOGY",color = colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()
		
#����ͼ
pdf(file="bubble.pdf", width=9, height=7)
bub=dotplot(kk,showCategory = showNum, orderBy = "GeneRatio",split="ONTOLOGY", color = colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bub)
dev.off()

#��ȡGO��Ϣ
go=data.frame(Category=GO$ONTOLOGY, ID=GO$ID, Term=GO$Description, Genes = gsub("/", ", ", GO$geneID), adj_pval = GO$p.adjust)
#��ȡ�����logFC
genelist <- data.frame(ID = rt$gene, logFC = rt$logFC)
row.names(genelist)=genelist[,1]
#����Ȧͼ����
circ <- circle_dat(go, genelist)
termNum =8         #�޶�GO��Ŀ
termNum=ifelse(nrow(go)<termNum,nrow(go),termNum)
geneNum=200        #�޶�������Ŀ
geneNum=ifelse(nrow(genelist)<geneNum, nrow(genelist), geneNum)
#����Ȧͼ
chord <- chord_dat(circ, genelist[1:geneNum,], go$Term[1:termNum])
pal= pal_lancet(alpha=0.5)(9)
pdf(file="GOcircos.pdf", width=11, height=11)
GOChord(chord, 
        ribbon.col = palette(pal),
        space = 0.021,           #����֮��ļ��
        gene.order = 'logFC',    #����logFCֵ�Ի�������
        gene.space = 0.25,       #��������ԲȦ����Ծ���
        gene.size = 5,           #�����������С 
        border.size = NA,       #������ϸ
        process.label = 6)       #GO�����С
dev.off()


