
#install.packages("digest")
#install.packages("GOplot")
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library("ggsci")
library("scales")
library(GOplot)
setwd("E:\\My_projects\\LUSC\\20RiskPyroptosis\\GOcir_KEGGcir")              #���ù���Ŀ¼
#���ù���Ŀ¼
inputFile="diff.txt"
pvalueFilter=0.05           #pֵ��������
qvalueFilter=0.05           #�������pֵ��������
rt=read.table(inputFile,sep="\t",header=T,check.names=F)        #��ȡ�����ļ�

#��������ת��Ϊ����id
genes=as.vector(rt[,1])
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=cbind(rt,entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        #ȥ������idΪNA�Ļ���
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#kegg��������
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =1, qvalueCutoff =1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$gene[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
#���渻�����
write.table(KEGG,file="KEGG.txt",sep="\t",quote=F,row.names = F)





ego=read.table("KEGG.txt", header = T,sep="\t",check.names=F)      #��ȡkegg��������ļ�
go=data.frame(Category = "All",ID = ego$ID,Term = ego$Description, Genes = gsub("/", ", ", ego$geneID), adj_pval = ego$p.adjust)

#��ȡ�����logFC�ļ�
id.fc <- read.table(inputFile, header = T,sep="\t",check.names=F)
genelist <- data.frame(ID = id.fc$gene, logFC = id.fc$logFC)
row.names(genelist)=genelist[,1]

circ <- circle_dat(go, genelist)
termNum = 3                                     #�޶�term��Ŀ
geneNum = nrow(genelist)                        #�޶�������Ŀ

chord <- chord_dat(circ, genelist[1:geneNum,], go$Term[1:termNum])
pdf(file="circ.pdf",width = 9,height = 8)
GOChord(chord, 
        space = 0.021,           #����֮��ļ��
        gene.order = 'logFC',    #����logFCֵ�Ի�������
        gene.space = 0.25,       #��������ԲȦ����Ծ���
        gene.size = 4,           #�����������С 
        border.size =NA,       #������ϸ
        process.label = 9)     #term�����С
dev.off()

termCol <- c("#223D6C","#D20A13","#FFD121","#088247","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")
pdf(file="KEGGcluster.pdf",width = 14,height = 10)
GOCluster(circ.gsym, 
          go$Term[1:termNum], 
          lfc.space = 0.2,                   #����������Ŀ�϶��С
          lfc.width = 1,                     #�仯������ԲȦ����
          term.col = termCol[1:termNum],     #�Զ���term����ɫ
          term.space = 0.2,                  #������term��Ŀ�϶��С
          term.width = 1)                    #����term��ԲȦ����
dev.off()          