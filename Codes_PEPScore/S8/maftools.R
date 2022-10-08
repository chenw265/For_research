
#if (!require("BiocManager"))
#    install.packages("BiocManager")
#BiocManager::install("maftools")


library(maftools)       #���ð�
setwd("E:\\My_projects\\LUSC\\11Maftools")      #���ù���Ŀ¼

#��ȡģ���ļ�����ȡ�����б�
geneRT=read.table("multiCox.txt", header=T, sep="\t", check.names=F, row.names=1)
gene=row.names(geneRT)

#��ȡ�ٴ������ļ�
clinical=read.table("clinical.txt", header=T, sep="\t", check.names=F, row.names=1)
outTab=rbind(Tumor_Sample_Barcode=colnames(clinical), clinical)
write.table(outTab, file="ann.txt", sep="\t", quote=F, col.names=F)

#�����ٲ�ͼ
pdf(file="oncoplot.pdf", width=8, height=7)
maf=read.maf(maf="new.maf", clinicalData="ann.txt")
oncoplot(maf=maf, clinicalFeatures=c("Age","Gender","Stage"), genes=gene)
dev.off()