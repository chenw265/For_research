
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install(c("GO.db", "preprocessCore", "impute", "limma"))

#install.packages(c("matrixStats", "Hmisc", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival")) 
#install.packages("WGCNA")


#引用包
library("limma")
library("WGCNA")
expFile="diffGeneExp.txt"      #差异基因的表达文件
traitFile="cluster.txt"
setwd("E:\\My_projects\\LUSC\\05WGCNAinDiff")       #设置工作目录

#读取文件,并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
#data=log2(data+1)
data=data[apply(data,1,sd)>0.01,]
datExpr0=t(data)

##正常和肿瘤数目
#group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
#group=sapply(strsplit(group,""), "[", 1)
#group=gsub("2", "1", group)
#conNum=length(group[group==1])       #正常组样品数目
#treatNum=length(group[group==0])     #肿瘤组样品数目

##性状及数目
traitData=read.table(traitFile, header=T, sep="\t", check.names=F)
traitData$C1 <- 0
traitData$C1[which(traitData$Cluster == "C1")] <- 1
traitData$C2 <- 0
traitData$C2[which(traitData$Cluster == "C2")] <- 1
rownames(traitData)=traitData[,1]
traitData=traitData[,3:ncol(traitData)]
conNum=nrow(traitData[which(traitData$C1 == 1),])       #C1组样品数
treatNum=nrow(traitData[which(traitData$C2 == 1),])    #C2组样品数

###检查缺失值
gsg = goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK){
	# Optionally, print the gene and sample names that were removed:
	if (sum(!gsg$goodGenes)>0)
	  	printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
	if (sum(!gsg$goodSamples)>0)
	  	printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
	# Remove the offending genes and samples from the data:
	datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

###样品聚类
sampleTree = hclust(dist(datExpr0), method = "average")
pdf(file = "1_sample_cluster.pdf", width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
###剪切线
abline(h = 20000, col = "red")
dev.off()

###删除剪切线以下的样品
clust = cutreeStatic(sampleTree, cutHeight = 20000, minSize = 10)
table(clust)
keepSamples = (clust==1)
datExpr0 = datExpr0[keepSamples, ]


###准备临床数据
#traitData=data.frame(Normal=c(rep(1,conNum),rep(0,treatNum)),
#                     Tumor=c(rep(0,conNum),rep(1,treatNum)))
#row.names(traitData)=colnames(data)
#fpkmSamples = rownames(datExpr0)
#traitSamples =rownames(traitData)
#sameSample=intersect(fpkmSamples,traitSamples)
#datExpr0=datExpr0[sameSample,]
#datTraits=traitData[sameSample,]

fpkmSamples = rownames(datExpr0)
traitSamples =rownames(traitData)
sameSample=intersect(fpkmSamples,traitSamples)
datExpr0=datExpr0[sameSample,]
datTraits=traitData[sameSample,]

###样品聚类
sampleTree2 = hclust(dist(datExpr0), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE)
pdf(file="2_sample_heatmap.pdf",width=15,height=12)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()

###power值散点图
enableWGCNAThreads()   #多线程工作
powers = c(1:20)       #幂指数范围1:20
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
pdf(file="3_scale_independence.pdf",width=9,height=5)
par(mfrow = c(1,2))
cex1 = 0.9
###拟合指数与power值散点图
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red") #可以修改
###平均连通性与power值散点图
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


###邻接矩阵转换
sft #查看最佳power值
softPower =sft$powerEstimate #最佳power值
adjacency = adjacency(datExpr0, power = softPower)
softPower

###TOM矩阵
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

###基因聚类
geneTree = hclust(as.dist(dissTOM), method = "average");
pdf(file="4_gene_clustering.pdf",width=12,height=9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()


###动态剪切模块识别
minModuleSize = 25      #模块基因数目
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
pdf(file="5_Dynamic_Tree.pdf",width=8,height=6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()


###相似模块聚类
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average")
pdf(file="6_Clustering_module.pdf",width=7,height=6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.25 #剪切高度可修改
abline(h=MEDissThres, col = "red")
dev.off()


###相似模块合并
merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
pdf(file="7_merged_dynamic.pdf", width = 9, height = 6)
plotDendroAndColors(geneTree, mergedColors,"Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()
moduleColors = mergedColors
table(moduleColors)
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs


###模块与性状数据热图
nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
pdf(file="8_Module_trait.pdf",width=5.5,height=5.5)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(5, 10, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

nTop = 50
###输出每个模块的基因
for (mod in 1:nrow(table(moduleColors))){  
	modules = names(table(moduleColors))[mod]
	probes = colnames(datExpr0)
	inModule = (moduleColors == modules)
	modGenes = probes[inModule]
	write.table(modGenes, file =paste0("9_",modules,"_genes.txt"),sep="\t",row.names=F,col.names=F,quote=F)
}

###输出每个模块的网络输入文件
for(mod in 1:nrow(table(moduleColors))){
	modules = names(table(moduleColors))[mod]
	probes = colnames(datExpr0)
	inModule = (moduleColors == modules)
	modProbes = probes[inModule]
	modGenes = modProbes
	modTOM = TOM[inModule, inModule]
	dimnames(modTOM) = list(modProbes, modProbes)
	outEdge = paste0("10_",modules , "_network.txt")
	outNode = paste0("10_",modules, "_nodes.txt")
	cyt = exportNetworkToCytoscape(modTOM,
		edgeFile = outEdge,
		nodeFile = outNode,
		weighted = TRUE,
		threshold = 0.3,
		nodeNames = modProbes,
		altNodeNames = modGenes,
		nodeAttr = moduleColors[inModule])
}

