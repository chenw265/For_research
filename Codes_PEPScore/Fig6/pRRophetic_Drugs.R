#安装依赖
#BiocManager::install(c("car", "ridge", "preprocessCore", "genefilter", "sva"))

#将下载的包放在当前工作路径内
#install.packages("pRRophetic_0.5.tar.gz", repos = NULL, dependencies = TRUE)

#Drugs
#681640	(5Z)-7-Oxozeaenol	17-AAG	5-Fluorouracil	A-443654	A-770041
#AC220	Afatinib	AG-014699	AICAR	AKT inhibitor VIII	AMG-706	AP-24534	
#AR-42	AS601245	AS605240	AT-7519	ATRA	AUY922	Axitinib	AZ628	AZD6482	
#AZD7762	AZD8055	BAY 61-3606	Belinostat	Bexarotene	BEZ235	BHG712	BI-2536	
#Bicalutamide	BIRB 0796	BIX02189	Bleomycin	Bleomycin (50 uM)	BMS345541	
#BMS-509744	BMS-536924	BMS-708163	BMS-754807	Bortezomib	Bosutinib	
#Bryostatin 1	BX-795	BX-912	CAL-101	Camptothecin	CAY10603	CCT007093	
#CCT018159	CEP-701	Cetuximab	CGP-082996	CGP-60474	CH5424802	CHIR-99021	
#CI-1040	Cisplatin	CMK	CP466722	CP724714	Crizotinib	CUDC-101	CX-5461	
#Cyclopamine	Cytarabine	Dabrafenib	Dasatinib	DMOG	Docetaxel	Doxorubicin	
#EHT 1864	EKB-569	Elesclomol	Embelin	Epothilone B	Erlotinib	Etoposide	
#EX-527	FH535	FK866	FMK	Foretinib	FR-180204	FTI-277	GDC0449	GDC0941	Gefitinib	
#Gemcitabine	Genentech Cpd 10	GNF-2	GSK1070916	GSK1904529A	GSK2126458	
#GSK269962A	GSK429286A	GSK-650394	GSK690693	GW 441756	GW-2580	GW843682X	
#HG-5-113-01	HG-5-88-01	HG-6-64-1	I-BET-762	Imatinib	IOX2	IPA-3	Ispinesib 
#Mesylate	JNJ-26854165	JNK Inhibitor VIII	JNK-9L	JQ1	JQ12	JW-7-24-1	
#JW-7-52-1	KIN001-055	KIN001-102	KIN001-135	KIN001-236	KIN001-244	
#KIN001-260	KIN001-266	KIN001-270	KU-55933	Lapatinib	LAQ824	Lenalidomide	
#LFM-A13	Linifanib 	Lisitinib	LY317615	Masitinib	Methotrexate	MG-132	
#Midostaurin	Mitomycin C	MK-2206	MLN4924	MP470	MPS-1-IN-1	MS-275	
#Navitoclax	NG-25	Nilotinib	NPK76-II-72-1	NSC-207895	NSC-87877	NU-7441	
#Nutlin-3a (-)	Obatoclax Mesylate	Olaparib	OSI-027	OSI-930	OSU-03012	
#PAC-1	Paclitaxel	Parthenolide	Pazopanib	PD-0325901	PD-0332991	PD-173074	
#PF-4708671	PF-562271	PFI-1	PHA-665752	PHA-793887	Phenformin	PI-103	
#PIK-93	piperlongumine	PLX4720	Pyrimethamine	QL-VIII-58	QL-X-138	QL-XI-92	
#QL-XII-47	QL-XII-61	QS11	Rapamycin	RDEA119	RO-3306	Roscovitine	rTRAIL	
#Ruxolitinib	Salubrinal	Saracatinib	SB 216763	SB 505124	SB52334	SB590885	
#selumetinib	SGC0946	Shikonin	SL 0101-1	SN-38	SNX-2112	Sorafenib	STF-62247	
#S-Trityl-L-cysteine	Sunitinib	T0901317	TAE684	TAK-715	Talazoparib	
#Tamoxifen	Temozolomide	Temsirolimus	TG101348	TGX221	Thapsigargin	
#THZ-2-102-1	THZ-2-49	Tipifarnib	Tivozanib	TL-1-85	TL-2-105	TPCA-1	
#Trametinib	Tubastatin A	TW 37	UNC0638	UNC1215	Veliparib	Vinblastine	
#Vinorelbine	VNLG/124	Vorinostat	VX-11e	VX-680	VX-702	WH-4-023	
#WZ-1-84	WZ3105	XAV939	XL-184	XMD11-85h	XMD13-2	XMD14-99	XMD15-27	
#XMD8-85	XMD8-92	Y-39983	YK 4-279	YM155	YM201636	ZG-10	Zibotentan	
#Z-LLNle-CHO	ZM-447439	ZSTK474


#引用包
library(limma)
library(ggpubr)
library(ggplot2)
library(pRRophetic)
set.seed(17980122)
expFile="symbol.txt"   #表达文件的输入
riskFile="risk.TCGA.txt"   #风险输入文件
setwd("E:\\My_projects\\LUSC\\16DrugsSensitive")  #工作路径设定

#读取文件,并对输入文件进行整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.5,]

#删掉正常样品
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]
data=t(data)
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
data=t(avereps(data))

#预测药物敏感性
drug="Erlotinib"   #药物名称
data(PANCANCER_IC_Tue_Aug_9_15_28_57_2016)
data(cgp2016ExprRma)
DrugNames=unique(drugData2016$Drug.name) #获取数据库所含有药物名称
senstivity=pRRopheticPredict(data, drug, selection = 1, dataset = "cgp2016")
senstivity=senstivity[senstivity!="NaN"]
#senstivity[senstivity>quantile(senstivity,0.99)]=quantile(senstivity,0.99)

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#风险文件和药物敏感性合并
sameSample=intersect(row.names(risk), names(senstivity))
risk=risk[sameSample, "Risk", drop=F]
senstivity=senstivity[sameSample]
rt=cbind(risk, senstivity)

#设置比较组
rt$Risk=factor(rt$Risk, levels = c("low", "high"))
type=levels(factor(rt[,"Risk"]))
comp=combn(type,2)
my_comparisons=list()
for (i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#绘制箱线图
boxplot=ggboxplot(rt, x="Risk", y="senstivity", color="Risk",
                  xlab = "Risk",
                  ylab = paste0(drug, " senstivity (IC50)"),
                  legend.title="Risk",
                  palette = c("#00468BFF", "#ED0000FF"),
                  add = "jitter"
)+
  stat_compare_means(comparisons = my_comparisons)
pdf(file = paste0(drug, ".pdf"), width = 5, height = 5)
print(boxplot)
dev.off()

#多药物循环画图
drugs=c("Erlotinib", "Cetuximab", "Cisplatin", "Etoposide", "Vinblastine", "Gemcitabine", "Docetaxel", "Paclitaxel")
for (i in 1:length(drugs)) {
  drug=drugs[i]
  senstivity=pRRopheticPredict(data, drug, selection = 1, dataset = "cgp2016")
  senstivity=senstivity[senstivity!="NaN"]
  #senstivity[senstivity>quantile(senstivity,0.99)]=quantile(senstivity,0.99)
  
  #读取风险文件
  risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
  
  #风险文件和药物敏感性合并
  sameSample=intersect(row.names(risk), names(senstivity))
  risk=risk[sameSample, "Risk", drop=F]
  senstivity=senstivity[sameSample]
  rt=cbind(risk, senstivity)
  
  #设置比较组
  rt$Risk=factor(rt$Risk, levels = c("low", "high"))
  type=levels(factor(rt[,"Risk"]))
  comp=combn(type,2)
  my_comparisons=list()
  for (i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
  
  #绘制箱线图
  boxplot=ggboxplot(rt, x="Risk", y="senstivity", fill="Risk",
                    xlab = "Risk",
                    ylab = paste0(drug, " senstivity (IC50)"),
                    legend.title="Risk",
                    palette = c("#00468BFF", "#ED0000FF"),
                    size = 0.5,
                    outlier.shape = NA,
                    notch = T,
                    add = "jitter",
                    add.params = list(size=1.5, color="Risk"
                    )
  )+
    stat_compare_means(comparisons = my_comparisons)
  pdf(file = paste0(drug, ".pdf"), width = 3, height = 3.5)
  print(boxplot)
  dev.off()
}




















#筛选敏感性药物
drugs=c("Erlotinib", "Cetuximab")

system.time({ 
  cl <- makeCluster(3)  
  results <- parLapply(cl,drugs,
                       function(x){
                         library(pRRophetic) 
                         predictedPtype=pRRopheticPredict(
                           testMatrix=data,
                           drug=x,
                           tissueType = "all", 
                           batchCorrect = "eb",
                           selection=1,
                           dataset = "cgp2016")
                         return(predictedPtype)
                       }) # lapply的并行版本
  res.df <- do.call('rbind',results) # 整合结果
  stopCluster(cl) # 关闭集群
})





ggboxplot(rt, x="Risk", y="senstivity", color="Risk",
          xlab = "Risk",
          ylab = paste0(drug, " senstivity (IC50)"),
          legend.title="Risk",
          palette =c("#00468BFF", "#ED0000FF") ,
          add = "jitter"
)+
  stat_compare_means(comparisons = my_comparisons)


c("#3C5488FF", "#E64B35FF")
c("#0099B4FF", "#AD002AFF")