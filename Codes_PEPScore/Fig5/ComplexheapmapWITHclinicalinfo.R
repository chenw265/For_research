#引用包
library(limma)
library(reshape2)
library(ggpubr)
library(tidyr)

riskFile="risk.TCGA.txt"            #风险文件
cliFile="clinical.txt"
immFile="CIBERSORT-Results.txt"     #免疫细胞浸润结果文件
pFilter=0.05                        #免疫细胞浸润结果的过滤条件
setwd("E:\\My_projects\\LUSC\\12ImmuneRelated\\02ImmuneCor")      #设置工作目录

#读取免疫细胞结果文件，并对数据进行整理
immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)
#immune=immune[immune[,"P-value"]<pFilter,]
data=as.matrix(immune[,1:(ncol(immune)-3)])

#删除正常样品
group=sapply(strsplit(row.names(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[group==0,]
row.names(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(data))
data=avereps(data)

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(data), row.names(risk))
rt=cbind(data[sameSample,,drop=F], risk[sameSample,"Risk",drop=F])
rt=rt[order(rt$Risk, decreasing=T),]
conNum=nrow(rt[rt$Risk=="low",])
treatNum=nrow(rt[rt$Risk=="high",])
group=data.frame(risk$Risk, row.names = rownames(risk))

#数据加工
cn=colnames(rt)
rt1=data.frame(ID=rownames(rt), rt, check.names = F)
dat=gather(rt1, "cellType","proportion", -ID, -Risk)

# 绘制ComplexHeatmap图
library(ComplexHeatmap)
library(RColorBrewer)
library(scales)
library(dplyr)

dat_complex <- dat %>% arrange(Risk) %>%
  pivot_wider(id_cols = ID, names_from = cellType, values_from = proportion)

## 加载临床信息
sample_clinical=read.table(cliFile, header=T, sep="\t", check.names=F)
rownames(sample_clinical)=sample_clinical[,1]
sameSample=intersect(row.names(risk), row.names(sample_clinical))
sample_clinical=cbind(sample_clinical[sameSample,,drop=F], risk[sameSample,"Risk",drop=F])
sample_clinical=sample_clinical[order(sample_clinical$Risk, decreasing=T),]

color_anno <- list()
color_anno[["Age"]] <- hue_pal()(length(levels(factor(sample_clinical$Age))))[1:length(levels(factor(sample_clinical$Age)))]
names(color_anno[["Age"]]) <- levels(factor(sample_clinical$Age))

color_anno[["Gender"]] <- hue_pal()(length(levels(factor(sample_clinical$Gender))))[1:length(levels(factor(sample_clinical$Gender)))]
names(color_anno[["Gender"]]) <- levels(factor(sample_clinical$Gender))

color_anno[["Race"]] <- hue_pal()(length(levels(factor(sample_clinical$Race))))[1:length(levels(factor(sample_clinical$Race)))]
names(color_anno[["Race"]]) <- levels(factor(sample_clinical$Race))

color_anno[["Stage"]] <- hue_pal()(length(levels(factor(sample_clinical$Stage))))[1:length(levels(factor(sample_clinical$Stage)))]
names(color_anno[["Stage"]]) <- levels(factor(sample_clinical$Stage))

color_anno[["Anatomic neoplasm subdivision"]] <- hue_pal()(length(levels(factor(sample_clinical$`Anatomic neoplasm subdivision`))))[1:length(levels(factor(sample_clinical$`Anatomic neoplasm subdivision`)))]
names(color_anno[["Anatomic neoplasm subdivision"]]) <- levels(factor(sample_clinical$`Anatomic neoplasm subdivision`))

color_anno[["Smoking"]] <- hue_pal()(length(levels(factor(sample_clinical$Smoking))))[1:length(levels(factor(sample_clinical$Smoking)))]
names(color_anno[["Smoking"]]) <- levels(factor(sample_clinical$Smoking))

color_anno[["Neoadjuvant therapy"]] <- hue_pal()(length(levels(factor(sample_clinical$`Neoadjuvant therapy`))))[1:length(levels(factor(sample_clinical$`Neoadjuvant therapy`)))]
names(color_anno[["Neoadjuvant therapy"]]) <- levels(factor(sample_clinical$`Neoadjuvant therapy`))

color_anno[["Group"]] <- hue_pal()(length(levels(factor(sample_clinical$Risk))))[1:length(levels(factor(sample_clinical$Risk)))]
names(color_anno[["Group"]]) <- levels(factor(sample_clinical$Risk))


ht_annotation <- 
  HeatmapAnnotation(Age = anno_simple(x = sample_clinical$Age, col = color_anno[["Age"]], border = FALSE)) %v%
  HeatmapAnnotation(Gender = anno_simple(x = sample_clinical$Gender, col = color_anno[["Gender"]], border = FALSE)) %v% 
  HeatmapAnnotation(Race = anno_simple(x = sample_clinical$Race, col = color_anno[["Race"]], border = FALSE)) %v%
  HeatmapAnnotation(Stage = anno_simple(x = sample_clinical$Stage, col = color_anno[["Stage"]], border = FALSE)) %v% 
  HeatmapAnnotation('Anatomic neoplasm subdivision' = anno_simple(x = sample_clinical$`Anatomic neoplasm subdivision`, col = color_anno[["Anatomic neoplasm subdivision"]], border = FALSE)) %v%
  HeatmapAnnotation(Smoking = anno_simple(x = sample_clinical$Smoking, col = color_anno[["Smoking"]], border = FALSE)) %v%
  HeatmapAnnotation('Neoadjuvant therapy' = anno_simple(x = sample_clinical$`Neoadjuvant therapy`, col = color_anno[["Neoadjuvant therapy"]], border = FALSE)) %v%
  HeatmapAnnotation(Group = anno_simple(x = sample_clinical$Risk, col = color_anno[["Group"]], border = FALSE)) %v%
  HeatmapAnnotation('Cell Type' = anno_barplot(x = dat_complex %>% select(levels(factor(dat$cellType))), 
                                             gp = gpar(fill = hue_pal()(22)[1:22], col = hue_pal()(22)[1:22]),
                                             bar_width = 1, 
                                             border = FALSE,
                                             height = unit(5, "cm")),
                    annotation_name_rot = 90
  )
lgd_list = list(
  Legend(labels = levels(factor(sample_clinical$Age)), title = "Age", 
         legend_gp = gpar(fill = color_anno[["Age"]])),
  Legend(labels = levels(factor(sample_clinical$Gender)), title = "Gender", 
         legend_gp = gpar(fill = color_anno[["Gender"]])),
  Legend(labels = levels(factor(sample_clinical$Race)), title = "Race", 
         legend_gp = gpar(fill = color_anno[["Race"]])),
  Legend(labels = levels(factor(sample_clinical$Stage)), title = "Stage", 
         legend_gp = gpar(fill = color_anno[["Stage"]])),
  Legend(labels = levels(factor(sample_clinical$`Anatomic neoplasm subdivision`)), title = "Anatomic neoplasm subdivision", 
         legend_gp = gpar(fill = color_anno[["Anatomic neoplasm subdivision"]])),
  Legend(labels = levels(factor(sample_clinical$Smoking)), title = "Smoking", 
         legend_gp = gpar(fill = color_anno[["Smoking"]])),
  Legend(labels = levels(factor(sample_clinical$`Neoadjuvant therapy`)), title = "Neoadjuvant therapy", 
         legend_gp = gpar(fill = color_anno[["Neoadjuvant therapy"]])),
  Legend(labels = levels(factor(sample_clinical$Risk)), title = "Group", 
         legend_gp = gpar(fill = color_anno[["Group"]])),
  Legend(labels = levels(factor(dat$cellType)), title = "Cell Type", 
         legend_gp = gpar(fill = hue_pal()(22)[1:22]))
)

pdf("complexheatmap.pdf", width = 8, height = 12)
draw(ht_annotation, heatmap_legend_list = lgd_list, heatmap_legend_side = "right")
dev.off()

