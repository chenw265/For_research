#install.packages("igraph")


library(igraph)        #引用包
input="10_turquoise_network.txt"      #输入文件
setwd("E:\\My_projects\\LUSC\\05WGCNAinDiff\\network")     #设置工作目录

#读取输入文件
data=read.table(input, header=T, sep="\t", check.names=F)
g=graph.data.frame(data)

#计算每个节点邻接基因的数据
node.num = table(c(data[,1], data[,2]))
h = hist(node.num,breaks=10,plot=F)
node.num.new = as.numeric(cut(node.num, h$breaks))+2
names(node.num.new) = names(node.num)
node.num = node.num.new

#准备网络参数
V(g)$size = node.num[match(names(components(g)$membership),names(node.num))]
V(g)$shape = "circle"
V(g)$lable.cex = 0.5
V(g)$color = gsub("(.*?)\\_(.*?)\\_.*", "\\2", input)
vertex.frame.color = NA
E(g)$color = "grey"
E(g)$weight = data$weight
V(g)$color =  "#FDAF91FF"

#绘制网络图
pdf(file="network.pdf", width=7, height=7)
plot(g,layout=layout_nicely,
     vertex.label.cex=V(g)$lable.cex,
     edge.width = E(g)$weight,
     edge.arrow.size=0,
     vertex.label.color="black",
     vertex.frame.color=vertex.frame.color,
     edge.color=E(g)$color,
     vertex.label.font=2,
     vertex.size=V(g)$size,
     edge.curved=0)
dev.off()
