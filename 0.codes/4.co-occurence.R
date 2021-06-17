# 加载包
library(igraph)
library(psych)
library(readxl)

# 读取otu-sample矩阵，行为sample，列为otu
# otu = read.table("otu_table.txt",head=T,row.names=1)

OTUs4cooccurance <- read_excel("../2.outputfiles/OTUs4cooccurance.xlsx", skip=0)
otu <- t(OTUs4cooccurance[,-1])
colnames(otu) <- OTUs4cooccurance[,1]$`OTU ID`

# 计算OTU间两两相关系数矩阵
# 数据量小时可以用psych包corr.test求相关性矩阵，数据量大时，可应用WGCNA中corAndPvalue, 但p值需要借助其他函数矫正
occor = corr.test(otu,use="pairwise",method="spearman",adjust="fdr",alpha=.05)
occor.r = occor$r # 取相关性矩阵R值
occor.p = occor$p # 取相关性矩阵p值

# 确定物种间存在相互作用关系的阈值，将相关性R矩阵内不符合的数据转换为0
occor.r[occor.p>0.1|abs(occor.r)<0.5] = 0 


# 构建igraph对象
igraph = graph_from_adjacency_matrix(occor.r,mode="undirected",weighted=TRUE,diag=FALSE)
igraph
# NOTE:可以设置weighted=NULL,但是此时要注意此函数只能识别相互作用矩阵内正整数，所以应用前请确保矩阵正确。
# 可以按下面命令转换数据
# occor.r[occor.r!=0] = 1
# igraph = graph_from_adjacency_matrix(occor.r,mode="undirected",weighted=NULL,diag=FALSE)

# 是否去掉孤立顶点，根据自己实验而定
# remove isolated nodes，即去掉和所有otu均无相关性的otu 可省略，前期矩阵已处理过
bad.vs = V(igraph)[degree(igraph) == 0]
igraph = delete.vertices(igraph, bad.vs)
igraph

# 将igraph weight属性赋值到igraph.weight
igraph.weight = E(igraph)$weight

# 做图前去掉igraph的weight权重，因为做图时某些layout会受到其影响
E(igraph)$weight = NA


# 简单出图
# 设定随机种子数，后续出图都从同一随机种子数出发，保证前后出图形状相对应
set.seed(123)
plot(igraph,main="Co-occurrence network",vertex.frame.color=NA,vertex.label=NA,edge.width=1,
     vertex.size=5,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))

# 如果构建网络时，weighted=NULL,此步骤不能统计
sum(igraph.weight>0)# number of postive correlation
sum(igraph.weight<0)# number of negative correlation

# set edge color，postive correlation 设定为red, negative correlation设定为blue
E.color = igraph.weight
E.color = ifelse(E.color>0, "red",ifelse(E.color<0, "blue","grey"))
E(igraph)$color = as.character(E.color)

# 改变edge颜色后出图
set.seed(123)
plot(igraph,main="Co-occurrence network",vertex.frame.color=NA,vertex.label=NA,edge.width=1,
     vertex.size=5,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))

# 可以设定edge的宽 度set edge width，例如将相关系数与edge width关联
E(igraph)$width = abs(igraph.weight)*4

# 改变edge宽度后出图
set.seed(123)
plot(igraph,main="Co-occurrence network",vertex.frame.color=NA,vertex.label=NA,
     vertex.size=5,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))


# 添加OTU注释信息，如分类单元和丰度
# 另外可以设置vertices size, vertices color来表征更多维度的数据
# 注意otu_pro.txt文件为我随机产生的数据，因此网络图可能不会产生特定的模式或规律。
# otu_pro = read.table("D:/Jiangnan Univeristy/ZhangzhiJunfen/20210612/3.new_outputs/otu_pro.txt",head=T,row.names=1)
otu_pro <- read.delim("../1.data/otu_pro.txt", row.names=1)
# set vertices size
igraph.size = otu_pro[V(igraph)$name,] # 筛选对应OTU属性
# colnames(igraph.size) <- c("test", "abundance", "phylum")
igraph.size1 = log((igraph.size$sum)*100) # 原始数据是什么，为什么*100再取e对数
V(igraph)$size = igraph.size1

# set vertices color
igraph.col = otu_pro[V(igraph)$name,]
as.factor(igraph.col$phylum)
levels(igraph.col$phylum)
levels(igraph.col$phylum) = c("green","deeppink","deepskyblue","yellow","brown","pink","gray","cyan","peachpuff") # 直接修改levles可以连值全部对应替换
#levels(igraph.col$phylum) = c("green","yellow","brown","pink","gray") # 直接修改levles可以连值全部对应替换
V(igraph)$color = as.character(igraph.col$color)

set.seed(123)
plot(igraph,main="Co-occurrence network",vertex.frame.color=NA,vertex.label=NA,
     edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))




# 老婆大人亲笔
# set vertices type.convert(geiguanqijieluanlaidaoluanheiheihei)









set.seed(123)
plot(igraph,main="Co-occurrence network",layout=layout_with_kk,vertex.frame.color=NA,vertex.label=NA,
     edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))

set.seed(123)
plot(igraph,main="Co-occurrence network",layout=layout.fruchterman.reingold,vertex.frame.color=NA,vertex.label=NA,
     edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))

set.seed(123)
plot(igraph,main="Co-occurrence network",layout=layout.fruchterman.reingold, vertex.frame.color=NA,vertex.label=V(igraph)$name,
     edge.lty=1,edge.curved=FALSE,margin=c(0,0,0,0))