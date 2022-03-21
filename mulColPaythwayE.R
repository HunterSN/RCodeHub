rm(list = ls())
options(warn= -1)
#warnings('off')
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyr)
library(enrichplot)
library(DOSE)
library(openxlsx)
library(UpSetR)
library(flexclust)
library(proxy)
library(ggplot2)
library(factoextra)



gene = read.csv('HCCHD.csv',header = T)
#gene2EID = gene

SYM2EID = function(num){
  lie = as.character(gene[,num])
  keytypes(org.Hs.eg.db)
  gene.df <- bitr(lie, fromType = "SYMBOL", #fromType是指你的数据ID类型是属于哪一类的
                  toType = c("ENTREZID"), #toType是指你要转换成哪种ID类型，可以写多种，也可以只写一种
                  OrgDb = "org.Hs.eg.db")#Orgdb是指对应的注释包是哪个
  #gene2EID[,num] = as.character(gene.df[,2])
  head(gene.df)
  return(gene.df)
}




#####KEGG多列#########
kegg = as.data.frame(c())
catch = SYM2EID(1)[,2]

ekegg = enrichKEGG(
  gene = catch,
  organism = "hsa",
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  use_internal_data = FALSE
)
ekegg = DOSE::setReadable(ekegg,OrgDb = 'org.Hs.eg.db',keyType = 'ENTREZID')## 将 Gene ID 转换为 symbol
catch = as.data.frame(ekegg)
#catch = catch[which(catch$qvalue < 0.2),]
#catch = catch[which(catch$Count > fivenum(catch$Count)[4]),]

catch = unite(catch, "ID_Description", ID, Description,remove = T)
kegg = cbind(catch$ID_Description,catch$Count)
kegg = as.data.frame(kegg)
names(kegg) = c('ID_Description','Count')
for (m in 2:ncol(gene)) {
  print(m)
  catch = SYM2EID(m)[,2]
  ekegg = enrichKEGG(
    gene = catch,
    organism = "hsa",
    keyType = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    qvalueCutoff = 0.2,
    use_internal_data = FALSE
  )
  ekegg = DOSE::setReadable(ekegg,OrgDb = 'org.Hs.eg.db',keyType = 'ENTREZID')## 将 Gene ID 转换为 symbol
  catch = as.data.frame(ekegg)
  #catch = catch[which(catch$Count > fivenum(catch$Count)[4]),]
  catch = unite(catch, "ID_Description", ID, Description,remove = T)
  catch2 = as.data.frame(cbind(catch$ID_Description,catch$Count))
  names(catch2) = c('ID_Description','Count')
  kegg = merge(kegg,catch2,by = 'ID_Description', all = T)
}

names(kegg) = c('ID_Description',names(gene))

#kegg[is.na(kegg)] = 0
library(openxlsx)
write.xlsx(kegg, file = "kegg.xlsx", colNames = TRUE)




#####GOMF多列#########

goMF = as.data.frame(c())

catch = SYM2EID(1)[,2]

egoMF = enrichGO(
  gene = catch,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "MF",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE,
  pool = FALSE
)
egoMF = DOSE::setReadable(egoMF,OrgDb = 'org.Hs.eg.db',keyType = 'ENTREZID')## 将 Gene ID 转换为 symbol
catch = as.data.frame(egoMF)
#catch = catch[which((catch$pvalue < 0.01) & (catch$p.adjust < 0.01)),]
catch = catch[which(catch$Count > fivenum(catch$Count)[3]),]

catch = unite(catch, "ID_Description", ID, Description,remove = T)


goMF = cbind(catch$ID_Description,catch$Count)
goMF = as.data.frame(goMF)
names(goMF) = c('ID_Description','Count')




for (m in 2:ncol(gene)) {
  print(m)
  catch = SYM2EID(m)[,2]
  egoMF = enrichGO(
    gene = catch,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "MF",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    minGSSize = 10,
    maxGSSize = 500,
    readable = FALSE,
    pool = FALSE
  )
  egoMF = DOSE::setReadable(egoMF,OrgDb = 'org.Hs.eg.db',keyType = 'ENTREZID')## 将 Gene ID 转换为 symbol
  catch = as.data.frame(egoMF)
  #catch = catch[which((catch$pvalue < 0.01) & (catch$p.adjust < 0.01)),]
  catch = catch[which(catch$Count > fivenum(catch$Count)[3]),]
  catch = unite(catch, "ID_Description", ID, Description,remove = T)
  catch2 = as.data.frame(cbind(catch$ID_Description,catch$Count))
  names(catch2) = c('ID_Description','Count')
  goMF = merge(goMF,catch2,by = 'ID_Description', all = T)
}

names(goMF) = c('ID_Description',names(gene))


library(openxlsx)
write.xlsx(goMF, file = "goMF.xlsx", colNames = TRUE)

#####GOBP多列#########

goBP = as.data.frame(c())

catch = SYM2EID(1)[,2]

egoBP = enrichGO(
  gene = catch,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE,
  pool = FALSE
)
egoBP = DOSE::setReadable(egoBP,OrgDb = 'org.Hs.eg.db',keyType = 'ENTREZID')## 将 Gene ID 转换为 symbol
catch = as.data.frame(egoBP)
#catch = catch[which((catch$pvalue < 0.01) & (catch$p.adjust < 0.01)),]
catch = catch[which(catch$Count > fivenum(catch$Count)[3]),]

catch = unite(catch, "ID_Description", ID, Description,remove = T)


goBP = cbind(catch$ID_Description,catch$Count)
goBP = as.data.frame(goBP)
names(goBP) = c('ID_Description','Count')




for (m in 2:ncol(gene)) {
  print(m)
  catch = SYM2EID(m)[,2]
  egoBP = enrichGO(
    gene = catch,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    minGSSize = 10,
    maxGSSize = 500,
    readable = FALSE,
    pool = FALSE
  )
  egoBP = DOSE::setReadable(egoBP,OrgDb = 'org.Hs.eg.db',keyType = 'ENTREZID')## 将 Gene ID 转换为 symbol
  catch = as.data.frame(egoBP)
  #catch = catch[which((catch$pvalue < 0.01) & (catch$p.adjust < 0.01)),]
  catch = catch[which(catch$Count > fivenum(catch$Count)[3]),]
  catch = unite(catch, "ID_Description", ID, Description,remove = T)
  catch2 = as.data.frame(cbind(catch$ID_Description,catch$Count))
  names(catch2) = c('ID_Description','Count')
  goBP = merge(goBP,catch2,by = 'ID_Description', all = T)
}

names(goBP) = c('ID_Description',names(gene))


library(openxlsx)
write.xlsx(goBP, file = "goBP.xlsx", colNames = TRUE)





#####GOCC多列#########

goCC = as.data.frame(c())

catch = SYM2EID(1)[,2]

egoCC = enrichGO(
  gene = catch,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "CC",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE,
  pool = FALSE
)
egoCC = DOSE::setReadable(egoCC,OrgDb = 'org.Hs.eg.db',keyType = 'ENTREZID')## 将 Gene ID 转换为 symbol
catch = as.data.frame(egoCC)
#catch = catch[which((catch$pvalue < 0.01) & (catch$p.adjust < 0.01)),]
catch = catch[which(catch$Count > fivenum(catch$Count)[3]),]

catch = unite(catch, "ID_Description", ID, Description,remove = T)


goCC = cbind(catch$ID_Description,catch$Count)
goCC = as.data.frame(goCC)
names(goCC) = c('ID_Description','Count')




for (m in 2:ncol(gene)) {
  print(m)
  catch = SYM2EID(m)[,2]
  egoCC = enrichGO(
    gene = catch,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "CC",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    minGSSize = 10,
    maxGSSize = 500,
    readable = FALSE,
    pool = FALSE
  )
  egoCC = DOSE::setReadable(egoCC,OrgDb = 'org.Hs.eg.db',keyType = 'ENTREZID')## 将 Gene ID 转换为 symbol
  catch = as.data.frame(egoCC)
  #catch = catch[which((catch$pvalue < 0.01) & (catch$p.adjust < 0.01)),]
  catch = catch[which(catch$Count > fivenum(catch$Count)[3]),]
  catch = unite(catch, "ID_Description", ID, Description,remove = T)
  catch2 = as.data.frame(cbind(catch$ID_Description,catch$Count))
  names(catch2) = c('ID_Description','Count')
  goCC = merge(goCC,catch2,by = 'ID_Description', all = T)
}

names(goCC) = c('ID_Description',names(gene))


library(openxlsx)
write.xlsx(goCC, file = "goCC.xlsx", colNames = TRUE)




#####GOALL多列#########

goALL = as.data.frame(c())

catch = SYM2EID(1)[,2]

egoALL = enrichGO(
  gene = catch,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE,
  pool = FALSE
)
egoALL = DOSE::setReadable(egoALL,OrgDb = 'org.Hs.eg.db',keyType = 'ENTREZID')## 将 Gene ID 转换为 symbol
catch = as.data.frame(egoALL)
#catch = catch[which((catch$pvalue < 0.01) & (catch$p.adjust < 0.01)),]
catch = catch[which(catch$Count > fivenum(catch$Count)[3]),]

catch = unite(catch, "ID_Description", ID, Description,remove = T)


goALL = cbind(catch$ID_Description,catch$Count)
goALL = as.data.frame(goALL)
names(goALL) = c('ID_Description','Count')




for (m in 2:ncol(gene)) {
  print(m)
  catch = SYM2EID(m)[,2]
  egoALL = enrichGO(
    gene = catch,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "ALL",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    minGSSize = 10,
    maxGSSize = 500,
    readable = FALSE,
    pool = FALSE
  )
  egoALL = DOSE::setReadable(egoALL,OrgDb = 'org.Hs.eg.db',keyType = 'ENTREZID')## 将 Gene ID 转换为 symbol
  catch = as.data.frame(egoALL)
  #catch = catch[which((catch$pvalue < 0.01) & (catch$p.adjust < 0.01)),]
  catch = catch[which(catch$Count > fivenum(catch$Count)[3]),]
  catch = unite(catch, "ID_Description", ID, Description,remove = T)
  catch2 = as.data.frame(cbind(catch$ID_Description,catch$Count))
  names(catch2) = c('ID_Description','Count')
  goALL = merge(goALL,catch2,by = 'ID_Description', all = T)
}

names(goALL) = c('ID_Description',names(gene))


library(openxlsx)
write.xlsx(goALL, file = "goALL.xlsx", colNames = TRUE)













######富集可视化#######

barplot(ekegg, showCategory = 10)
dotplot(ekegg, showCategory = 10)
#plotGOgraph(ekegg)
#goplot(ekegg)
emapplot(ekegg, showCategory = 30)
cnetplot(ekegg, showCategory = 5)

library(enrichplot)
library(DOSE)

ekegg = DOSE::setReadable(ekegg,OrgDb = 'org.Hs.eg.db',keyType = 'ENTREZID')## 将 Gene ID 转换为 symbol
ego = DOSE::setReadable(ego,OrgDb = 'org.Hs.eg.db',keyType = 'ENTREZID')## 将 Gene ID 转换为 symbol
ekegg



#Bar plot#####
#oragene = enrichDGN(SYM2EID(2)[,2])
barplot(ekegg,showCategory = 20)
## 该函数默认参数为：
## enrichDGN(gene, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe,
##   minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2,
##   readable = FALSE)

#Dot Plot#####
#gseagene <- gseNCG(SYM2EID(2)[,2], nPerm=10000)
p1 <- dotplot(ekegg, showCategory=30) + ggtitle("dotplot for ORA")
p1
#p2 <- dotplot(gseagene, showCategory=30) + ggtitle("dotplot for GSEA")
#plot_grid(p1, p2, ncol=2)

#Gene-Concept Network######
cnetplot(ekegg,categorySize="pvalue",colorEdge = T)
cnetplot(ekegg,categorySize="geneNum",colorEdge = T)
cnetplot(ekegg,circular=T,colorEdge=T)## 圆形布局，给线条上色
#oragnx <- setReadable(oragene, 'org.Hs.eg.db', 'ENTREZID')  ## 将 Gene ID 转换为 symbol
#hn = as.character(as.vector(gene$HN))
#cnetplot(hn,showCategory = 5, foldChange = NULL, layout = "kk")

#Enrichment Map#####
emapplot(ekegg,pie_scale=1.5,layout = "kk")
#Heatmap######
heatplot(ekegg)

#UpSet Plot#####

upsetplot(ekegg)
upsetplot(ego)

oragene = enrichDGN(SYM2EID(2)[,2])
upsetplot(oragene)

library(UpSetR)
upset(fromList(go), 
      nsets = 15, 
      sets = c("HN","baizhu","banxia","cheqianzi","chenpi","danshen","danggui","dangshen","ezhu","fuling","gancao","huangqi","niuxi","qianshi","qumai"), 
      mainbar.y.label = "intersection size", 
      sets.x.label = "Numbers of GO enrichment gene", 
      main.bar.color = "#2a83a2", sets.bar.color = "#3b7960",
      mb.ratio = c(0.5, 0.5),
      order.by = "freq", 
      decreasing = c(TRUE,FALSE),
      queries = list(list(query = intersects, params = list('niuxi','banxia'),active = T),
                     list(query = intersects, params = list("cheqianzi", "gancao"), active=T),
                     list(query = intersects, params = list("HN","baizhu","banxia","cheqianzi","chenpi","danshen","danggui","dangshen","ezhu","fuling","gancao","huangqi","niuxi","qianshi","qumai"), active=T)
      )
)
#intersection size title y轴标题大小
#intersection size tick labels=1.3 y轴刻度标签大小
#set size title 左侧柱状图标题大小 
#set size tick labels 左侧柱状图刻度标签大小 
#set names 左侧柱状图分类标签大小 
#numbers above bars 柱状图上面数字的大小
upset(fromList(kegg), 
      nsets = 15, 
      sets = c("HN","baizhu","banxia","cheqianzi","chenpi","danshen","danggui","dangshen","ezhu","fuling","gancao","huangqi","niuxi","qianshi","qumai"), 
      mainbar.y.label = "intersection size", 
      sets.x.label = "Numbers of GO enrichment gene", 
      main.bar.color = "#2a83a2", sets.bar.color = "#3b7960",
      mb.ratio = c(0.5, 0.5),
      order.by = "freq", 
      decreasing = c(TRUE,FALSE),
      queries = list(list(query = intersects, params = list('niuxi','HN'),active = T),
                     list(query = intersects, params = list("HN","baizhu","banxia","cheqianzi","chenpi","danshen","danggui","dangshen","ezhu","fuling","gancao","huangqi","niuxi","qianshi","qumai"), active=T)
      )
)

#####HCT#####
#A way
library(flexclust)
library(proxy)
library(ggplot2)
library(factoextra)

nutrient = go
rownames(nutrient) = nutrient[,1]
nutrient = nutrient[-1]
for(i in 1:ncol(nutrient)){ nutrient[ , i] <- as.integer(nutrient[ , i]) }
nutrient[is.na(nutrient)] = 0
d<-dist(nutrient, method = 'euclidean')

fit<-hclust(d,method='ward.D')

plot(fit,hang=-1,cex=0.8)#hang的作用在于让名字显示得更加规整



g = plot(fit, labels = NULL, hang = 0.1,
         axes = TRUE, frame.plot = FALSE, ann = TRUE,
         main = "Cluster Dendrogram",
         sub = NULL, xlab = NULL, ylab = "Height")

#B way
library(stats)
nutrient
nutrient = kegg
rownames(nutrient) = nutrient[,1]
nutrient = nutrient[-1]
for(i in 1:ncol(nutrient)){ nutrient[ , i] <- as.integer(nutrient[ , i]) }
nutrient[is.na(nutrient)] = 0

means <- sapply(nutrient,mean);SD <- sapply(nutrient,sd)
nutrientScale <- scale(nutrient,center=means,scale=SD)  #scale 标准化函数
Dist <- dist(nutrientScale,method="euclidean")
heatmap(as.matrix(Dist),labRow = F,labCol = F)

ClusteModel <- hclust(Dist,method = "ward.D")
result <- cutree(ClusteModel,k=3)
result = as.data.frame(result)
plot(ClusteModel)


#####heatmap#####
library(pheatmap)
nutrient = kegg
rownames(nutrient) = nutrient[,1]
nutrient = nutrient[-1]
for(i in 1:ncol(nutrient)){ nutrient[ , i] <- as.integer(nutrient[ , i]) }
nutrient[is.na(nutrient)] = 0

kegghm =  pheatmap(nutrient, #表达数据
                   cluster_rows = T,#行聚类
                   cluster_cols = T,#列聚类
                   annotation_legend=TRUE, # 显示样本分类
                   show_rownames = T,# 显示行名
                   show_colnames = T,# 显示列名
                   scale = "row", #对行标准化
                   color =colorRampPalette(c("#8854d0", "#ffffff","#fa8231"))(100), # 热图基准颜色
                   cutree_cols = 3, #分割
                   cutree_rows = 3, #分割
                   treeheight_row = 50,  #树高
                   treeheight_col = 30,  #树高
                   clustering_distance_rows = 'euclidean', # 计算聚类间距的算法，可选'correlation', 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski'
                   clustering_method = 'ward.D' # 聚类方法, 可选'ward', 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'
)

kegghm

nutrient = go
rownames(nutrient) = nutrient[,1]
nutrient = nutrient[-1]
for(i in 1:ncol(nutrient)){ nutrient[ , i] <- as.integer(nutrient[ , i]) }
nutrient[is.na(nutrient)] = 0

gohm =  pheatmap(nutrient, #表达数据
                 cluster_rows = T,#行聚类
                 cluster_cols = T,#列聚类
                 annotation_legend=TRUE, # 显示样本分类
                 show_rownames = T,# 显示行名
                 show_colnames = T,# 显示列名
                 scale = "row", #对行标准化
                 color =colorRampPalette(c("#8854d0", "#ffffff","#fa8231"))(100), # 热图基准颜色
                 #cutree_cols = 3, #分割
                 #cutree_rows = 3, #分割
                 treeheight_row = 50,  #树高
                 treeheight_col = 30,  #树高
                 clustering_distance_rows = 'euclidean', # 计算聚类间距的算法，可选'correlation', 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski'
                 clustering_method = 'ward.D' # 聚类方法, 可选'ward', 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'
)


gohm

#####和弦图######

library(circlize)

mat = kegg
rownames(mat) = mat[,1]
mat = mat[-1]
for(i in 1:ncol(mat)){ mat[ , i] <- as.numeric(mat[ , i]) }
mat[is.na(mat)] = 0

chordDiagram(mat, 
             grid.col = grid_color, 
             grid.border = NULL, 
             transparency = alpha_usr, 
             row.col = NULL, 
             column.col = NULL, 
             order = NULL, 
             directional = conf$extra$directional, 
             direction.type = conf$extra$direction_type, 
             diffHeight = convert_height(2, "mm"), 
             reduce = 1e-05, xmax = NULL, self.link = 2, 
             symmetric = FALSE, 
             keep.diagonal = FALSE, 
             preAllocateTracks = NULL, 
             annotationTrack = c("name", "grid", "axis"), 
             annotationTrackHeight = convert_height(c(conf$extra$dist_name,conf$extra$width_circle), "mm"), link.border = NA, 
             link.lwd = par("lwd"), 
             link.lty = par("lty"), 
             link.sort = FALSE, 
             link.decreasing = TRUE, 
             link.largest.ontop = FALSE, 
             link.visible = conf$extra$link_visible, 
             link.rank = NULL, 
             link.overlap = FALSE, 
             scale = conf$extra$scale, 
             group = NULL, 
             big.gap = 10, 
             small.gap = 1)

chordDiagram(mat)

library(explore)
explore(go)


