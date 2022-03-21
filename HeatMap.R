rm(list = ls())
library(pheatmap)

#####通用函数#####
quling = function(souce){
  nutrient = souce
  rownames(nutrient) = nutrient[,1]
  nutrient = nutrient[-1]
  for(i in 1:ncol(nutrient)){ nutrient[ , i] <- as.integer(nutrient[ , i]) }
  nutrient[is.na(nutrient)] = 0
  nutrient = nutrient[which(rowSums(nutrient) > 0),]
  return(nutrient)
}


tongyongHM = function(tongyong,nhang,nlie,hangjulei,liejulie,hangming){
  nutrient = tongyong
  rownames(nutrient) = nutrient[,1]
  nutrient = nutrient[-1]
  for(i in 1:ncol(nutrient)){ nutrient[ , i] <- as.integer(nutrient[ , i]) }
  nutrient[is.na(nutrient)] = 0
  nutrient = nutrient[which(rowSums(nutrient) > 0),]
  
  tongyonghm =  pheatmap(nutrient, #表达数据
                     cluster_rows = hangjulei,#行聚类
                     cluster_cols = liejulie,#列聚类
                     annotation_legend=TRUE, # 显示样本分类
                     show_rownames = hangming,# 显示行名
                     show_colnames = T,# 显示列名
                     scale = "row", #对行标准化
                     #color =colorRampPalette(c("#8854d0", "#ffffff","#fa8231"))(100), # 热图基准颜色
                     color =colorRampPalette(c("navy", "white", "firebrick3"))(100),
                     cutree_cols = nlie, #分割
                     cutree_rows = nhang, #分割
                     treeheight_row = 50,  #树高
                     treeheight_col = 30,  #树高
                     border=F,#框线
                     clustering_distance_rows = 'euclidean', # 计算聚类间距的算法，可选'correlation', 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski'
                     clustering_method = 'ward.D' # 聚类方法, 可选'ward', 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'
  )
  return(tongyonghm)
}

tongyongHMNum = function(tongyong,nhang,nlie,hangjulei,liejulie,hangming){
  nutrient = tongyong
  rownames(nutrient) = nutrient[,1]
  nutrient = nutrient[-1]
  for(i in 1:ncol(nutrient)){ nutrient[ , i] <- as.integer(nutrient[ , i]) }
  nutrient[is.na(nutrient)] = 0
  nutrient = nutrient[which(rowSums(nutrient) > 0),]
  
  tongyonghm =  pheatmap(nutrient, #表达数据
                         cluster_rows = hangjulei,#行聚类
                         cluster_cols = liejulie,#列聚类
                         annotation_legend=TRUE, # 显示样本分类
                         show_rownames = hangming,# 显示行名
                         show_colnames = T,# 显示列名
                         scale = "row", #对行标准化
                         #color =colorRampPalette(c("#8854d0", "#ffffff","#fa8231"))(100), # 热图基准颜色
                         color =colorRampPalette(c("navy", "white", "firebrick3"))(100),
                         cutree_cols = nlie, #分割
                         cutree_rows = nhang, #分割
                         treeheight_row = 50,  #树高
                         treeheight_col = 30,  #树高
                         border=F,#框线
                         display_numbers = TRUE,number_color = "blue",
                         clustering_distance_rows = 'euclidean', # 计算聚类间距的算法，可选'correlation', 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski'
                         clustering_method = 'ward.D' # 聚类方法, 可选'ward', 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'
  )
  return(tongyonghm)
}


#####单个#####
library(pheatmap)
kegg = read.csv('kegg_dis_s2.csv',header = T)
q0 = quling(kegg)
tongyongHM(kegg,7,5,T,T,T)
tongyongHMNum(kegg,5,5,F,T,T)

go = read.csv('go.csv',header = T)
q0 = quling(go)
tongyongHM(go,4,5,T,T,T)

target = read.csv('gsft.csv',header = T)
q0 = quling(target)
tongyongHM(target,2,4,F,T,T)


mol = read.csv('gsfmol.csv',header = T)
q0 = quling(mol)
tongyongHM(mol,2,4,F,F,F)




