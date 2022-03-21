rm(list = ls())

library(flexclust)
library(proxy)
library(ggplot2)
library(factoextra)
nutrient = read.csv('gscy.csv',header = T)


rownames(nutrient) = nutrient[,1]
nutrient = nutrient[-1]

#可视化层次聚类的过程
#nutrient.scale<-scale(nutrient)
#d<-dist(nutrient.scale, method = 'euclidean')
d<-dist(nutrient, method = 'euclidean')

fit<-hclust(d,method='ward.D')

plot(fit,hang=-1,cex=0.8)#hang的作用在于让名字显示得更加规整



g = plot(fit, labels = NULL, hang = 0.1,
     axes = TRUE, frame.plot = FALSE, ann = TRUE,
     main = "Cluster Dendrogram",
     sub = NULL, xlab = "herb", ylab = "Height")




heatmap(as.matrix(d))

result = cutree(fit,k=5) #cutree函数提取每个样本所属的类别
result
temp = cmdscale(d, k=2)#cmdscale数据降维
x = temp[,1]

y = temp[,2]

gsffenlei = read.csv('GSFfenlei.csv',header = T)
rownames(gsffenlei) = gsffenlei[,1]
gsffenlei = gsffenlei[,-1]

library(ggplot2)
library(ggrepel)
library(ggsci)
p = ggplot(data.frame(x,y),aes(x,y))
Clusters = factor(result)
Classification = gsffenlei
m = p+geom_point(size=5,alpha=0.8,aes(colour = Clusters,shape = Classification)) +
  geom_text_repel(aes(label = rownames(nutrient)))+
  theme_classic(base_size = 10)+
  scale_color_jama()+
  ggtitle('Distribution of Clustering Results') +theme(plot.title = element_text(size = 13,hjust = 0.5))+
  ylab('')+
  xlab("")
m + ggtitle('')

n = fviz_dend(fit, k = 6, 
          cex = 0.8, 
          k_colors = c("#374E55FF", "#DF8F44FF", "#00A1D5FF", "#B24745FF","#79AF97FF","#6A6599FF"),
          color_labels_by_k = TRUE, 
          rect = TRUE,
          lwd = 0.5,
          ggtheme = theme_classic()
)
n+ ggtitle('')

library(grid)

grid.newpage()
pushViewport(viewport(layout = grid.layout(1,2)))
vplayout = function(x,y)viewport(layout.pos.row = x,layout.pos.col = y)
#print(a,vp = vplayout(1,1:2))
print(n+ ggtitle('A') +theme(plot.title = element_text(size = 20,hjust = 0.5)),vp = vplayout(1,1))
print(m+ ggtitle('B') +theme(plot.title = element_text(size = 20,hjust = 0.5)),vp = vplayout(1,2))






