library(ggpubr)
library(ggthemes)

deg.data = read.csv("vlp.csv",header = T)
head(deg.data)


#对差异基因较正后p值（adj.p.val一列）进行log10转换
deg.data$logP = -log10(deg.data$adj.P.Val)

#绘制基本热图
ggscatter(deg.data, x = 'logFC', y = 'logP') + theme_base()

#新加一列Group
deg.data$Group = 'not-significant'
#将adj.P.Val小于0.05，logFC大于2对基因设置为显著上调基因
#将adj.P.Val小于0.05，logFC小于2对基因设置为显著下调基因
logFCcut = 1.1
deg.data$Group[which((deg.data$adj.P.Val < 0.05) & (deg.data$logFC > logFCcut))] = "up-regulated"
deg.data$Group[which((deg.data$adj.P.Val < 0.05) & (deg.data$logFC < -logFCcut))] = "down-regulated"
#查看上调和下调基因数目
table(deg.data$Group)

#绘制新的火山图
ggscatter(deg.data,x = 'logFC', y = 'logP', color = 'Group', size = 1) + theme_base()


#改变火山图颜色（palette）和点的大小（size）
ggscatter(deg.data, x = 'logFC', y = 'logP',
          color = "Group",
          palette = c('green','gray','red'),
          size = 1) + theme_base()

#为火山图添加logP分界线（geom_hline）和logFC分界线（geom_vline)
ggscatter(deg.data, x = 'logFC', y = 'logP',
          color = "Group",
          palette = c('green','gray','red'),
          size = 1) + theme_base() +
  geom_hline(yintercept = 1.30, linetype = 'dashed') +
  geom_vline(xintercept = c(-logFCcut,logFCcut), linetype = 'dashed')

#新加一列Label
deg.data$Label = ""
#对差异表达基因的P值进行从小到大排序
deg.data = deg.data[order(deg.data$adj.P.Val),]
#高表达的基因中，选择adj.P.Val最小的10个
up.genes = head(deg.data$genesym[which(deg.data$Group == 'up-regulated')],10)
#低表达的基因中，选择adj.P.Val最小的10个
down.genes = head(deg.data$genesym[which(deg.data$Group == 'down-regulated')],10)
#将up.genes和down.genes合并，并加入到Label中
deg.top10.genes = c(as.character(up.genes),as.character(down.genes))
deg.data$Label[match(deg.top10.genes,deg.data$genesym)] = deg.top10.genes




#改变火山图点的颜色和坐标轴标注，使图片更美观
ggscatter(deg.data, x = 'logFC', y = 'logP',
          color = "Group",
          palette = c('#2f5688','#BBBBBB','#CC0000'),
          size = 1,
          label = deg.data$Label,
          font.label = 14,
          repel = T,
          xlab = 'log2FoldChange',
          ylab = '-log10(Adjust P-value)',) + theme_base() +
  geom_hline(yintercept = 1.3, linetype = 'dashed') +
  geom_vline(xintercept = c(-logFCcut,logFCcut), linetype = 'dashed')



write.csv(deg.data,"vlpout.csv",quote = F,row.names = F,fileEncoding = 'GBK')



