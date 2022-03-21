library(ggplot2)
library(ggsci)

gsfc = read.csv('gsftarget.csv',header = T, stringsAsFactors =  F)

myLabel = as.vector(gsfc$target)   ## 转成向量，否则图例的标签可能与实际顺序不一致
myLabel = paste(myLabel, "(", round(gsfc$pinshu/ sum(gsfc$pinshu) * 100, 2), "%)        ", sep = "")   ## 用 round() 对结果保留两位小数


p = ggplot(gsfc,aes(x = '',y = pinshu,fill = target)) + 
  geom_bar(width = 1,stat = 'identity') +
  coord_polar(theta = 'y')+
  labs(x = "", y = "", title = "") + 
  theme(axis.ticks = element_blank()) + 
  theme(legend.title = element_blank(), legend.position = "right") + 
  scale_fill_discrete(breaks = gsfc$target, labels = myLabel) ## 将原来的图例标签换成现在的myLabel
p
p + theme(axis.text.x = element_blank())+ 
  #geom_text(aes(y = pinshu/2 + c(0, cumsum(pinshu)[-length(pinshu)]), x = sum(pinshu)/20, label = myLabel), size = 5)## 在图中加上百分比：x 调节标签到圆心的距离, y 调节标签的左右位置
  theme(panel.grid=element_blank()) +    ## 去掉白色圆框和中间的坐标线
  theme(panel.border=element_blank())+   ## 去掉最外层正方形的框框  
  scale_fill_d3(palette = c("category20c"))

p +theme(axis.text.x = element_blank())+ 
  #geom_text(aes(y = pinshu/2 + c(0, cumsum(pinshu)[-length(pinshu)]), x = sum(pinshu)/20, label = myLabel), size = 5)## 在图中加上百分比：x 调节标签到圆心的距离, y 调节标签的左右位置
  theme(panel.grid=element_blank()) +    ## 去掉白色圆框和中间的坐标线
  theme(panel.border=element_blank())+
  scale_fill_npg()




p + scale_fill_npg()+theme_bw()+   #选择Nature调色盘 
  theme(axis.text.x = element_text(size=rel(1.0),angle = 45,hjust = 1,color ="black"), #字体倾斜
        panel.grid =element_blank())

p + scale_fill_d3(palette = c("category20c"))+theme_bw()+ #选择D3调色盘
  theme(axis.text.x = element_text(size=rel(1.0),angle = 45,hjust = 1,color ="black"),
        panel.grid =element_blank(),#删除网格线
        plot.title = element_text(hjust = 0.5))#使标题居中

