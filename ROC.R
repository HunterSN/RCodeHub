rm(list = ls())
library(plotROC)
library(tidyverse)
library(ggplot2)
library(ggsci)

mypal <- pal_npg("nrc", alpha = 0.7)(9)
mypal
library("scales")
show_col(mypal)

aa = read.csv('gsffrocduib.csv',header = T)


names(aa)

#####封装函数####
plotROC <- function(.data, predict_col, target, group, positive=1, all=TRUE){
  if(!(require(tidyverse) & require(plotROC))){
    stop("--> tidyverse and plotROC packages are required..")
  } 
  
  predict_col <- enquo(predict_col)
  target <- enquo(target)
  group  <- enquo(group)
  
  predictN <- quo_name(predict_col)
  groupN   <- quo_name(group)
  
  df <- .data %>% dplyr::select(!! predict_col, !! target, !! group) %>%
    mutate(targetN = ifelse(!! target == positive, 1, 0)) %>% as.data.frame()
  if (all){
    df2 <- df 
    df2[, groupN] <- "ALL"
    
    df <- rbind(df, df2)
  }
  p  <- df %>%  ggplot(aes_string(m = predictN, 
                                  d = "targetN",
                                  color = groupN)) + geom_roc(show.legend = TRUE, labels=FALSE)
  p <- p + ggpubr::theme_classic2()
  
  ng <- levels(factor(df[, groupN]))
  if(length(ng) == 3){
    auc <- calc_auc(p)$AUC
    names(auc) <- ng
    auc <- base::sort(auc, decreasing = TRUE)
    p <- p + annotate("text", x = .75, y = .25, 
                      label = paste(names(auc)[1], " AUC =", round(auc[1], 3), "\n",
                                    names(auc)[2], " AUC =", round(auc[2], 3), "\n",
                                    names(auc)[3], " AUC =", round(auc[3], 3), "\n"),
                      size = 6)
  }
  
  p + xlab("1 - Specificity") + ylab("Sensitivity") + 
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
}

# 参数1：提供数据框
# 参数2：提供预测数值列
# 参数3：提供二分类信息列（尽量为0-1，字符也可以）
# 参数4：提供一个组别
# 参数5：这里1表示成功，如果target是success和failure，可以知道positive="success"
# 注意，这里只有3条曲线绘制时才会给出AUC在图上，可以修改函数进行自定义

#plotROC(longtest, predict_col = M, target = D, group = name, positive = 2)


#####单条#####
basicplot <- ggplot(aa, aes(d = cut, m = -m2socre )) + geom_roc(n.cuts = 0, labelsize = 5, labelround = 2)
basicplot
styledplot <- basicplot + style_roc()
styledplot
direct_label(basicplot, labels = "Biomarker", nudge_y = -.1) + style_roc()

######多条#######

aa = read.csv('gsffrocduib.csv',header = T)
####m2####
aa$Dm2socre = -aa$Dm2socre
longtest <- melt_roc(aa, "cut", c("m2socre","Bm2socre","Cm2socre","Dm2socre"))
longtest$name[which(longtest$name == "Bm2socre")] = "ScoreC Comparison1"
longtest$name[which(longtest$name == "Cm2socre")] = "ScoreC Comparison2"
longtest$name[which(longtest$name == "Dm2socre")] = "ScoreC Comparison3"
longtest$name[which(longtest$name == "m2socre")] = "ScoreC"
head(longtest)
p2 = ggplot(longtest, aes(d = D, m = -M, color = name)) + 
  geom_roc(n.cuts = 0, labelsize = 5, labelround = 2) + 
  style_roc() +
  ggsci::scale_color_npg()
p2
auc<-calc_auc(p2)
head(auc)
p2 = p2+annotate("text",x = .70, y = .35, ## 注释text的位置
           label = paste("AUC of ScoreC =", round(calc_auc(p2)$AUC[1], 3))) +
  annotate("text",x = .70, y = .25, ## 注释text的位置)
           label=paste("AUC of ScoreC Comparison1 =", round(calc_auc(p2)$AUC[2], 3)))+
  annotate("text",x = .70, y = .15, ## 注释text的位置)
           label=paste("AUC of ScoreC Comparison2 =", round(calc_auc(p2)$AUC[3], 3)))+
  annotate("text",x = .70, y = .05, ## 注释text的位置)
           label=paste("AUC of ScoreC Comparison3 =", round(calc_auc(p2)$AUC[4], 3)))+
  ggtitle('M2') +theme(plot.title = element_text(size = 17,hjust = 0.5))+
  theme(legend.title=element_blank())

p2

####m1####
#aa$Dm1socre = -aa$Dm1socre
longtest <- melt_roc(aa, "cut", c("m1socre","Bm1socre","Cm1socre","Dm1socre"))
longtest$name[which(longtest$name == "Bm1socre")] = "ScoreC Comparison1"
longtest$name[which(longtest$name == "Cm1socre")] = "ScoreC Comparison2"
longtest$name[which(longtest$name == "Dm1socre")] = "ScoreC Comparison3"
longtest$name[which(longtest$name == "m1socre")] = "ScoreC"
head(longtest)
p1 = ggplot(longtest, aes(d = D, m = -M, color = name)) + 
  geom_roc(n.cuts = 0, labelsize = 5, labelround = 2) + 
  style_roc() +
  ggsci::scale_color_npg()
p1
auc<-calc_auc(p1)
head(auc)
p1 = p1+annotate("text",x = .70, y = .35, ## 注释text的位置
           label = paste("AUC of ScoreC =", round(calc_auc(p1)$AUC[1], 3))) +
  annotate("text",x = .70, y = .25, ## 注释text的位置)
           label=paste("AUC of ScoreC Comparison1 =", round(calc_auc(p1)$AUC[2], 3)))+
  annotate("text",x = .70, y = .15, ## 注释text的位置)
           label=paste("AUC of ScoreC Comparison2 =", round(calc_auc(p1)$AUC[3], 3)))+
  annotate("text",x = .70, y = .05, ## 注释text的位置)
           label=paste("AUC of ScoreC Comparison3 =", round(calc_auc(p1)$AUC[4], 3)))+
  ggtitle('M1') +theme(plot.title = element_text(size = 17,hjust = 0.5))+
  theme(legend.title=element_blank())
p1
####m7####
aa$Dm7socre = -aa$Dm7socre
longtest <- melt_roc(aa, "cut", c("m7socre","Bm7socre","Cm7socre","Dm7socre"))
longtest$name[which(longtest$name == "Bm7socre")] = "ScoreC Comparison1"
longtest$name[which(longtest$name == "Cm7socre")] = "ScoreC Comparison2"
longtest$name[which(longtest$name == "Dm7socre")] = "ScoreC Comparison3"
longtest$name[which(longtest$name == "m7socre")] = "ScoreC"
head(longtest)
p7 = ggplot(longtest, aes(d = D, m = -M, color = name)) + 
  geom_roc(n.cuts = 0, labelsize = 5, labelround = 2) + 
  style_roc() +
  ggsci::scale_color_npg()
p7
auc<-calc_auc(p7)
head(auc)
p7 =p7+annotate("text",x = .70, y = .35, ## 注释text的位置
           label = paste("AUC of ScoreC =", round(calc_auc(p7)$AUC[1], 3))) +
  annotate("text",x = .70, y = .25, ## 注释text的位置)
           label=paste("AUC of ScoreC Comparison1 =", round(calc_auc(p7)$AUC[2], 3)))+
  annotate("text",x = .70, y = .15, ## 注释text的位置)
           label=paste("AUC of ScoreC Comparison2 =", round(calc_auc(p7)$AUC[3], 3)))+
  annotate("text",x = .70, y = .05, ## 注释text的位置)
           label=paste("AUC of ScoreC Comparison3 =", round(calc_auc(p7)$AUC[4], 3)))+
  ggtitle('M7') +theme(plot.title = element_text(size = 17,hjust = 0.5))+
  theme(legend.title=element_blank())
p7
####m0####
aa$Dm0socre = -aa$Dm0socre
longtest <- melt_roc(aa, "cut", c("m0socre","Bm0socre","Cm0socre","Dm0socre"))
longtest$name[which(longtest$name == "Bm0socre")] = "ScoreC Comparison1"
longtest$name[which(longtest$name == "Cm0socre")] = "ScoreC Comparison2"
longtest$name[which(longtest$name == "Dm0socre")] = "ScoreC Comparison3"
longtest$name[which(longtest$name == "m0socre")] = "ScoreC"
head(longtest)
p0 = ggplot(longtest, aes(d = D, m = -M, color = name)) + 
  geom_roc(n.cuts = 0, labelsize = 5, labelround = 2) + 
  style_roc() +
  ggsci::scale_color_npg()
p0
auc<-calc_auc(p0)
head(auc)
p0 = p0 +annotate("text",x = .70, y = .35, ## 注释text的位置
           label = paste("AUC of ScoreC =", round(calc_auc(p0)$AUC[1], 3))) +
  annotate("text",x = .70, y = .25, ## 注释text的位置)
           label=paste("AUC of ScoreC Comparison1 =", round(calc_auc(p0)$AUC[2], 3)))+
  annotate("text",x = .70, y = .15, ## 注释text的位置)
           label=paste("AUC of ScoreC Comparison2 =", round(calc_auc(p0)$AUC[3], 3)))+
  annotate("text",x = .70, y = .05, ## 注释text的位置)
           label=paste("AUC of ScoreC Comparison3 =", round(calc_auc(p0)$AUC[4], 3)))+
  ggtitle('M0') +theme(plot.title = element_text(size = 17,hjust = 0.5))+
  theme(legend.title=element_blank())
p0

####m6####
aa$Dm6socre = -aa$Dm6socre
longtest <- melt_roc(aa, "cut", c("m6socre","Bm6socre","Cm6socre","Dm6socre"))
longtest$name[which(longtest$name == "Bm6socre")] = "ScoreC Comparison1"
longtest$name[which(longtest$name == "Cm6socre")] = "ScoreC Comparison2"
longtest$name[which(longtest$name == "Dm6socre")] = "ScoreC Comparison3"
longtest$name[which(longtest$name == "m6socre")] = "ScoreC"
head(longtest)
p6 = ggplot(longtest, aes(d = D, m = -M, color = name)) + 
  geom_roc(n.cuts = 0, labelsize = 5, labelround = 2) + 
  style_roc() +
  ggsci::scale_color_npg()
p6
auc<-calc_auc(p6)
head(auc)
p6 = p6+annotate("text",x = .70, y = .35, ## 注释text的位置
           label = paste("AUC of ScoreC =", round(calc_auc(p6)$AUC[1], 3))) +
  annotate("text",x = .70, y = .25, ## 注释text的位置)
           label=paste("AUC of ScoreC Comparison1 =", round(calc_auc(p6)$AUC[2], 3)))+
  annotate("text",x = .70, y = .15, ## 注释text的位置)
           label=paste("AUC of ScoreC Comparison2 =", round(calc_auc(p6)$AUC[3], 3)))+
  annotate("text",x = .70, y = .05, ## 注释text的位置)
           label=paste("AUC of ScoreC Comparison3 =", round(calc_auc(p6)$AUC[4], 3)))+
  ggtitle('M6') +theme(plot.title = element_text(size = 17,hjust = 0.5))+
  theme(legend.title=element_blank())
p6
####HN####
aa$DHNsocre = -aa$DHNsocre
longtest <- melt_roc(aa, "cut", c("HNsocre","BHNsocre","CHNsocre","DHNsocre"))
longtest$name[which(longtest$name == "BHNsocre")] = "ScoreC Comparison1"
longtest$name[which(longtest$name == "CHNsocre")] = "ScoreC Comparison2"
longtest$name[which(longtest$name == "DHNsocre")] = "ScoreC Comparison3"
longtest$name[which(longtest$name == "HNsocre")] = "ScoreC"
head(longtest)
phn = ggplot(longtest, aes(d = D, m = -M, color = name)) + 
  geom_roc(n.cuts = 0, labelsize = 5, labelround = 2) + 
  style_roc() +
  ggsci::scale_color_npg()
phn
auc<-calc_auc(phn)
head(auc)
phn = phn+annotate("text",x = .70, y = .35, ## 注释text的位置
           label = paste("AUC of ScoreC =", round(calc_auc(phn)$AUC[1], 3))) +
  annotate("text",x = .70, y = .25, ## 注释text的位置)
           label=paste("AUC of ScoreC Comparison1 =", round(calc_auc(phn)$AUC[2], 3)))+
  annotate("text",x = .70, y = .15, ## 注释text的位置)
           label=paste("AUC of ScoreC Comparison2 =", round(calc_auc(phn)$AUC[3], 3)))+
  annotate("text",x = .70, y = .05, ## 注释text的位置)
           label=paste("AUC of ScoreC Comparison3 =", round(calc_auc(phn)$AUC[4], 3)))+
  ggtitle('HN') +theme(plot.title = element_text(size = 17,hjust = 0.5))+
  theme(legend.title=element_blank())
phn

par(mfrow = c(2,3))
p2
p1
p7
p0
p6
phn
p2 + ggtitle('A') +theme(plot.title = element_text(size = 20,hjust = 0.5))
library(grid)

grid.newpage()
pushViewport(viewport(layout = grid.layout(2,3)))
vplayout = function(x,y)viewport(layout.pos.row = x,layout.pos.col = y)
#print(a,vp = vplayout(1,1:2))
print(p2+ ggtitle('M2') +theme(plot.title = element_text(size = 17,hjust = 0.5)),vp = vplayout(1,1))
print(p1+ ggtitle('M1') +theme(plot.title = element_text(size = 17,hjust = 0.5)),vp = vplayout(1,2))
print(p7+ ggtitle('M7') +theme(plot.title = element_text(size = 17,hjust = 0.5)),vp = vplayout(1,3))
print(p0+ ggtitle('M0') +theme(plot.title = element_text(size = 17,hjust = 0.5)),vp = vplayout(2,1))
print(p6+ ggtitle('M6') +theme(plot.title = element_text(size = 17,hjust = 0.5)),vp = vplayout(2,2))
print(phn+ ggtitle('HN') +theme(plot.title = element_text(size = 17,hjust = 0.5)),vp = vplayout(2,3))

