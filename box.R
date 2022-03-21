rm(list = ls())
library(ggplot2)
library(tidyverse)
library(ggsci)
mol = read.csv('mol2.csv',header = T)

p = ggplot(mol,aes(x = class, y = num,fill = class)) +geom_boxplot(alpha=0.7) + 
  scale_y_continuous(name = "ScoreC")+
  scale_x_discrete(name = "") +
  ggtitle("Boxplot of ScoreC in Different Model") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, face =  "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11)) 
p

p2 = p + scale_fill_nejm()
p2 +ggtitle("")

