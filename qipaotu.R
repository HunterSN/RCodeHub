library(ggplot2)
library(ggtech)
library(ggsci)


HNpathway = read.csv('gao_kegg.csv',header = T)

x= HNpathway$PValue
y = factor(HNpathway$Pathway, levels = HNpathway$Pathway)

pr = ggplot(HNpathway, aes(x,y)) + geom_point() + geom_point(aes(size=Count,color=-1*log10(PValue)))+
  scale_color_gradient(low="green",high ="red") + 
  labs(x="PValue",y="Pathway name",title="Pathway enrichment")
pr



pr = ggplot(HNpathway, aes(PValue,Pathway)) + geom_point() + geom_point(aes(size=Count,color=-1*log10(PValue)))+
  scale_color_gradient(low="green",high ="red")
pr = pr + labs(color=expression(-log[10](PValue)),size="Count",  
               x="PValue",y="Pathway name",title="Pathway enrichment") +
  theme_bw()

pr
