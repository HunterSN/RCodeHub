rm(list = ls())

library(clusterProfiler)
library(org.Hs.eg.db)
gene = read.csv('genes.csv',header = F)

gene = as.character(gene[,1])

keytypes(org.Hs.eg.db)

gene.df <- bitr(gene, fromType = "ENTREZID", #fromType是指你的数据ID类型是属于哪一类的
                toType = c("ENSEMBL", "ENTREZID","SYMBOL"), #toType是指你要转换成哪种ID类型，可以写多种，也可以只写一种
                OrgDb = "org.Hs.eg.db")#Orgdb是指对应的注释包是哪个
head(gene.df)


write.csv(gene.df,"genesChanged.csv",quote = F,row.names = F,fileEncoding = 'GBK')
