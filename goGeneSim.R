rm(list = ls())

library(clusterProfiler)
library(org.Hs.eg.db)
gene = read.csv('genes.csv',header = F)

gene = as.character(gene[,1])

keytypes(org.Hs.eg.db)

gene.df <- bitr(gene, fromType = "SYMBOL", #fromType是指你的数据ID类型是属于哪一类的
                toType = c("ENSEMBL", "ENTREZID"), #toType是指你要转换成哪种ID类型，可以写多种，也可以只写一种
                OrgDb = "org.Hs.eg.db")#Orgdb是指对应的注释包是哪个
head(gene.df)




library(GOSemSim)
d <- godata('org.Hs.eg.db', ont="MF", computeIC=FALSE)


genesim = gene.df[,-2]
genesim = genesim[!duplicated(genesim$SYMBOL),]
num = unique(as.vector(genesim$ENTREZID))
genetolist = t(combn(num,2)) 
rn = nrow(genetolist)
genetolist = as.data.frame(genetolist)



for (l in (1:rn)) {
  
  catch = geneSim(genetolist[l,1],genetolist[l,2],semData=d, measure = "Wang", drop = "IEA", combine = "BMA")
  genetolist[l,3] = catch[1]
  print(l)
}

names(genetolist) = c("geneA","geneB","goSim")

write.csv(genetolist,'gosim.csv',quote = F,row.names = F,fileEncoding = 'GBK')

write.csv(gene.df,'genedf.csv',quote = F,row.names = F,fileEncoding = 'GBK')

