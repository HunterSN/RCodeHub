rm(list = ls())

library(clusterProfiler)
library(org.Hs.eg.db)
gene = read.csv('genes.csv',header = F)
disgenes = read.csv('disgenes.csv',header = F)


######gene转换#####
gene = as.character(gene[,1])

keytypes(org.Hs.eg.db)

gene.df <- bitr(gene, fromType = "SYMBOL", #fromType是指你的数据ID类型是属于哪一类的
                toType = c("ENSEMBL", "ENTREZID"), #toType是指你要转换成哪种ID类型，可以写多种，也可以只写一种
                OrgDb = "org.Hs.eg.db")#Orgdb是指对应的注释包是哪个
head(gene.df)

######disGene转换######

geneDis = as.character(disgenes[,1])

keytypes(org.Hs.eg.db)

geneDis.df <- bitr(geneDis, fromType = "SYMBOL", #fromType是指你的数据ID类型是属于哪一类的
                toType = c("ENSEMBL", "ENTREZID"), #toType是指你要转换成哪种ID类型，可以写多种，也可以只写一种
                OrgDb = "org.Hs.eg.db")#Orgdb是指对应的注释包是哪个
head(geneDis.df)





library(GOSemSim)
d <- godata('org.Hs.eg.db', ont="MF", computeIC=FALSE)


genesim = gene.df[,-2]
genesim = genesim[!duplicated(genesim$SYMBOL),]


geneDissim = geneDis.df[,-2]
geneDissim = geneDissim[!duplicated(geneDissim$SYMBOL),]


genelist = as.data.frame(c())
catch = as.data.frame(geneDissim$ENTREZID)
for (n in 1:nrow(genesim)) {
  catch[,2] = genesim[n,2]
  genelist = rbind(genelist,catch)
  catch = as.data.frame(geneDissim$ENTREZID)
  print(n)
}


genetolist = as.data.frame(genelist)
rn = nrow(genetolist)
for (l in (1:rn)) {
  
  catch = geneSim(genetolist[l,1],genetolist[l,2],semData=d, measure = "Wang", drop = "IEA", combine = "BMA")
  genetolist[l,3] = catch[1]
  print(l)
}

names(genetolist) = c("geneDisE","geneE","goSim")
names(geneDissim) = c('geneDisS','geneDisE')
names(genesim) = c('geneS','geneE')

geneSimOut =merge(genetolist, geneDissim, all.x = T)
geneSimOut =merge(geneSimOut, genesim, all.x = T)
disgenesType = read.csv('disGenestype.csv',header = T)
names(disgenesType) = c('geneDisS','disType')
geneSimOut =merge(geneSimOut, disgenesType, all.x = T)


gene = read.csv('genes.csv',header = F)
for (m in 1:nrow(gene)) {
  catch = subset(geneSimOut, geneSimOut$geneS == gene[m,1] & geneSimOut$disType == "CHD" & geneSimOut$goSim != "NA")
  gene[m,2] = sum(as.numeric(catch$goSim))
  catch = subset(geneSimOut, geneSimOut$geneS == gene[m,1] & geneSimOut$disType == "HY" & geneSimOut$goSim != "NA")
  gene[m,3] = sum(as.numeric(catch$goSim))
  print(m)
}
names(gene) = c('gene','CHDSim','HYSim')



write.csv(gene,'geneSimOut.csv',quote = F,row.names = T,fileEncoding = 'GBK')
