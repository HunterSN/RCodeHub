rm(list = ls())

par(mfrow = c(1,1))


target = read.csv('target.csv',header = T,stringsAsFactors = F)
symbolMol = read.csv('symbol.csv', header = T, stringsAsFactors = F)
symbolScreening = read.csv('symbolScreening.csv', header = T, stringsAsFactors = F)
ch1 = read.csv('ch1.csv',header = T)
hy = read.csv('hy.csv',header = T)
ch3 = read.csv('ch3.csv',header = T)
######条形图func####
ggplotBar = function(tableOfImputCol){
  library(ggplot2)
  sthf = as.data.frame(tableOfImputCol)
  tiao2 = ggplot(sthf,aes(x = sthf$Var1,y = sthf$Freq,fill = Var1))+geom_col( position="dodge") + 
    theme(panel.grid=element_blank()) +    ## 去掉白色圆框和中间的坐标线
    theme(panel.border=element_blank()) +
    theme(legend.title=element_blank()) +
    theme(legend.position = "none") + #删除图例
    theme(text = element_text(family = 'STSong')) +
    labs(title = "", x="", y="")  #修改坐标轴及标题
  tiao2 
}


####基因转换func#########

geneChanged = function(filename){
  library(clusterProfiler)
  library(org.Hs.eg.db)
  gene = filename
  
  gene = as.character(gene[,1])
  
  keytypes(org.Hs.eg.db)
  
  gene.df <- bitr(gene, fromType = "ENTREZID", #fromType是指你的数据ID类型是属于哪一类的
                  toType = c("ENSEMBL", "ENTREZID","SYMBOL"), #toType是指你要转换成哪种ID类型，可以写多种，也可以只写一种
                  OrgDb = "org.Hs.eg.db")#Orgdb是指对应的注释包是哪个
  head(gene.df)
  return(gene.df)
}

####表格合并#####

Mol_disSymbol = merge(symbolMol, symbolScreening, all = T)

CHD_Mol = Mol_disSymbol[which(Mol_disSymbol$chd == 1),]
HY_Mol = Mol_disSymbol[which(Mol_disSymbol$hy == 1),]
CHD_HY_Mol = merge(CHD_Mol, HY_Mol, all = T)


ggplotBar(table(CHD_HY_Mol$herb))
ggplotBar(table(CHD_Mol$herb))
ggplotBar(table(HY_Mol$herb))








#######ch1探针转换为GeneSymbol#######


ch1Note <-data.table::fread("ch1GSE42148_family.soft",skip ="ID")
ch1Symbol =merge(ch1, ch1Note, all.x = T)
#library(stringr)
#grep("Homo sapiens", ch1SymbolCol$Description)
#ch1SymbolCol = ch1Symbol[grep("Homo sapiens", ch1SymbolCol$Description),] 

ch1Symbol = ch1Symbol[!duplicated(ch1Symbol), ]
ch1Symbol = subset(hySymbol, ch1Symbol$GeneName!="NA")


#######ch1火山图#######

library(ggpubr)
library(ggthemes)

ch1.data = ch1Symbol

#对差异基因较正后p值（adj.p.val一列）进行log10转换
ch1.data$logP = -log10(ch1.data$adj.P.Val)


#新加一列Group
ch1.data$Group = 'not-significant'
#将adj.P.Val小于0.05，logFC大于2对基因设置为显著上调基因
#将adj.P.Val小于0.05，logFC小于2对基因设置为显著下调基因
logFCcut = 1
ch1.data$Group[which((ch1.data$adj.P.Val < 0.05) & (ch1.data$logFC > logFCcut))] = "up-regulated"
ch1.data$Group[which((ch1.data$adj.P.Val < 0.05) & (ch1.data$logFC < -logFCcut))] = "down-regulated"
#查看上调和下调基因数目
table(ch1.data$Group)

#新加一列Label
ch1.data$Label = ""
#对差异表达基因的P值进行从小到大排序
ch1.data = ch1.data[order(ch1.data$adj.P.Val),]
#高表达的基因中，选择adj.P.Val最小的10个
up.genes = head(ch1.data$GeneName[which(ch1.data$Group == 'up-regulated')],10)
#低表达的基因中，选择adj.P.Val最小的10个
down.genes = head(ch1.data$GeneName[which(ch1.data$Group == 'down-regulated')],10)
#将up.genes和down.genes合并，并加入到Label中
ch1deg.top10.genes = c(as.character(up.genes),as.character(down.genes))
ch1.data$Label[match(ch1deg.top10.genes,ch1.data$Label)] = ch1deg.top10.genes




#改变火山图点的颜色和坐标轴标注，使图片更美观
ggscatter(ch1.data, x = 'logFC', y = 'logP',
          color = "Group",
          palette = c('#2f5688','#BBBBBB','#CC0000'),
          size = 1,
          label = ch1.data$Label,
          font.label = 14,
          repel = T,
          xlab = 'log2FoldChange',
          ylab = '-log10(Adjust P-value)',) + theme_base() +
  geom_hline(yintercept = logFCcut, linetype = 'dashed') +
  geom_vline(xintercept = c(-logFCcut,logFCcut), linetype = 'dashed')

#########hy探针转换为GeneSymbol#######
load(file = 'GPL13825_probe_ID.Rdata')

hyNote <-data.table::fread("hyGSE76845_family.soft",skip ="ID")

names(hyNote)[5] = names(probe2ID)[1]

hyNoteSymbol =merge(hyNote, probe2ID,all.x = T)

hyNoteSymbolhe = hyNoteSymbol[,c(2, 14, 15)]

hySymbol =merge(hy, hyNoteSymbolhe, all.x = T)
hySymbol = hySymbol[!duplicated(hySymbol), ]
hySymbol = subset(hySymbol, hySymbol$gene_Symble!="NA")
###########hy火山图#######

  
library(ggpubr)
library(ggthemes)

hy.data = hySymbol

#对差异基因较正后p值（adj.p.val一列）进行log10转换
hy.data$logP = -log10(hy.data$adj.P.Val)

hy.data$Group = 'not-significant'

logFCcut = 3
hy.data$Group[which((hy.data$adj.P.Val < 0.05) & (hy.data$logFC > logFCcut))] = "up-regulated"
hy.data$Group[which((hy.data$adj.P.Val < 0.05) & (hy.data$logFC < -logFCcut))] = "down-regulated"

hyBrowsing = subset(hy.data, hy.data$trans_biotype == "protein_coding" )
hyBrowsingOut = subset(hyBrowsing, hyBrowsing$Group == "up-regulated" | hyBrowsing$Group == "down-regulated" )

hyBrowsingOut$lable = ''
#高表达的基因中，选择adj.P.Val最小的10个
up.genes = head(hyBrowsingOut$gene_Symble[which(hyBrowsingOut$Group == 'up-regulated')],10)
#低表达的基因中，选择adj.P.Val最小的10个
down.genes = head(hyBrowsingOut$gene_Symble[which(hyBrowsingOut$Group == 'down-regulated')],10)
#将up.genes和down.genes合并，并加入到Label中
deg.top10.genes = c(as.character(up.genes),as.character(down.genes))
hyBrowsingOut$lable[match(deg.top10.genes,hyBrowsingOut$gene_Symble)] = deg.top10.genes

table(hyBrowsingOut$Group)

#改变火山图点的颜色和坐标轴标注，使图片更美观
ggscatter(hyBrowsing, x = 'logFC', y = 'logP',
          color = "Group",
          palette = c('#2f5688','#BBBBBB','#CC0000'),
          size = 1,
          font.label = 14,
          repel = T,
          xlab = 'log2FoldChange',
          ylab = '-log10(Adjust P-value)', title = 'Genes related to hypertension') + theme_base() +
  geom_hline(yintercept = 1.30, linetype = 'dashed') +
  geom_vline(xintercept = c(-logFCcut,logFCcut), linetype = 'dashed')





########ch3火山图#######

library(ggpubr)
library(ggthemes)

ch3.data = ch3

#对差异基因较正后p值（adj.p.val一列）进行log10转换
ch3.data$logP = -log10(ch3.data$P.Val)


#新加一列Group
ch3.data$Group = 'not-significant'

logFCcut = 2
ch3.data$Group[which((ch3.data$P.Val < 0.05) & (ch3.data$logFC > logFCcut))] = "up-regulated"
ch3.data$Group[which((ch3.data$P.Val < 0.05) & (ch3.data$logFC < -logFCcut))] = "down-regulated"

ch3BrowsingOut = subset(ch3.data, ch3.data$Group == "up-regulated" | ch3.data$Group == "down-regulated")
library(stringr)
grep("non-protein", ch3BrowsingOut$Gene.title)
ch3BrowsingOut = ch3BrowsingOut[-grep("non-protein", ch3BrowsingOut$Gene.title),] 


table(ch3BrowsingOut$Group)
ch3BrowsingOut$lable = ''
#高表达的基因中，选择adj.P.Val最小的10个
up.genes = head(ch3BrowsingOut$Gene.symbol[which(ch3BrowsingOut$Group == 'up-regulated')],10)
#低表达的基因中，选择adj.P.Val最小的10个
down.genes = head(ch3BrowsingOut$Gene.symbol[which(ch3BrowsingOut$Group == 'down-regulated')],10)
#将up.genes和down.genes合并，并加入到Label中
deg.top10.genes = c(as.character(up.genes),as.character(down.genes))
ch3BrowsingOut$lable[match(deg.top10.genes,ch3BrowsingOut$Gene.symbol)] = deg.top10.genes



#改变火山图点的颜色和坐标轴标注，使图片更美观
ggscatter(ch3.data, x = 'logFC', y = 'logP',
          color = "Group",
          palette = c('#2f5688','#BBBBBB','#CC0000'),
          size = 1,
          font.label = 14,
          repel = T,
          xlab = 'log2FoldChange',
          ylab = '-log10(Adjust P-value)',title = 'Genes related to CHD') + theme_base() +
  geom_hline(yintercept = 1.3, linetype = 'dashed') +
  geom_vline(xintercept = c(-logFCcut,logFCcut), linetype = 'dashed')



########最后整理######

CHDGenes = as.data.frame(ch3BrowsingOut$Gene.symbol)
CHDGenes[,2] = 'CHD'
names(CHDGenes) = c('genes', 'dis')
CHDGenes = CHDGenes[!duplicated(CHDGenes), ]


HYGenes = as.data.frame(hyBrowsingOut$gene_Symble)
HYGenes[,2] = "HY"
names(HYGenes) = c('genes', 'dis')
HYGenes = HYGenes[!duplicated(HYGenes), ]


disGenes = merge(CHDGenes, HYGenes, all = T)
disGenes = disGenes[!duplicated(disGenes), ]
disGenes = disGenes[-which(disGenes$genes == ""), ]

library(openxlsx)
write.xlsx(disGenes, file = "disGenes.xlsx", colNames = TRUE)
write.xlsx(CHD_HY_Mol, file = "CHD_HY_Mol.xlsx", colNames = TRUE)
write.xlsx(hyBrowsingOut, file = "hygenes.xlsx", colNames = TRUE)
write.xlsx(ch3BrowsingOut, file = "CHDgenes.xlsx", colNames = TRUE)

