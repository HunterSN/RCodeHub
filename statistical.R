rm(list = ls())
library(readxl)

kegg = read.csv('kegg.csv',header = T)
nutrient = read.csv('keggHCT.csv',header = T)

#kegg = read.csv('kegg_dis_s2.csv',header = T)
#nutrient = read.csv('kegg_dis_s2HCT.csv',header = T)




bp = read_xlsx('goBP.xlsx')
cc = read_xlsx('goCC.xlsx')
mf = read_xlsx('goMF.xlsx')
goAll = read_xlsx('goALL.xlsx')

herbIndex = read.csv('KNNML_0.527.csv',header = T)

keggS = subset(kegg, kegg$HCCHD != 'NA')

bpS = subset(bp, bp$HCCHD != 'NA') 
ccS = subset(cc, cc$HCCHD != 'NA') 
mfS = subset(mf, mf$HCCHD != 'NA') 
goAllS = subset(goAll, goAll$HCCHD != 'NA') 

####kegg#####
herbIndex$kegg = 0
for (i in 1:nrow(herbIndex)) {
  print(i)
  herbA = as.data.frame(keggS[,which(names(keggS) == as.character(herbIndex$lhs[i]))])
  herbB = as.data.frame(keggS[,which(names(keggS) == as.character(herbIndex$rhs[i]))])
  herbIndex$kegg[i] = sum(as.numeric(herbA[which(herbA != 'NA',),])) + sum(as.numeric(herbB[which(herbB != 'NA',),]))
  
}


####bp#####
herbIndex$bp = 0
for (i in 1:nrow(herbIndex)) {
  print(i)
  herbA = as.data.frame(bpS[,which(names(bpS) == as.character(herbIndex$lhs[i]))])
  herbB = as.data.frame(bpS[,which(names(bpS) == as.character(herbIndex$rhs[i]))])
  herbIndex$bp[i] = sum(as.numeric(herbA[which(herbA != 'NA',),])) + sum(as.numeric(herbB[which(herbB != 'NA',),]))
  
}


####cc#####
herbIndex$cc = 0
for (i in 1:nrow(herbIndex)) {
  print(i)
  herbA = as.data.frame(ccS[,which(names(ccS) == as.character(herbIndex$lhs[i]))])
  herbB = as.data.frame(ccS[,which(names(ccS) == as.character(herbIndex$rhs[i]))])
  herbIndex$cc[i] = sum(as.numeric(herbA[which(herbA != 'NA',),])) + sum(as.numeric(herbB[which(herbB != 'NA',),]))
  
}
####mf#####
herbIndex$mf = 0
for (i in 1:nrow(herbIndex)) {
  print(i)
  herbA = as.data.frame(mfS[,which(names(mfS) == as.character(herbIndex$lhs[i]))])
  herbB = as.data.frame(mfS[,which(names(mfS) == as.character(herbIndex$rhs[i]))])
  herbIndex$mf[i] = sum(as.numeric(herbA[which(herbA != 'NA',),])) + sum(as.numeric(herbB[which(herbB != 'NA',),]))
  
}
####goAll#####
herbIndex$goAll = 0
for (i in 1:nrow(herbIndex)) {
  print(i)
  herbA = as.data.frame(goAllS[,which(names(goAllS) == as.character(herbIndex$lhs[i]))])
  herbB = as.data.frame(goAllS[,which(names(goAllS) == as.character(herbIndex$rhs[i]))])
  herbIndex$goAll[i] = sum(as.numeric(herbA[which(herbA != 'NA',),])) + sum(as.numeric(herbB[which(herbB != 'NA',),]))
  
}



####T检验####
tT = function(a){
  y1 = herbIndex[which(as.numeric(herbIndex$KNNwh) > 0.514,),which(names(herbIndex) == a)]
  y2 = herbIndex[which(as.numeric(herbIndex$KNNwh) < 0.514,),which(names(herbIndex) == a)]
  t.test(y1,y2)
  
}

tT('kegg')
tT('bp')
tT('cc')
tT('mf')
tT('goAll')





####HCT########
library(flexclust)
library(proxy)
library(ggplot2)
library(factoextra)

quling = function(souce){
  nutrient = souce
  rownames(nutrient) = nutrient[,1]
  nutrient = nutrient[-1]
  for(i in 1:ncol(nutrient)){ nutrient[ , i] <- as.integer(nutrient[ , i]) }
  nutrient[is.na(nutrient)] = 0
  nutrient = nutrient[which(rowSums(nutrient) > 0),]
  return(nutrient)
}


nutrient = quling(nutrient)
d<-dist(nutrient, method = 'euclidean')

fit<-hclust(d,method='ward.D')

k = 5

result = cutree(fit,k = k) #cutree函数提取每个样本所属的类别
result
result = as.data.frame(result)
result = cbind(rownames(result), result)
names(result) = c(names(keggS)[1], 'class')
keggSR = merge(keggS, result)



catch = herbIndex$num
catchANOVA = as.data.frame(c())
classANOVA = as.data.frame(c())
for (n in 1:k) {
  print(n)
  catchV = c()
  for (m in 1:nrow(herbIndex)) {
    herbA = as.data.frame(subset(keggSR, keggSR$class == n, select = as.character(herbIndex$lhs[m])))
    herbB = as.data.frame(subset(keggSR, keggSR$class == n, select = as.character(herbIndex$rhs[m])))
    catchV = c(catchV, sum(as.numeric(herbA[which(herbA != 'NA',),])) + sum(as.numeric(herbB[which(herbB != 'NA',),])))
    catchANOVA[1,1] = sum(as.numeric(herbA[which(herbA != 'NA',),])) + sum(as.numeric(herbB[which(herbB != 'NA',),]))
    catchANOVA[1,2] = n
    catchANOVA[1,3] = m
    classANOVA = rbind(classANOVA, catchANOVA)
  }
  catchV = as.data.frame(catchV)
  catch = cbind(catch, catchV)
}

names(catch) = c('num', 1:k)

herbIndex = merge(herbIndex, catch)
names(classANOVA) = c('couple','class','num')
typeC = as.data.frame(cbind(herbIndex$num, herbIndex$type))
names(typeC) = c('num', 'type')
classANOVA = merge(classANOVA, typeC)

#######方差######

head(classANOVA)
summary(classANOVA)

#散点图分布#
library(ggpubr)
ggline(classANOVA, x = 'class', y = 'couple',
       add = c("mean_se", "jitter"), 
       order = 1:k,
       ylab = "couple", xlab = "class")

#正态性检验#
group_data = split(classANOVA$couple, classANOVA$class)
unlist(lapply(group_data, function(x){
  shapiro.test(x)$p.value
}))

#QQ图#
library(car)
qqPlot(group_data[[1]])


#方差齐性#
library(car)
leveneTest(as.numeric(classANOVA$couple),as.factor(classANOVA$class), center = median)


#检测离群点#
outlierTest(lm(couple ~ class, data = classANOVA))

#单因素方差分析#
#方差齐
summary(aov(couple ~ class, data = classANOVA))

#方差不齐
oneway.test(couple ~ class, data = classANOVA, var.equal = F)
kruskal.test(couple ~ class, data = classANOVA)
kruskal.test(couple ~ type, data = classANOVA)

#双因素方差
#方差齐
data = classANOVA
summary(aov(couple ~ class * type, data = classANOVA))





#####草药HCT######

library(flexclust)
library(proxy)
library(ggplot2)
library(factoextra)

quling = function(souce){
  nutrient = souce
  rownames(nutrient) = nutrient[,1]
  nutrient = nutrient[-1]
  for(i in 1:ncol(nutrient)){ nutrient[ , i] <- as.integer(nutrient[ , i]) }
  nutrient[is.na(nutrient)] = 0
  nutrient = nutrient[which(rowSums(nutrient) > 0),]
  return(nutrient)
}

nutrient2T = t(nutrient)
nutrient2T = quling(nutrient2T)
d2<-dist(nutrient2T, method = 'euclidean')

fit2<-hclust(d2,method='ward.D')




k2 = 5

result2 = cutree(fit2,k = k2) #cutree函数提取每个样本所属的类别
result2
resultA = as.data.frame(result2)
resultA = as.data.frame(cbind(rownames(resultA), resultA))

names(resultA) = c('lhs', 'herbAclass')
resultB = resultA
names(resultB) = c('rhs','herbBclass')


herbIndex = merge(herbIndex, resultA)
herbIndex = merge(herbIndex, resultB)


i = 1
guiyi = as.data.frame(scale(nutrient),scale = T)
herbA = as.data.frame(guiyi[,which(names(guiyi) == as.character(herbIndex[i,1]))])
herbB = as.data.frame(guiyi[,which(names(guiyi) == as.character(herbIndex[i,2]))])
(length(which(herbA >= 0,)) + length(which(herbB <=0,)))/nrow(guiyi)
#herbIndex$keggNum[i] = length(which(herbA != 'NA',)) + length(which(herbB != 'NA',))


herbIndex$keggNum = 0
for (i in 1:nrow(herbIndex)) {
  print(i)
  herbA = as.data.frame(keggS[,which(names(keggS) == as.character(herbIndex[i,1]))])
  herbB = as.data.frame(keggS[,which(names(keggS) == as.character(herbIndex[i,2]))])
  herbIndex$keggNum[i] = length(which(herbA > sum(nutrient)/(nrow(nutrient)*ncol(nutrient)),)) + length(which(herbB > sum(nutrient)/(nrow(nutrient)*ncol(nutrient)),))
  
}




y1 = herbIndex[which(herbIndex$herbAclass != herbIndex$herbBclass & as.numeric(herbIndex$KNNwh) > 0.527 & as.numeric(herbIndex$keggNum >= mean(herbIndex$keggNum)),),which(names(herbIndex) == 'keggNum')]
y2 = herbIndex[which(herbIndex$herbAclass != herbIndex$herbBclass& as.numeric(herbIndex$KNNwh) > 0.527 & as.numeric(herbIndex$keggNum < mean(herbIndex$keggNum)),),which(names(herbIndex) == 'keggNum')]
t.test(y1,y2)

y1 = herbIndex[which(herbIndex$herbAclass == herbIndex$herbBclass & as.numeric(herbIndex$KNNwh) > 0.527,),which(names(herbIndex) == 'keggNum')]
y2 = herbIndex[which(herbIndex$herbAclass != herbIndex$herbBclass& as.numeric(herbIndex$KNNwh) > 0.527,),which(names(herbIndex) == 'keggNum')]
t.test(y1,y2)

library(pander)
library(dplyr)
library(Epi)

saveSame = length(which(herbIndex$herbAclass == herbIndex$herbBclass & as.numeric(herbIndex$KNNwh) >= 0.527,))
saveDiff = length(which(herbIndex$herbAclass != herbIndex$herbBclass& as.numeric(herbIndex$KNNwh) >= 0.527,))
deleSame = length(which(herbIndex$herbAclass == herbIndex$herbBclass & as.numeric(herbIndex$KNNwh) < 0.527,))
deleDiff = length(which(herbIndex$herbAclass != herbIndex$herbBclass& as.numeric(herbIndex$KNNwh) < 0.527,))

chisqTabel = matrix(c(saveSame, saveDiff, deleSame,deleDiff), byrow = T ,nrow = 2)
rownames(chisqTabel) = c('Save', 'Dele')
colnames(chisqTabel) = c('Same', 'Diff')
chisqTabel
pander(addmargins(chisqTabel))
twoby2(chisqTabel)
chisq.test(chisqTabel, correct = T)
#mcnemar.test(chisqTabel)





