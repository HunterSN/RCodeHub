rm(list = ls())

library(proxy)


#disIndex = read.csv('disIndex.csv', header = F)
#prescr = read.csv('C&H.csv', header = T)
#herbCoup = read.csv('herbCoup.csv', header = T)


disIndex = read.csv('HNdisIndex.csv', header = F)
prescr = read.csv('HNC&H.csv', header = T)
herbCoup = read.csv('HNherbCoup.csv', header = T,fileEncoding = 'GBK')




for (lie in 2:ncol(disIndex)) {
  print(lie)
  for (hang in 2:nrow(disIndex)) {
    a = subset(prescr, prescr$num == disIndex[1, lie] & prescr$herbs == as.character(disIndex[hang, 1]))
    disIndex[hang,lie] = nrow(a)
  }
}



nutrient = disIndex
rownames(nutrient) = nutrient[,1]
nutrient = nutrient[-1]
names(nutrient) = nutrient[1,]
nutrient = nutrient[-1,]




#####欧式距离######
dEuc<-dist(nutrient, method = 'euclidean')
dEuc = as.matrix(dEuc)

for (i in 1:nrow(herbCoup)) {
  print(i)
  herbCoup$Euc[i] = as.character(subset(dEuc, rownames(dEuc) == herbCoup[i,1], select = as.character(herbCoup[i,2])))
}

#####马氏距离######
dMan<-dist(nutrient, method = 'manhattan')
dMan = as.matrix(dMan)

for (i in 1:nrow(herbCoup)) {
  print(i)
  herbCoup$Man[i] = as.character(subset(dMan, rownames(dMan) == herbCoup[i,1], select = as.character(herbCoup[i,2])))
}


#####堪培拉距离######
dCan<-dist(nutrient, method = 'canberra')
dCan = as.matrix(dCan)

for (i in 1:nrow(herbCoup)) {
  print(i)
  herbCoup$Can[i] = as.character(subset(dCan, rownames(dCan) == herbCoup[i,1], select = as.character(herbCoup[i,2])))
}

#####二元变量的距离######
dBin<-dist(nutrient, method = 'binary')
dBin = as.matrix(dBin)

for (i in 1:nrow(herbCoup)) {
  print(i)
  herbCoup$Bin[i] = as.character(subset(dBin, rownames(dBin) == herbCoup[i,1], select = as.character(herbCoup[i,2])))
}

write.csv(herbCoup,"HNherbCoupOut.csv",quote = F,row.names = F,fileEncoding = 'GBK')


save.image(file = "HNDistence.RData")
