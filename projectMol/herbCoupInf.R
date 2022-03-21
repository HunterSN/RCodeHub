rm(list = ls())

#######Apriori#####
library(Matrix)
library(arules) #apriori
library(arulesViz)

a = read.csv('C&H.csv',header = T,stringsAsFactors = F)#读取数据
b_apriori = split(a$herb,a$num)#按处方分类组成list
c_apriori = apriori(as(b_apriori,'transactions'),parameter = list(supp = 0,conf = 0,maxlen = 2,minlen = 2)) #apriori计算
summary(c_apriori)
output =inspect(c_apriori)

shai = output
shai$lhs =  gsub('[{]',"",shai$lhs)
shai$lhs =  gsub('[}]',"",shai$lhs)

shai$rhs =  gsub('[{]',"",shai$rhs)
shai$rhs =  gsub('[}]',"",shai$rhs)

aa = subset(shai, (shai$lhs != 'other')&(shai$rhs != 'other'))
library(openxlsx)
write.xlsx(aa, file = "shaitapsapriori.xlsx", colNames = TRUE)

######PMI######

a = read.csv('C&H.csv',header = T,stringsAsFactors = F)  #读取数据
b = split(a$herb,a$num)  #按处方分类组成list
herb = unique(as.vector(a$herb)) #提取草药
ch = cbind(aa$lhs,aa$rhs) 
hangshu = nrow(ch)
caoyaohuxinxi = data.frame()
nrow(ch)
chufang = length(b)

#正式代码
for (i in 1:nrow(ch)) {
  print(i)
  countherb1 = 0
  countherb2 = 0
  tongshi = 0
  for (m in b) {
    if (is.element(ch[i,1],m)) {
      countherb1 = countherb1 + 1
    }
  }
  for (m in b) {
    if (is.element(ch[i,2],m)) {
      countherb2 = countherb2 + 1
    }
  }
  for (m in b) {
    if (is.element(ch[i,1],m)&&is.element(ch[i,2],m)) {
      tongshi = tongshi + 1
    }
  }
  huxinxi = abs(log((tongshi/chufang)/((countherb1/chufang)*(countherb2/chufang))))
  xx = ch[i,1]
  yy = ch[i,2]
  caoyaohuxinxi[i,1] = xx
  caoyaohuxinxi[i,2] = yy
  caoyaohuxinxi[i,3] = countherb1
  caoyaohuxinxi[i,4] = countherb2
  caoyaohuxinxi[i,5] = tongshi
  caoyaohuxinxi[i,6] = huxinxi
}

names(caoyaohuxinxi) = c("草药A","草药B","A频数","B频数","共用频数","PMI")

write.csv(caoyaohuxinxi,"m1pmi.csv",quote = F,row.names = F,fileEncoding = 'GBK')


aa$PMI = caoyaohuxinxi$PMI


######药对信息整合#######

molInf = read.csv('molInf.csv', header = T, stringsAsFactors = F)

HM = read.csv('HM.csv',header = T,stringsAsFactors = F)
HMI = merge(HM, molInf, all.x = T)
HMI = subset(HMI, HMI$OB!="NA")
HMI = subset(HMI, HMI$DisRwr!="#NUM!")
herbIndex = as.data.frame(cbind(aa$lhs, aa$rhs))

i = 1

for (i in 1:nrow(herbIndex)) {
  print(i)
  
  catch = subset(HMI, HMI$herb == herbIndex[i,1] | HMI$herb == herbIndex[i,2])
  
  herbIndex[i,3] = sum(catch$OB)
  herbIndex[i,4] = sum(as.numeric(catch$DisRwr))
  herbIndex[i,5] = sum(as.numeric(catch$HYgoSim), as.numeric(catch$CHDgoSim))
  
}
names(herbIndex) = c('lhs', 'rhs', 'OB', 'DisRwr','goSim')
aaOut = merge(aa, herbIndex, all.x = T)

######症状Jaccard######


HM = read.csv('HM.csv',header = T,stringsAsFactors = F)
molDis = read.csv('molDis.csv', header = T, stringsAsFactors = F)

herbIndex = as.data.frame(cbind(aaOut$lhs, aaOut$rhs))


for (i in 1:nrow(herbIndex)) {
  print(i)
  #catch = subset(HM, HM$herb == herbIndex[i,1] | HM$herb == herbIndex[i,2])
  catch = subset(HMI, HMI$herb == herbIndex[i,1] | HMI$herb == herbIndex[i,2])
  catchInter = catch[which(duplicated(catch$mol)),]
  catchUnion = catch[which(!duplicated(catch$mol)),]
  
  inter = sum(catchInter$dis)
  
  union = sum(catchUnion$dis)
  
  herbIndex[i,3] = inter/union
}

aaOut$Jaccard = herbIndex[,3]
names(herbIndex) = c('lhs', 'rhs', 'Jaccard')
aaOut = merge(aaOut, herbIndex, all.x = T)

write.csv(aaOut,"herbCoupInf.csv",quote = F,row.names = F,fileEncoding = 'GBK')
