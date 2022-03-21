rm(list = ls())

#######Apriori#####
library(Matrix)
library(arules) #apriori
library(arulesViz)

a = read.csv('HNC&H.csv',header = T,stringsAsFactors = F)#读取数据
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
write.xlsx(aa, file = "HNshaitapsapriori.xlsx", colNames = TRUE)

######PMI######

a = read.csv('HNC&H.csv',header = T,stringsAsFactors = F)  #读取数据
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

write.csv(caoyaohuxinxi,"HNpmi.csv",quote = F,row.names = F,fileEncoding = 'GBK')


aa$PMI = caoyaohuxinxi$PMI


######药对信息整合#######

molInf = read.csv('HNmolInf.csv', header = T, stringsAsFactors = F)

HM = read.csv('HNHM.csv',header = T,stringsAsFactors = F)
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
  herbIndex[i,5] = sum(as.numeric(catch$HNgoSim))
  
}
names(herbIndex) = c('lhs', 'rhs', 'OB', 'DisRwr','HNgoSim')
aaOut = merge(aa, herbIndex, all.x = T)

######症状Jaccard######


HM = read.csv('HNHM.csv',header = T,stringsAsFactors = F)
molDis = read.csv('HNmolDis.csv', header = T, stringsAsFactors = F)

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

write.csv(aaOut,"HNherbCoupInf.csv",quote = F,row.names = F,fileEncoding = 'GBK')

######验证######

setInfHN = read.csv('HNherbCoupInf.csv', header = T,fileEncoding = 'GBK')

setInfHN$Blift = 0
cutN = 3
setInfHN[which(setInfHN$lift >= as.numeric(quantile(setInfHN$lift)[cutN])),]$Blift = 1




setInfHN$type = setInfHN$Blift

#######LR######

library(pROC)
library(modEvA)




library(class)

v = 10 #交叉次数

#分割原始数据（对原始数据均匀分成v重，以类别标签返回
#通过随机sample方法，如sample(5) --> 3 2 1 4 5,以sample产生的数据作为索引下标，提取cut类别标签数据。
#这样做到来随机抽取类别数据
grps = cut(1:nrow(setInfHN), v, labels = FALSE)[sample(1:nrow(setInfHN))]
print(grps)
# 对每份数据分别运行ML函数 对随机数据的1/v作为测试集，剩下的作为训练集
# 对于lapply函数，后面要加上传递给function的参数data, cl, k. lapply的返回值为list

setInfHN$type = as.factor(setInfHN$type)
setInfHN$type = as.factor(setInfHN$type)

predHN = lapply(1:v,function(i, setInfHN){
  omit = which(grps == i)
  pcl = predict.glm(glm(type ~ OB + DisRwr + HNgoSim + Jaccard + PMI,setInfHN[-omit,],family=binomial(link="logit")),newdata = setInfHN[omit,])
},setInfHN)
#整合测试结果
wh = unlist(predHN)
#grps 顺序已经被打乱，重新从小到大排序
setInfHN$grps = grps
knnT = as.data.frame(table(wh))
setInfHNF = setInfHN[order(setInfHN$grps),]

setInfHNFLG = setInfHNF


#混淆矩阵
MLLG= table(setInfHNFLG$type, wh, dnn =c("真实值","预测值"))
MLLG = as.data.frame(MLSVM)

HNLG_roc = roc(setInfHNFLG$type, as.numeric(wh),levels = c('0', '1'),direction = "<")

plot(LG_roc)
LG_roc$auc
coords(HNLG_roc, "best")


plot(smooth(HNLG_roc), col = mypal[6],
     print.auc=TRUE, #display pAUC value on the plot with following options:
     
     print.auc.pattern=paste("AUC of HNLR =", round(HNLG_roc$auc[1],3)), 
     print.auc.x=0.7,print.auc.y=0.5,
     lwd=lwdS,
     main="ROC of HNLR"
)



setInfHNFLG$wh = as.numeric(wh)

write.csv(setInfHNFLG,"HNLR.csv",quote = F,row.names = F,fileEncoding = 'GBK')


save.image(file = "HNTestML.RData")

######KNN######


KS = 4 #M#an E#uc
setInf = defSet(1)

library(kknn)
library(pROC)
library(modEvA)


library(class)

v = 10 #交叉次数

#分割原始数据（对原始数据均匀分成v重，以类别标签返回
#通过随机sample方法，如sample(5) --> 3 2 1 4 5,以sample产生的数据作为索引下标，提取cut类别标签数据。
#这样做到来随机抽取类别数据
grps = cut(1:nrow(setInf), v, labels = FALSE)[sample(1:nrow(setInf))]
print(grps)
# 对每份数据分别运行ML函数 对随机数据的1/v作为测试集，剩下的作为训练集
# 对于lapply函数，后面要加上传递给function的参数data, cl, k. lapply的返回值为list
pred = lapply(1:v,function(i, setInf){
  omit = which(grps == i)
  pcl = fitted(kknn(type ~ OB + DisRwr + HNgoSim + Jaccard + Man, setInf[-omit,], setInf[omit,], k = KS ,distance = 2 ))
  
},setInf)
#整合测试结果
KNNwh = unlist(pred)
#grps 顺序已经被打乱，重新从小到大排序
setInf$grps = grps
knnT = as.data.frame(table(KNNwh))
setInfF = setInf[order(setInf$grps),]

setInfFKNN = setInfF
knn_roc = roc(setInfFKNN$type, as.numeric(KNNwh),levels = c('0', '1'),direction = "<")
plot(knn_roc)
knn_roc$auc
coords(knn_roc, "best")
aupr=AUC(obs=setInfFKNN$type,pred=as.numeric(KNNwh),curve = "ROC", simplif=TRUE, main = "ROC curve")
aupr=AUC(obs=setInfFKNN$type,pred=as.numeric(KNNwh),curve = "PR", simplif=TRUE, main = "PR curve")


