rm(list = ls())
library(proxy)
library(dplyr)
library(kknn)
library(pROC)
library(modEvA)
library(class)


defSet = function(a){
  setInf = read.csv('gsffo.csv', header = T)
 
  #setInf[which(setInf$support >= as.numeric(quantile(setInf$support)[cutN])),]$Bsupport = 1
  #setInf$rmrcut = setInf$BRwr
  
  #setInf$rmrcut = setInf$BgoSim
  
  #setInf$rmrcut = setInf$BOB
  
  #setInf$rmrcut = setInf$BJac
  
  #setInf$rmrcut = setInf$BPMI
  
  #setInf$rmrcut = setInf$Blift
  
  
  return(setInf)
}


inf = read.csv('gsffo.csv', header = T)
modelinf = read.csv('MOL_target.csv',header = T)
wsx = read.csv('wsx.csv',header = T)

#确定化合物有效分割点
fencenga = quantile(inf$comrmrlog)
print(fencenga)
inf$rmrcut[which(inf$comrmrlog >= as.numeric(fencenga[4]))] = 1
inf$rmrcut[which(is.na(inf$rmrcut))] = 0


#####距离#########



molIndex = inf[,c(1,2)] 
molInf = merge(molIndex, modelinf, all.y = T)


for (i in 1:ncol(wsx)) {
  catch = matrix(nrow = nrow(molInf))  %>% as.data.frame() 
  names(catch) = names(wsx)[i]
  catch[which(molInf$target %in% wsx[,i]),1] = 1
  molInf = cbind(molInf, catch)
}



which(names(molInf) %in% names(wsx))



catch2 = molInf$MOL_ID %>% as.data.frame()
names(catch2) = 'MOL_ID'
catch2 = cbind(catch2, dplyr::mutate_all(molInf[,which(names(molInf) %in% names(wsx))],as.numeric))
#catch2 = distinct(catch2)
catch2[is.na(catch2)] = 0



nutrient = inf$MOL_ID %>% as.data.frame()




for (n in 1:ncol(wsx)) {
  catch = matrix(nrow = nrow(inf))  %>% as.data.frame()
  names(catch) = names(wsx)[n]
  for (m in 1:nrow(inf)) {
    catchT = subset(catch2, catch2$MOL_ID == inf$MOL_ID[m])
    catch[m,1] = catchT[,which(names(catchT) == names(wsx)[n])] %>% sum()
  }
  nutrient = cbind(nutrient, catch)
}


rownames(nutrient) = nutrient[,1]
nutrient = nutrient[-1]


distence = inf$MOL_ID %>% as.data.frame()
names(distence) = 'MOL_ID'
rwrScreening = subset(inf$MOL_ID, inf$rmrcut == 1) %>% as.matrix()
#####欧式距离######
dEuc<-dist(nutrient, method = 'euclidean')
dEuc = as.matrix(dEuc)



for (i in 1:nrow(distence)) {
  catchB = subset(dEuc, rownames(dEuc) == distence[i,1]) %>% as.data.frame()
  distence$Euc[i] = catchB[,which(names(catchB) %in% rwrScreening)] %>% sum()
  distence$EucT[i] = catchB[1,] %>% sum()
}

#####马氏距离######
dMan<-dist(nutrient, method = 'manhattan')
dMan = as.matrix(dMan)

for (i in 1:nrow(distence)) {
  catchB = subset(dMan, rownames(dMan) == distence[i,1]) %>% as.data.frame()
  distence$Man[i] = catchB[,which(names(catchB) %in% rwrScreening)] %>% sum()
  distence$ManT[i] = catchB[1,] %>% sum()
}



#####堪培拉距离######
dCan<-dist(nutrient, method = 'canberra')
dCan = as.matrix(dCan)


for (i in 1:nrow(distence)) {
  catchB = subset(dCan, rownames(dCan) == distence[i,1]) %>% as.data.frame()
  distence$Can[i] = catchB[,which(names(catchB) %in% rwrScreening)] %>% sum()
  distence$CanT[i] = catchB[1,] %>% sum()
}


#####二元变量的距离######
dBin<-dist(nutrient, method = 'binary')
dBin = as.matrix(dBin)

for (i in 1:nrow(distence)) {
  catchB = subset(dBin, rownames(dBin) == distence[i,1]) %>% as.data.frame()
  distence$Bin[i] = catchB[,which(names(catchB) %in% rwrScreening)] %>% sum()
  distence$BinT[i] = catchB[1,] %>% sum()
}

distence[is.na(distence)] = 0

#write.csv(distence,"distenceOut.csv",quote = F,row.names = F,fileEncoding = 'GBK')





molFeather = merge(inf, distence, all.x = T)

rm(list=ls()[which(ls()!='outT' & ls()!='molFeather')])



####KNN交叉验证#####

KS = 19
setInf = molFeather

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
  pcl = fitted(kknn(rmrcut ~ ob + jacaverage + HNgoSim +  ManT, setInf[-omit,], setInf[omit,], k = KS ,distance = 2 ))
  
},setInf)
#整合测试结果
wh = unlist(pred)
#grps 顺序已经被打乱，重新从小到大排序
setInf$grps = grps
knnT = as.data.frame(table(wh))
setInfF = setInf[order(setInf$grps),]

setInfFKNN = setInfF
knn_roc = roc(setInfFKNN$rmrcut, as.numeric(wh),levels = c('0', '1'),direction = "<")
plot(knn_roc)
knn_roc$auc
coords(knn_roc, "best")
aupr=AUC(obs=setInfFKNN$rmrcut,pred=as.numeric(wh),curve = "ROC", simplif=TRUE, main = "ROC curve")
aupr=AUC(obs=setInfFKNN$rmrcut,pred=as.numeric(wh),curve = "PR", simplif=TRUE, main = "PR curve")



knn_roc$levels











####SVM交叉验证#####
setInf = molFeather

library(e1071)
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
  pcl = predict(svm(rmrcut ~ ob + jacaverage + HNgoSim +  ManT, data = setInf[-omit,],  type = 'C-classification', kernel = 'radial'),newdata = setInf[omit,])
  
},setInf)
#整合测试结果
wh = unlist(pred)
#grps 顺序已经被打乱，重新从小到大排序
setInf$grps = grps
knnT = as.data.frame(table(wh))
setInfF = setInf[order(setInf$grps),]

setInfFSVM = setInfF

#混淆矩阵
MLSVM= table(setInfFSVM$rmrcut, wh, dnn =c("真实值","预测值"))
MLSVM = as.data.frame(MLSVM)
SVM_roc = roc(setInfFSVM$rmrcut, as.numeric(wh),levels = c('0', '1'),direction = "<")
plot(SVM_roc)


####GBRT交叉验证#####


setInf = molFeather

library(gbm)

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
  mol_GBRTF = gbm(rmrcut ~ ob + jacaverage + HNgoSim +  ManT,
                  distribution = "bernoulli",
                  data = setInf[-omit,],
                  var.monotone = NULL,
                  n.trees = 100,
                  interaction.depth = 1,
                  n.minobsinnode = 10,
                  shrinkage = 0.001,
                  bag.fraction = 0.5,
                  train.fraction = 1.0,
                  cv.folds=0,
                  keep.data = TRUE,
                  verbose = "CV",
                  class.stratify.cv=NULL,
                  n.cores = NULL)
  
  
  pcl = predict(mol_GBRTF,setInf[omit,],gbm.perf(mol_GBRTF,method = 'OOB'))
  
},setInf)
#整合测试结果
wh = unlist(pred)
#grps 顺序已经被打乱，重新从小到大排序
setInf$grps = grps
knnT = as.data.frame(table(wh))
setInfF = setInf[order(setInf$grps),]

setInfFGBRT = setInfF

#混淆矩阵
MLGBRT= table(setInfFGBRT$rmrcut, wh, dnn =c("真实值","预测值"))
MLGBRT = as.data.frame(MLSVM)

GBRT_roc = roc(setInfFGBRT$rmrcut, as.numeric(wh),levels = c('0', '1'),direction = "<")
plot(GBRT_roc)


####Bayes交叉验证#####

setInf = molFeather

library(klaR)
library(e1071)
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

setInf$rmrcut = as.factor(setInf$rmrcut)
setInf$rmrcut = as.factor(setInf$rmrcut)

pred = lapply(1:v,function(i, setInf){
  omit = which(grps == i)
  pcl = predict(NaiveBayes(rmrcut ~ ob + jacaverage + HNgoSim +  ManT,setInf[-omit,]),newdata = setInf[omit,])
  pre_BayesF = as.data.frame(pcl)[,2]
},setInf)
#整合测试结果
wh = unlist(pred)
#grps 顺序已经被打乱，重新从小到大排序
setInf$grps = grps
knnT = as.data.frame(table(wh))
setInfF = setInf[order(setInf$grps),]

setInfFB = setInfF

#混淆矩阵
MLB= table(setInfFB$rmrcut, wh, dnn =c("真实值","预测值"))
MLB = as.data.frame(MLSVM)



B_roc = roc(setInfFB$rmrcut, as.numeric(wh),levels = c('0', '1'),direction = "<")

plot(B_roc)


####LG交叉验证#####

setInf = molFeather

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

setInf$rmrcut = as.factor(setInf$rmrcut)
setInf$rmrcut = as.factor(setInf$rmrcut)

pred = lapply(1:v,function(i, setInf){
  omit = which(grps == i)
  pcl = predict.glm(glm(rmrcut ~ ob + jacaverage + HNgoSim +  ManT,setInf[-omit,],family=binomial(link="logit")),newdata = setInf[omit,])
},setInf)
#整合测试结果
wh = unlist(pred)
#grps 顺序已经被打乱，重新从小到大排序
setInf$grps = grps
knnT = as.data.frame(table(wh))
setInfF = setInf[order(setInf$grps),]

setInfFLG = setInfF


#混淆矩阵
MLLG= table(setInfFLG$rmrcut, wh, dnn =c("真实值","预测值"))
MLLG = as.data.frame(MLSVM)


LG_roc = roc(setInfFLG$rmrcut, as.numeric(wh),levels = c('0', '1'),direction = "<")


plot(LG_roc)
LG_roc$auc
coords(LG_roc, "best")
aupr=AUC(obs=setInfFLG$rmrcut,pred=as.numeric(wh),curve = "ROC", simplif=TRUE, main = "ROC curve")
aupr=AUC(obs=setInfFLG$rmrcut,pred=as.numeric(wh),curve = "PR", simplif=TRUE, main = "PR curve")


setInfFLG$wh = as.numeric(wh)








#####smooth####
#knn_roc = roc(setInfFKNN$rmrcut, as.numeric(wh),levels = c('0', '1'),direction = "<")
#SVM_roc = roc(setInfFSVM$rmrcut, as.numeric(wh),levels = c('0', '1'),direction = "<")
#GBRT_roc = roc(setInfFGBRT$rmrcut, as.numeric(wh),levels = c('0', '1'),direction = "<")
B_roc = roc(setInfFB$rmrcut, as.numeric(wh),levels = c('0', '1'),direction = "<")
#LG_roc = roc(setInfFLG$rmrcut, as.numeric(wh),levels = c('0', '1'),direction = "<")

library(ggplot2)
library(ggsci)
library("scales")
mypal <- pal_npg("nrc", alpha = 0.7)(9)
mypal
show_col(mypal)




plot(knn_roc)
lines(SVM_roc)
lines(GBRT_roc)
lines(B_roc)
lines(LG_roc)

lwdS = 3

plot(smooth(knn_roc), col = mypal[1],
     print.auc=TRUE, #display pAUC value on the plot with following options:
     
     print.auc.pattern=paste("AUC of KNN =", round(knn_roc$auc[1],3)), 
     print.auc.x=0.7,print.auc.y=0.5,
     lwd=lwdS,
     main="ROC of models"
)
legend("bottomright", legend=c("KNN", "SVM", "GBDT",'Bayes','LR'), col=mypal, lwd=lwdS,)


plot.roc(smooth(SVM_roc),
         add=T,  # 增加曲线
         lwd=lwdS,
         print.auc=TRUE, print.auc.pattern=paste("AUC of SVM =", round(SVM_roc$auc[1],3)),print.auc.x=0.7,print.auc.y=0.4,
         smooth = T,col = mypal[2])  
plot.roc(smooth(GBRT_roc),
         add=T,  # 增加曲线
         lwd=lwdS,
         print.auc=TRUE, print.auc.pattern=paste("AUC of GBDT =", round(GBRT_roc$auc[1],3)),print.auc.x=0.7,print.auc.y=0.3,
         smooth = T,col = mypal[3])  
plot.roc(smooth(B_roc),
         add=T,  # 增加曲线
         lwd=lwdS,
         print.auc=TRUE,print.auc.pattern=paste("AUC of Bayes =", round(B_roc$auc[1],3)), print.auc.x=0.7,print.auc.y=0.2,
         smooth = T,col = mypal[4])  
plot.roc(smooth(LG_roc),
         add=T,  # 增加曲线
         lwd=lwdS,
         print.auc=TRUE, print.auc.pattern=paste("AUC of LR =", round(LG_roc$auc[1],3)),print.auc.x=0.7,print.auc.y=0.1,
         smooth = T,col = mypal[5])  


auc = c(as.numeric(knn_roc$auc[1]), SVM_roc$auc[1], GBRT_roc$auc[1], B_roc$auc[1], LG_roc$auc[1])

coords(LG_roc, "best")



#####unsmooth#####
plot(knn_roc, col = mypal[1],
     print.auc=TRUE, 
     
     print.auc.pattern=paste("AUC of KNN =", round(knn_roc$auc[1],3)), 
     print.auc.x=0.7,print.auc.y=0.5,
     lwd=lwdS,
     main="ROC of models"
)
legend("bottomright", legend=c("KNN", "SVM", "GBDT",'Bayes','LR'), col=mypal, lwd=lwdS,)


plot.roc(SVM_roc,
         add=T,  # 增加曲线
         lwd=lwdS,
         print.auc=TRUE, print.auc.pattern=paste("AUC of SVM =", round(SVM_roc$auc[1],3)),print.auc.x=0.7,print.auc.y=0.4,
         smooth = F,col = mypal[2])  
plot.roc(GBRT_roc,
         add=T,  # 增加曲线
         lwd=lwdS,
         print.auc=TRUE, print.auc.pattern=paste("AUC of GBDT =", round(GBRT_roc$auc[1],3)),print.auc.x=0.7,print.auc.y=0.3,
         smooth = F,col = mypal[3])  
plot.roc(B_roc,
         add=T,  # 增加曲线
         lwd=lwdS,
         print.auc=TRUE,print.auc.pattern=paste("AUC of Bayes =", round(B_roc$auc[1],3)), print.auc.x=0.7,print.auc.y=0.2,
         smooth = F,col = mypal[4])  
#plot.roc(LG_roc,
#         add=T,  # 增加曲线
#         lwd=lwdS,
#         print.auc=TRUE, print.auc.pattern=paste("AUC of LR =", round(LG_roc$auc[1],3)),print.auc.x=0.7,print.auc.y=0.1,
#        smooth = F,col = mypal[5])  







#write.csv(setInfFLG,"LRML.csv",quote = F,row.names = F,fileEncoding = 'GBK')


#save.image(file = "ML913.RData")




rm(list = ls())
library(proxy)
library(dplyr)
library(kknn)
library(pROC)
library(modEvA)
library(class)


load('ML913.RData')

library(plotROC)
library(tidyverse)
library(ggplot2)
library(ggsci)

mypal <- pal_npg("nrc", alpha = 0.7)(9)
mypal
library("scales")
show_col(mypal)






plot(knn_roc, col = mypal[1],
     print.auc=TRUE, 
     
     print.auc.pattern=paste("AUC of KNN =", round(knn_roc$auc[1],3)), 
     print.auc.x=0.7,print.auc.y=0.5,
     lwd=lwdS,
     main="ROC of models"
)
legend("bottomright", legend=c("KNN", "SVM", "GBDT",'BN'), col=mypal, lwd=lwdS,)


plot.roc(SVM_roc,
         add=T,  # 增加曲线
         lwd=lwdS,
         print.auc=TRUE, print.auc.pattern=paste("AUC of SVM =", round(SVM_roc$auc[1],3)),print.auc.x=0.7,print.auc.y=0.4,
         smooth = F,col = mypal[2])  
plot.roc(GBRT_roc,
         add=T,  # 增加曲线
         lwd=lwdS,
         print.auc=TRUE, print.auc.pattern=paste("AUC of GBDT =", round(GBRT_roc$auc[1],3)),print.auc.x=0.7,print.auc.y=0.3,
         smooth = F,col = mypal[3])  
plot.roc(B_roc,
         add=T,  # 增加曲线
         lwd=lwdS,
         print.auc=TRUE,print.auc.pattern=paste("AUC of Bayes =", round(B_roc$auc[1],3)), print.auc.x=0.7,print.auc.y=0.2,
         smooth = F,col = mypal[4])  