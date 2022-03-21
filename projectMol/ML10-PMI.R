rm(list = ls())


defSet = function(a){
  setInf = read.csv('herbCoupInf10.csv', header = T)
  #setInf = iris
  setInf$BRwr = 0
  setInf$BgoSim = 0
  setInf$BOB = 0
  setInf$BJac = 0
  setInf$BPMI = 0
  setInf$Blift = 0
  setInf$Bsupport = 0
  
  cutN = 3
  as.numeric(quantile(setInf$DisRwr)[cutN])
  
  
  setInf[which(setInf$DisRwr >= as.numeric(quantile(setInf$DisRwr)[cutN])),]$BRwr = 1
  setInf[which(setInf$goSim >= as.numeric(quantile(setInf$goSim)[cutN])),]$BgoSim = 1
  setInf[which(setInf$OB >= as.numeric(quantile(setInf$OB)[cutN])),]$BOB = 1
  setInf[which(setInf$Jaccard >= as.numeric(quantile(setInf$Jaccard)[cutN])),]$BJac = 1
  setInf[which(setInf$PMI >= as.numeric(quantile(setInf$PMI)[cutN])),]$BPMI = 1
  
  setInf[which(setInf$lift >= as.numeric(quantile(setInf$lift)[cutN])),]$Blift = 1
  
  #setInf[which(setInf$support >= as.numeric(quantile(setInf$support)[cutN])),]$Bsupport = 1
  #setInf$type = setInf$BRwr
  
  #setInf$type = setInf$BgoSim
  
  #setInf$type = setInf$BOB
  
  #setInf$type = setInf$BJac
  
  #setInf$type = setInf$BPMI
  
  setInf$type = setInf$Blift
  
  #setInf$type = setInf$Bsupport
  
  return(setInf)
}




####KNN交叉验证#####

KS = 11
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
  pcl = fitted(kknn(type ~ OB + DisRwr + goSim + Jaccard, setInf[-omit,], setInf[omit,], k = KS ,distance = 2 ))
  
},setInf)
#整合测试结果
wh = unlist(pred)
#grps 顺序已经被打乱，重新从小到大排序
setInf$grps = grps
knnT = as.data.frame(table(wh))
setInfF = setInf[order(setInf$grps),]

setInfFKNN = setInfF
knn_roc = roc(setInfFKNN$type, as.numeric(wh),levels = c('0', '1'),direction = "<")
plot(knn_roc)
knn_roc$auc
coords(knn_roc, "best")
aupr=AUC(obs=setInfFKNN$type,pred=as.numeric(wh),curve = "ROC", simplif=TRUE, main = "ROC curve")
aupr=AUC(obs=setInfFKNN$type,pred=as.numeric(wh),curve = "PR", simplif=TRUE, main = "PR curve")



####SVM交叉验证#####
setInf = defSet(1)

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
  pcl = predict(svm(type ~ OB + DisRwr + goSim + Jaccard, data = setInf[-omit,],  type = 'C-classification', kernel = 'radial'),newdata = setInf[omit,])
  
},setInf)
#整合测试结果
wh = unlist(pred)
#grps 顺序已经被打乱，重新从小到大排序
setInf$grps = grps
knnT = as.data.frame(table(wh))
setInfF = setInf[order(setInf$grps),]

setInfFSVM = setInfF

#混淆矩阵
MLSVM= table(setInfFSVM$type, wh, dnn =c("真实值","预测值"))
MLSVM = as.data.frame(MLSVM)



####GBRT交叉验证#####


setInf = defSet(1)

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
  mol_GBRTF = gbm(type ~ OB + DisRwr + goSim + Jaccard,
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
MLGBRT= table(setInfFGBRT$type, wh, dnn =c("真实值","预测值"))
MLGBRT = as.data.frame(MLSVM)



####Bayes交叉验证#####

setInf = defSet(1)

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

setInf$type = as.factor(setInf$type)
setInf$type = as.factor(setInf$type)

pred = lapply(1:v,function(i, setInf){
  omit = which(grps == i)
  pcl = predict(NaiveBayes(type ~ OB + DisRwr + goSim + Jaccard,setInf[-omit,]),newdata = setInf[omit,])
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
MLB= table(setInfFB$type, wh, dnn =c("真实值","预测值"))
MLB = as.data.frame(MLSVM)






####LG交叉验证#####

setInf = defSet(1)

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

setInf$type = as.factor(setInf$type)
setInf$type = as.factor(setInf$type)

pred = lapply(1:v,function(i, setInf){
  omit = which(grps == i)
  pcl = predict.glm(glm(type ~ OB + DisRwr + goSim + Jaccard,setInf[-omit,],family=binomial(link="logit")),newdata = setInf[omit,])
},setInf)
#整合测试结果
wh = unlist(pred)
#grps 顺序已经被打乱，重新从小到大排序
setInf$grps = grps
knnT = as.data.frame(table(wh))
setInfF = setInf[order(setInf$grps),]

setInfFLG = setInfF


#混淆矩阵
MLLG= table(setInfFLG$type, wh, dnn =c("真实值","预测值"))
MLLG = as.data.frame(MLSVM)


LG_roc = roc(setInfFLG$type, as.numeric(wh),levels = c('0', '1'),direction = "<")


plot(LG_roc)
LG_roc$auc
coords(LG_roc, "best")
aupr=AUC(obs=setInfFLG$type,pred=as.numeric(wh),curve = "ROC", simplif=TRUE, main = "ROC curve")
aupr=AUC(obs=setInfFLG$type,pred=as.numeric(wh),curve = "PR", simplif=TRUE, main = "PR curve")


setInfFLG$wh = as.numeric(wh)








#####smooth####
knn_roc = roc(setInfFKNN$type, as.numeric(wh),levels = c('0', '1'),direction = "<")
SVM_roc = roc(setInfFSVM$type, as.numeric(wh),levels = c('0', '1'),direction = "<")
GBRT_roc = roc(setInfFGBRT$type, as.numeric(wh),levels = c('0', '1'),direction = "<")
B_roc = roc(setInfFB$type, as.numeric(wh),levels = c('0', '1'),direction = "<")
LG_roc = roc(setInfFLG$type, as.numeric(wh),levels = c('0', '1'),direction = "<")

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
     main=""
)
legend("bottomright", legend=c("KNN", "SVM", "GBDT",'BN','LR'), col=mypal, lwd=lwdS,)


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
         print.auc=TRUE,print.auc.pattern=paste("AUC of BN =", round(B_roc$auc[1],3)), print.auc.x=0.7,print.auc.y=0.2,
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
     main="ROC curves of 5 machine learning models"
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
         print.auc=TRUE,print.auc.pattern=paste("AUC of BN =", round(B_roc$auc[1],3)), print.auc.x=0.7,print.auc.y=0.2,
         smooth = F,col = mypal[4])  
plot.roc(LG_roc,
         add=T,  # 增加曲线
         lwd=lwdS,
         print.auc=TRUE, print.auc.pattern=paste("AUC of LR =", round(LG_roc$auc[1],3)),print.auc.x=0.7,print.auc.y=0.1,
         smooth = F,col = mypal[5])  







#write.csv(setInfFLG,"LRML.csv",quote = F,row.names = F,fileEncoding = 'GBK')


save.image(file = "ML-PMI.RData")
