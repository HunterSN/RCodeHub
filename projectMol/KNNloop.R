rm(list = ls())



defSet = function(a){
  setInf = read.csv('HNherbCoupInf.csv', header = T,fileEncoding = 'GBK')
  #setInf = iris
  setInf$Blift = 0
  cutN = 3
  as.numeric(quantile(setInf$DisRwr)[cutN])
  setInf[which(setInf$lift >= as.numeric(quantile(setInf$lift)[cutN])),]$Blift = 1
  setInf$type = setInf$Blift
  return(setInf)
}




####KNN循环#####
setInf = defSet(1)
library(kknn)
library(pROC)
library(modEvA)


library(class)
knnKS = as.data.frame(c())
#for (KS in 1:nrow(setInf)) {
for (KS in 1:138) {
  v = 10 #交叉次数
  grps = cut(1:nrow(setInf), v, labels = FALSE)[sample(1:nrow(setInf))]
  pred = lapply(1:v,function(i, setInf){
    omit = which(grps == i)
    pcl = fitted(kknn(type ~ OB + DisRwr + HNgoSim + Jaccard + Euc, setInf[-omit,], setInf[omit,], k = KS ,distance = 2 ))
    
  },setInf)
  wh = unlist(pred)
  setInf$grps = grps
  knnT = as.data.frame(table(wh))
  setInfF = setInf[order(setInf$grps),]
  setInfFKNN = setInfF
  
  knn_roc = roc(setInfFKNN$type, as.numeric(wh),levels = c('0', '1'),direction = "<")
  knnKS[KS,1] = as.numeric(knn_roc$auc)
  t = table(setInfFKNN$type, as.numeric(wh))
  knnKS[KS,2] = sum(diag(t))/sum(t) #准确率
  print(KS)
}
knnKS$V3 = 1:nrow(knnKS)
plot(knnKS$V3,knnKS$V1)






#############################




####################
plot(knn_roc)
knn_roc$auc
coords(knn_roc, "best")
aupr=AUC(obs=setInfFKNN$type,pred=as.numeric(wh),curve = "ROC", simplif=TRUE, main = "ROC curve")
aupr=AUC(obs=setInfFKNN$type,pred=as.numeric(wh),curve = "PR", simplif=TRUE, main = "PR curve")




plot(SVM_roc)
SVM_roc$auc
coords(SVM_roc, "best")
aupr=AUC(obs=setInfFSVM$type,pred=as.numeric(wh),curve = "ROC", simplif=TRUE, main = "ROC curve")
aupr=AUC(obs=setInfFSVM$type,pred=as.numeric(wh),curve = "PR", simplif=TRUE, main = "PR curve")




plot(GBRT_roc)
GBRT_roc$auc
coords(GBRT_roc, "best")
aupr=AUC(obs=setInfFGBRT$type,pred=as.numeric(wh),curve = "ROC", simplif=TRUE, main = "ROC curve")
aupr=AUC(obs=setInfFGBRT$type,pred=as.numeric(wh),curve = "PR", simplif=TRUE, main = "PR curve")




plot(B_roc)
B_roc$auc
coords(B_roc, "best")
aupr=AUC(obs=setInfFB$type,pred=as.numeric(wh),curve = "ROC", simplif=TRUE, main = "ROC curve")
aupr=AUC(obs=setInfFB$type,pred=as.numeric(wh),curve = "PR", simplif=TRUE, main = "PR curve")




plot(LG_roc)
LG_roc$auc
coords(LG_roc, "best")
aupr=AUC(obs=setInfFLG$type,pred=as.numeric(wh),curve = "ROC", simplif=TRUE, main = "ROC curve")
aupr=AUC(obs=setInfFLG$type,pred=as.numeric(wh),curve = "PR", simplif=TRUE, main = "PR curve")












