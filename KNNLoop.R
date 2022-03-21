rm(list=ls()[which(ls()!='outT' & ls()!='molFeather')])

set.seed(1234)


par(mfrow = c(1,1))

####KNN循环jacaver_ManT_19#####
setInf = molFeather
library(kknn)
library(pROC)
library(modEvA)
library(class)

knnKS = as.data.frame(c())
#for (KS in 1:nrow(setInf)) {
for (KS in 1:100) {
  v = 10 #交叉次数
  grps = cut(1:nrow(setInf), v, labels = FALSE)[sample(1:nrow(setInf))]
  pred = lapply(1:v,function(i, setInf){
    omit = which(grps == i)
    pcl = fitted(kknn(rmrcut ~ ob + jacaverage + HNgoSim +  ManT, setInf[-omit,], setInf[omit,], k = KS ,distance = 2 ))
    
  },setInf)
  wh = unlist(pred)
  setInf$grps = grps
  knnT = as.data.frame(table(wh))
  setInfF = setInf[order(setInf$grps),]
  setInfFKNN = setInfF
  
  knn_roc = roc(setInfFKNN$rmrcut, as.numeric(wh),levels = c('0', '1'),direction = "<")
  knnKS[KS,1] = as.numeric(knn_roc$auc)
  t = table(setInfFKNN$rmrcut, as.numeric(wh))
  knnKS[KS,2] = sum(diag(t))/sum(t) #准确率
  print(KS)
}
knnKS$V3 = 1:nrow(knnKS)
plot(knnKS$V3,knnKS$V1)


####KNN循环M2jacaver_ManT_19#####
setInf = molFeather
library(kknn)
library(pROC)
library(modEvA)
library(class)

knnKS = as.data.frame(c())
#for (KS in 1:nrow(setInf)) {
for (KS in 1:100) {
  v = 10 #交叉次数
  grps = cut(1:nrow(setInf), v, labels = FALSE)[sample(1:nrow(setInf))]
  pred = lapply(1:v,function(i, setInf){
    omit = which(grps == i)
    pcl = fitted(kknn(rmrcut ~ ob + jacaverage + m2goSim +  ManT, setInf[-omit,], setInf[omit,], k = KS ,distance = 2 ))
    
  },setInf)
  m2 = unlist(pred)
  setInf$grps = grps
  knnT = as.data.frame(table(m2))
  setInfF = setInf[order(setInf$grps),]
  setInfFKNN = setInfF
  
  knn_roc = roc(setInfFKNN$rmrcut, as.numeric(m2),levels = c('0', '1'),direction = "<")
  knnKS[KS,1] = as.numeric(knn_roc$auc)
  t = table(setInfFKNN$rmrcut, as.numeric(m2))
  knnKS[KS,2] = sum(diag(t))/sum(t) #准确率
  print(KS)
}
knnKS$V3 = 1:nrow(knnKS)
plot(knnKS$V3,knnKS$V1)


####KNN循环M1jacaver_ManT_12#####
setInf = molFeather
library(kknn)
library(pROC)
library(modEvA)
library(class)

knnKS = as.data.frame(c())
#for (KS in 1:nrow(setInf)) {
for (KS in 1:100) {
  v = 10 #交叉次数
  grps = cut(1:nrow(setInf), v, labels = FALSE)[sample(1:nrow(setInf))]
  pred = lapply(1:v,function(i, setInf){
    omit = which(grps == i)
    pcl = fitted(kknn(rmrcut ~ ob + jacaverage + m1goSim +  ManT, setInf[-omit,], setInf[omit,], k = KS ,distance = 2 ))
    
  },setInf)
  m1 = unlist(pred)
  setInf$grps = grps
  knnT = as.data.frame(table(m1))
  setInfF = setInf[order(setInf$grps),]
  setInfFKNN = setInfF
  
  knn_roc = roc(setInfFKNN$rmrcut, as.numeric(m1),levels = c('0', '1'),direction = "<")
  knnKS[KS,1] = as.numeric(knn_roc$auc)
  t = table(setInfFKNN$rmrcut, as.numeric(m1))
  knnKS[KS,2] = sum(diag(t))/sum(t) #准确率
  print(KS)
}
knnKS$V3 = 1:nrow(knnKS)
plot(knnKS$V3,knnKS$V1)


####KNN循环M7jacaver_ManT_22#####
setInf = molFeather
library(kknn)
library(pROC)
library(modEvA)
library(class)

knnKS = as.data.frame(c())
#for (KS in 1:nrow(setInf)) {
for (KS in 1:100) {
  v = 10 #交叉次数
  grps = cut(1:nrow(setInf), v, labels = FALSE)[sample(1:nrow(setInf))]
  pred = lapply(1:v,function(i, setInf){
    omit = which(grps == i)
    pcl = fitted(kknn(rmrcut ~ ob + jacaverage + m7goSim +  ManT, setInf[-omit,], setInf[omit,], k = KS ,distance = 2 ))
    
  },setInf)
  m7 = unlist(pred)
  setInf$grps = grps
  knnT = as.data.frame(table(m7))
  setInfF = setInf[order(setInf$grps),]
  setInfFKNN = setInfF
  
  knn_roc = roc(setInfFKNN$rmrcut, as.numeric(m7),levels = c('0', '1'),direction = "<")
  knnKS[KS,1] = as.numeric(knn_roc$auc)
  t = table(setInfFKNN$rmrcut, as.numeric(m7))
  knnKS[KS,2] = sum(diag(t))/sum(t) #准确率
  print(KS)
}
knnKS$V3 = 1:nrow(knnKS)
plot(knnKS$V3,knnKS$V1)

####KNN循环M0jacaver_ManT_19#####
setInf = molFeather
library(kknn)
library(pROC)
library(modEvA)
library(class)

knnKS = as.data.frame(c())
#for (KS in 1:nrow(setInf)) {
for (KS in 1:100) {
  v = 10 #交叉次数
  grps = cut(1:nrow(setInf), v, labels = FALSE)[sample(1:nrow(setInf))]
  pred = lapply(1:v,function(i, setInf){
    omit = which(grps == i)
    pcl = fitted(kknn(rmrcut ~ ob + jacaverage + m0goSim +  ManT, setInf[-omit,], setInf[omit,], k = KS ,distance = 2 ))
    
  },setInf)
  m0 = unlist(pred)
  setInf$grps = grps
  knnT = as.data.frame(table(m0))
  setInfF = setInf[order(setInf$grps),]
  setInfFKNN = setInfF
  
  knn_roc = roc(setInfFKNN$rmrcut, as.numeric(m0),levels = c('0', '1'),direction = "<")
  knnKS[KS,1] = as.numeric(knn_roc$auc)
  t = table(setInfFKNN$rmrcut, as.numeric(m0))
  knnKS[KS,2] = sum(diag(t))/sum(t) #准确率
  print(KS)
}
knnKS$V3 = 1:nrow(knnKS)
plot(knnKS$V3,knnKS$V1)

####KNN循环jacaver_EucT_19#####
setInf = molFeather
library(kknn)
library(pROC)
library(modEvA)
library(class)

knnKS = as.data.frame(c())
#for (KS in 1:nrow(setInf)) {
for (KS in 1:100) {
  v = 10 #交叉次数
  grps = cut(1:nrow(setInf), v, labels = FALSE)[sample(1:nrow(setInf))]
  pred = lapply(1:v,function(i, setInf){
    omit = which(grps == i)
    pcl = fitted(kknn(rmrcut ~ ob + jacaverage + HNgoSim +  EucT, setInf[-omit,], setInf[omit,], k = KS ,distance = 2 ))
    
  },setInf)
  wh = unlist(pred)
  setInf$grps = grps
  knnT = as.data.frame(table(wh))
  setInfF = setInf[order(setInf$grps),]
  setInfFKNN = setInfF
  
  knn_roc = roc(setInfFKNN$rmrcut, as.numeric(wh),levels = c('0', '1'),direction = "<")
  knnKS[KS,1] = as.numeric(knn_roc$auc)
  t = table(setInfFKNN$rmrcut, as.numeric(wh))
  knnKS[KS,2] = sum(diag(t))/sum(t) #准确率
  print(KS)
}
knnKS$V3 = 1:nrow(knnKS)
plot(knnKS$V3,knnKS$V1)



####KNN循环jacsum_ManT_12#####
setInf = molFeather
library(kknn)
library(pROC)
library(modEvA)
library(class)

knnKS = as.data.frame(c())
#for (KS in 1:nrow(setInf)) {
for (KS in 1:100) {
  v = 10 #交叉次数
  grps = cut(1:nrow(setInf), v, labels = FALSE)[sample(1:nrow(setInf))]
  pred = lapply(1:v,function(i, setInf){
    omit = which(grps == i)
    pcl = fitted(kknn(rmrcut ~ ob + jacsum + HNgoSim +  ManT, setInf[-omit,], setInf[omit,], k = KS ,distance = 2 ))
    
  },setInf)
  wh = unlist(pred)
  setInf$grps = grps
  knnT = as.data.frame(table(wh))
  setInfF = setInf[order(setInf$grps),]
  setInfFKNN = setInfF
  
  knn_roc = roc(setInfFKNN$rmrcut, as.numeric(wh),levels = c('0', '1'),direction = "<")
  knnKS[KS,1] = as.numeric(knn_roc$auc)
  t = table(setInfFKNN$rmrcut, as.numeric(wh))
  knnKS[KS,2] = sum(diag(t))/sum(t) #准确率
  print(KS)
}
knnKS$V3 = 1:nrow(knnKS)
plot(knnKS$V3,knnKS$V1)

####KNN循环jacsum_EucT_11#####
setInf = molFeather
library(kknn)
library(pROC)
library(modEvA)
library(class)

knnKS = as.data.frame(c())
#for (KS in 1:nrow(setInf)) {
for (KS in 1:100) {
  v = 10 #交叉次数
  grps = cut(1:nrow(setInf), v, labels = FALSE)[sample(1:nrow(setInf))]
  pred = lapply(1:v,function(i, setInf){
    omit = which(grps == i)
    pcl = fitted(kknn(rmrcut ~ ob + jacsum + HNgoSim +  EucT, setInf[-omit,], setInf[omit,], k = KS ,distance = 2 ))
    
  },setInf)
  wh = unlist(pred)
  setInf$grps = grps
  knnT = as.data.frame(table(wh))
  setInfF = setInf[order(setInf$grps),]
  setInfFKNN = setInfF
  
  knn_roc = roc(setInfFKNN$rmrcut, as.numeric(wh),levels = c('0', '1'),direction = "<")
  knnKS[KS,1] = as.numeric(knn_roc$auc)
  t = table(setInfFKNN$rmrcut, as.numeric(wh))
  knnKS[KS,2] = sum(diag(t))/sum(t) #准确率
  print(KS)
}
knnKS$V3 = 1:nrow(knnKS)
plot(knnKS$V3,knnKS$V1)



rm(list=ls()[which(ls()!='outT' & ls()!='molFeather')])
