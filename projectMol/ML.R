rm(list = ls())
set.seed(123)
defSet = function(a){
    setInf = read.csv('herbCoupInf10.csv', header = T)
    #setInf = iris
 
    setInf$Blift = 0
 
    cutN = 3
    setInf[which(setInf$lift >= as.numeric(quantile(setInf$lift)[cutN])),]$Blift = 1
    
   
    setInf$type = setInf$Blift
    
    return(setInf)
}




#####KNN####

setInf = defSet(1)

library(kknn)
library(pROC)
library(modEvA)
index = sample(1:nrow(setInf), as.integer(0.7*nrow(setInf))); index; #70%随机抽样
table(setInf[index,ncol(setInf)])
table(setInf[,ncol(setInf)])


#参数设置
train = setInf[index, ] 
test = setInf[-index, ]
#test = read.csv('train.csv', header = T)
#
mol_KNN = kknn(type ~ OB + DisRwr + goSim + Jaccard + Man, train, test, k = 11 ,distance = 2 )
summary(mol_KNN)
#在测试集上预测
pre_kNN = fitted(mol_KNN)
pre_kNN2 = predict(mol_KNN, newdata = train)
pre_kNN

t = table(test$type, pre_kNN)
sum(diag(t))/sum(t) #准确率
#检验判类准确性（混淆矩阵
aKNN= table(test$type, pre_kNN, dnn =c("真实值","预测值"))
aKNN = as.data.frame(aKNN)

#绘制ROC曲线并计算AUC值

knn_roc = roc(test$type, as.numeric(pre_kNN),levels = c('0', '1'),direction = "<")
plot(knn_roc)
knn_roc$auc
coords(knn_roc, "best")
aupr=AUC(obs=test$type,pred=as.numeric(pre_kNN),curve = "ROC", simplif=TRUE, main = "ROC curve")
aupr=AUC(obs=test$type,pred=as.numeric(pre_kNN),curve = "PR", simplif=TRUE, main = "PR curve")

####KNN交叉验证#####
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
    pcl = fitted(kknn(type ~ OB + DisRwr + goSim + Jaccard + PMI, setInf[-omit,], setInf[omit,], k = 7 ,distance = 2 ))

},setInf)
#整合测试结果
wh = unlist(pred)
#grps 顺序已经被打乱，重新从小到大排序
setInf$grps = grps
knnT = as.data.frame(table(wh))
setInfF = setInf[order(setInf$grps),]

knn_roc = roc(setInfF$type, as.numeric(wh),levels = c('0', '1'),direction = "<")
plot(knn_roc)
knn_roc$auc
coords(knn_roc, "best")
aupr=AUC(obs=test$type,pred=as.numeric(pre_kNN),curve = "ROC", simplif=TRUE, main = "ROC curve")
aupr=AUC(obs=test$type,pred=as.numeric(pre_kNN),curve = "PR", simplif=TRUE, main = "PR curve")


#####SVM#####

setInf = defSet(1)

library(e1071)
library(pROC)
library(modEvA)

index = sample(1:nrow(setInf), as.integer(0.7*nrow(setInf))); index; #70%随机抽样
table(setInf[index,ncol(setInf)])
table(setInf[,ncol(setInf)])


SVMtrain = setInf[index, ] 
SVMtest = setInf[-index, ]
SVMtrain$type = as.factor(SVMtrain$type)
SVMtest$type = as.factor(SVMtest$type)

mol_SVM = svm(type ~ OB + DisRwr + goSim + Jaccard + PMI, data = SVMtrain, type = 'C-classification', kernel = 'radial')

#测试集
pre_svm = predict(mol_SVM, newdata = SVMtest)
obs_p_svm = data.frame(prob = pre_svm, obs = SVMtest$type)


#混淆矩阵
aSVM= table(SVMtest$type, pre_svm, dnn =c("真实值","预测值"))
aSVM = as.data.frame(aSVM)

#绘制ROC曲线并计算AUC值

SVM_roc = roc(SVMtest$type, as.numeric(pre_svm),levels = c('0', '1'),direction = "<")
plot(SVM_roc)
SVM_roc$auc
coords(SVM_roc, "best")
aupr=AUC(obs=SVMtest$type,pred=as.numeric(pre_svm),curve = "ROC", simplif=TRUE, main = "ROC curve")
aupr=AUC(obs=SVMtest$type,pred=as.numeric(pre_svm),curve = "PR", simplif=TRUE, main = "PR curve")

####SVM交叉验证#####
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
    pcl = predict(svm(type ~ OB + DisRwr + goSim + Jaccard + PMI, data = setInf[-omit,],  type = 'C-classification', kernel = 'radial'),newdata = setInf[omit,])
    
},setInf)
#整合测试结果
wh = unlist(pred)
#grps 顺序已经被打乱，重新从小到大排序
setInf$grps = grps
knnT = as.data.frame(table(wh))
setInfF = setInf[order(setInf$grps),]

#混淆矩阵
MLSVM= table(setInfF$type, wh, dnn =c("真实值","预测值"))
MLSVM = as.data.frame(MLSVM)

knn_roc = roc(setInfF$type, as.numeric(wh),levels = c('0', '1'),direction = "<")
plot(knn_roc)
knn_roc$auc
coords(knn_roc, "best")
aupr=AUC(obs=setInfF$type,pred=as.numeric(wh),curve = "ROC", simplif=TRUE, main = "ROC curve")
aupr=AUC(obs=setInfF$type,pred=as.numeric(wh),curve = "PR", simplif=TRUE, main = "PR curve")



#####GBRT#####

setInf = defSet(1)

library(gbm)

library(pROC)
library(modEvA)

index = sample(1:nrow(setInf), as.integer(0.7*nrow(setInf))); index; #70%随机抽样
table(setInf[index,ncol(setInf)])
table(setInf[,ncol(setInf)])


GBRTtrain = setInf[index, ] 
GBRTtest = setInf[-index, ]

mol_GBRT = gbm(type ~ OB + DisRwr + goSim + Jaccard + PMI,
    distribution = "bernoulli",
    data = GBRTtrain,
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

bestIter = gbm.perf(mol_GBRT,method = 'OOB')

pre_GBRT = predict(mol_GBRT, GBRTtest, bestIter)

#混淆矩阵
aGBRT= table(GBRTtest$type, pre_GBRT, dnn =c("真实值","预测值"))
aGBRT = as.data.frame(aGBRT)

#绘制ROC曲线并计算AUC值

GBRT_roc = roc(GBRTtest$type, as.numeric(pre_GBRT),levels = c('0', '1'),direction = "<")
plot(GBRT_roc)
GBRT_roc$auc
coords(GBRT_roc, "best")
aupr=AUC(obs=GBRTtest$type,pred=as.numeric(pre_GBRT),curve = "ROC", simplif=TRUE, main = "ROC curve")
aupr=AUC(obs=GBRTtest$type,pred=as.numeric(pre_GBRT),curve = "PR", simplif=TRUE, main = "PR curve")


####GBRT交叉验证#####
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
    mol_GBRTF = gbm(type ~ OB + DisRwr + goSim + Jaccard + PMI,
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

#混淆矩阵
MLGBRT= table(setInfF$type, wh, dnn =c("真实值","预测值"))
MLGBRT = as.data.frame(MLSVM)

knn_roc = roc(setInfF$type, as.numeric(wh),levels = c('0', '1'),direction = "<")
plot(knn_roc)
knn_roc$auc
coords(knn_roc, "best")
aupr=AUC(obs=setInfF$type,pred=as.numeric(wh),curve = "ROC", simplif=TRUE, main = "ROC curve")
aupr=AUC(obs=setInfF$type,pred=as.numeric(wh),curve = "PR", simplif=TRUE, main = "PR curve")




#####Bayes#####

setInf = defSet(1)

library(klaR)
library(e1071)
library(pROC)
library(modEvA)

index = sample(1:nrow(setInf), as.integer(0.7*nrow(setInf))); index; #70%随机抽样
table(setInf[index,ncol(setInf)])
table(setInf[,ncol(setInf)])


Bayestrain = setInf[index, ] 
Bayestest = setInf[-index, ]
Bayestrain$type = as.factor(Bayestrain$type)
Bayestest$type = as.factor(Bayestest$type)

mol_Bayes = NaiveBayes(type ~ OB + DisRwr + goSim + Jaccard + PMI,Bayestrain)
#mol_Bayes = naiveBayes(type ~ OB + HYRwr + CHDRwr + HYgoSim + CHDgoSim,Bayestrain)

pre_Bayes = predict(mol_Bayes, Bayestest)
pre_Bayes = as.data.frame(pre_Bayes)[,2]
#混淆矩阵
aBayes= table(Bayestest$type, pre_Bayes, dnn =c("真实值","预测值"))
aBayes = as.data.frame(aBayes)

#绘制ROC曲线并计算AUC值

Bayes_roc = roc(Bayestest$type, as.numeric(pre_Bayes),levels = c('1', '0'),direction = "<")
plot(Bayes_roc)
Bayes_roc$auc
coords(Bayes_roc, "best")
aupr=AUC(obs=Bayestest$type,pred=as.numeric(pre_Bayes),curve = "ROC", simplif=TRUE, main = "ROC curve")
aupr=AUC(obs=Bayestest$type,pred=as.numeric(pre_Bayes),curve = "PR", simplif=TRUE, main = "PR curve")



####Bayes交叉验证#####
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
    pcl = predict(NaiveBayes(type ~ OB + DisRwr + goSim + Jaccard + PMI,setInf[-omit,]),newdata = setInf[omit,])
    pre_BayesF = as.data.frame(pcl)[,2]
},setInf)
#整合测试结果
wh = unlist(pred)
#grps 顺序已经被打乱，重新从小到大排序
setInf$grps = grps
knnT = as.data.frame(table(wh))
setInfF = setInf[order(setInf$grps),]

#混淆矩阵
MLB= table(setInfF$type, wh, dnn =c("真实值","预测值"))
MLB = as.data.frame(MLSVM)

knn_roc = roc(setInfF$type, as.numeric(wh),levels = c('1', '0'),direction = "<")
plot(knn_roc)
knn_roc$auc
coords(knn_roc, "best")
aupr=AUC(obs=setInfF$type,pred=as.numeric(wh),curve = "ROC", simplif=TRUE, main = "ROC curve")
aupr=AUC(obs=setInfF$type,pred=as.numeric(wh),curve = "PR", simplif=TRUE, main = "PR curve")







#####LG#####

setInf = defSet(1)

library(pROC)
library(modEvA)

index = sample(1:nrow(setInf), as.integer(0.7*nrow(setInf))); index; #70%随机抽样
table(setInf[index,ncol(setInf)])
table(setInf[,ncol(setInf)])


LGtrain = setInf[index, ] 
LGtest = setInf[-index, ]
LGtrain$type = as.factor(LGtrain$type)
LGtest$type = as.factor(LGtest$type)

mol_LG = glm(type ~ OB + DisRwr + goSim + Jaccard + PMI,LGtrain,family=binomial(link="logit"))
summary(mol_LG)

pre_LG = predict.glm(mol_LG, type="response", newdata = LGtest)

#混淆矩阵
aLG= table(LGtest$type, pre_LG, dnn =c("真实值","预测值"))
aLG = as.data.frame(aLG)

#绘制ROC曲线并计算AUC值

LG_roc = roc(LGtest$type, as.numeric(pre_LG),levels = c('0', '1'),direction = "<")
plot(LG_roc)
LG_roc$auc
coords(LG_roc, "best")
aupr=AUC(obs=LGtest$type,pred=as.numeric(pre_LG),curve = "ROC", simplif=TRUE, main = "ROC curve")
aupr=AUC(obs=LGtest$type,pred=as.numeric(pre_LG),curve = "PR", simplif=TRUE, main = "PR curve")



####LG交叉验证#####
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
    pcl = predict.glm(glm(type ~ OB + DisRwr + goSim + Jaccard + PMI,setInf[-omit,],family=binomial(link="logit")),newdata = setInf[omit,])
},setInf)
#整合测试结果
wh = unlist(pred)
#grps 顺序已经被打乱，重新从小到大排序
setInf$grps = grps
knnT = as.data.frame(table(wh))
setInfF = setInf[order(setInf$grps),]

#混淆矩阵
MLLG= table(setInfF$type, wh, dnn =c("真实值","预测值"))
MLLG = as.data.frame(MLSVM)

knn_roc = roc(setInfF$type, as.numeric(wh),levels = c('0', '1'),direction = "<")
plot(knn_roc)
knn_roc$auc
coords(knn_roc, "best")
aupr=AUC(obs=setInfF$type,pred=as.numeric(wh),curve = "ROC", simplif=TRUE, main = "ROC curve")
aupr=AUC(obs=setInfF$type,pred=as.numeric(wh),curve = "PR", simplif=TRUE, main = "PR curve")














