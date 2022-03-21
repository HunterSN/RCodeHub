rm(list = ls())



library(meta)#
library(robvis)
library(ggplot2)
fengxian = read.csv('fengxian.csv',header = T)
rob_summary(data = fengxian, tool = "ROB1",overall = F,weighted = FALSE, colour = c("skyblue","orange","green"))
rob_traffic_light(data = fengxian, tool = "ROB1", colour = c("skyblue","orange"), psize = 10)+ggtitle("")+theme(plot.title = element_text(hjust=0.5))


dataeff = read.csv('effective.csv',header = T)
dataadve = read.csv('adverse.csv',header = T)
datacrp = read.csv('CRP.csv',header = T)
datahscrp = read.csv('hs-CRP.csv',header = T)
datail6 = read.csv('IL-6.csv',header = T)
datatnfa = read.csv('TNF-α.csv',header = T)


#####连续变量#####
metalianxu = function(mdata, mSM){
  par(mfrow = c(2,1))
  mm<-metacont(n1,m1,sd1,n2,m2,sd2,  #试验组（暴露组） 对照组（非暴露组）
               
               data=mdata,sm= mSM,
               
               studlab = num)
  return(mm)
}


#####连续变量时间亚组#####


metalianxuitime = function(mdata, mSM, mtime){
  par(mfrow = c(2,1))
  
  mdata$Group = "not"
  mdata$Group[which(mdata$time >= mtime)] = "over"
  mdata$Group[which(mdata$time < mtime)] = "down"
  
  mm<-metacont(n1,m1,sd1,n2,m2,sd2,  #试验组（暴露组） 对照组（非暴露组）
               data=mdata,sm= mSM,
               studlab = num,
               byvar = Group,
               print.byvar = F)
  return(mm)
}

#####连续变量时间亚组egger#####



timeyazuegger = function(mdata, mSM,mtime,mmethod){
  mdata$Group = "not"
  mdata$Group[which(mdata$time >= mtime)] = "over"
  mdata$Group[which(mdata$time < mtime)] = "down"
  
  par(mfrow = c(length(unique(mdata$Group)),1))
  
  for (i in 1:length(unique(mdata$Group))) {
    mdataa = subset(mdata, mdata$Group == unique(mdata$Group)[i])
    mm<-metacont(n1,m1,sd1,n2,m2,sd2,  #试验组（暴露组） 对照组（非暴露组）
                 
                 data=mdataa,sm= mSM,
                 
                 studlab = num)
    
    em = metabias(mm,method = mmethod, plotit = T, k.min = nrow(mdataa) - 1)
    title(main = paste(unique(mdata$Group)[i],as.character(em$p.value)))
    print(paste(unique(mdata$Group)[i],as.character(em$p.value)))
  }
}

#####检测连续变量时间亚组egger#####


#####分类变量######
metafenlei = function(mdata, mSM){
  par(mfrow = c(2,1))
  mm<-metabin(t1,n1,c1,n2,  #试验组（暴露组） 对照组（非暴露组）
              
              data=mdata,sm= mSM,
              
              studlab = num)
  
  return(mm)
}





####egger######
egger = function(mm,mkmin, mmethod, mtitle){
  em = metabias(mm,method = mmethod, plotit = T, k.min = mkmin)
  title(main = mtitle)
  return(em)
}



#####有效性######
eff = metafenlei(dataeff, "OR") ##"OR", "RD", "RR", "ASD", or "DOR"
eff
forest(eff)
funnel(eff)
egger(eff, nrow(dataeff)-1, 'linreg', "effective")

######不良反应#####
adv = metafenlei(dataadve, 'OR')##"OR", "RD", "RR", "ASD", or "DOR"
adv
forest(adv)
funnel(adv)
egger(adv, nrow(dataadve)-1, 'linreg', "adverse")

########CRP#####
crp = metalianxu(datacrp, "MD") #"MD", "SM", "SMD", or "ROM".
crp
forest(crp)
funnel(crp)
egger(crp, nrow(datacrp)-1, 'rank', "CRP")



#######hs-CRP#####
hscrp = metalianxu(datahscrp, "MD")#"MD", "SM", "SMD", or "ROM".
hscrp
forest(hscrp)
funnel(hscrp)
egger(hscrp, nrow(datahscrp)-1, 'rank', "hs-CRP")


bijiao = unique(datahscrp$time)
########IL-6######
il6 = metalianxu(datail6, "MD")#"MD", "SM", "SMD", or "ROM".
il6
forest(il6)
funnel(il6)
egger(il6, nrow(datail6)-1, 'linreg', "IL-6")

########TNF-α####
tnfa = metalianxu(datatnfa, "MD")#"MD", "SM", "SMD", or "ROM"
tnfa
forest(tnfa)
funnel(tnfa)
egger(tnfa, nrow(datatnfa)-1, 'rank', "TNF-α")


par(mfrow = c(2,3))

egger(crp, nrow(datacrp)-1, 'rank', "CRP")
egger(hscrp, nrow(datahscrp)-1, 'rank', "hs-CRP")
egger(il6, nrow(datail6)-1, 'linreg', "IL-6")
egger(tnfa, nrow(datatnfa)-1, 'rank', "TNF-alpha")
egger(adv, nrow(dataadve)-1, 'linreg', "Adverse reactions")

