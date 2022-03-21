rm(list = ls())



library(meta)#
library(robvis)

dataeff = read.csv('effective.csv',header = T)
dataadve = read.csv('adverse.csv',header = T)
datacrp = read.csv('CRP.csv',header = T)
datahscrp = read.csv('hs-CRP.csv',header = T)
datail6 = read.csv('IL-6.csv',header = T)
datatnfa = read.csv('TNF-α.csv',header = T)
datail10 = read.csv('IL-10.csv',header = T)

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
eff = metafenlei(dataeff, "RD") ##"OR", "RD", "RR", "ASD", or "DOR"
eff
forest(eff)
funnel(eff)
egger(eff, nrow(dataeff)-1, 'linreg', "effective")

######不良反应#####
adv = metafenlei(dataadve, 'RD')##"OR", "RD", "RR", "ASD", or "DOR"
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

######CRP亚组都不好####
crpy = metalianxuitime(datacrp, "MD", 16)
forest(crpy)

forest(crpy)
crpy = metalianxuitime(datacrp, "MD", 24)
forest(crpy)
forest(crpy)
funnel(crpy)
bijiao = unique(datacrp$time)


timeyazuegger(datatnfa, "MD", 16,  'linreg')###down可

#######hs-CRP#####
hscrp = metalianxu(datahscrp, "SMD")#"MD", "SM", "SMD", or "ROM".
hscrp
forest(hscrp)
funnel(hscrp)
egger(hscrp, nrow(datahscrp)-1, 'rank', "hs-CRP")
######hs-CRP亚组####
hscrpy = metalianxuitime(datahscrp, "SMD", 48)####
forest(hscrpy)
hscrpy = metalianxuitime(datahscrp, "SMD", 8)####
forest(hscrpy)


bijiao = unique(datahscrp$time)
########IL-6######
il6 = metalianxu(datail6, "MD")#"MD", "SM", "SMD", or "ROM".
il6
forest(il6)
funnel(il6)
egger(il6, nrow(datail6)-1, 'linreg', "IL-6")
######IL-6亚组均不好####
bijiao = unique(datail6$time)
il6y = metalianxuitime(datail6, "SMD", 48)
forest(il6y)
il6y = metalianxuitime(datail6, "SMD", 4)
forest(il6y)
il6y = metalianxuitime(datail6, "SMD", 24)
forest(il6y)
il6y = metalianxuitime(datail6, "SMD", 2)
forest(il6y)

il6y = metalianxuitime(datail6, "SMD", 12)
forest(il6y)
il6y = metalianxuitime(datail6, "SMD", 52)
forest(il6y)
il6y = metalianxuitime(datail6, "SMD", 16)
forest(il6y)
il6y = metalianxuitime(datail6, "SMD", 8)
forest(il6y)

il6y = metalianxuitime(datail6, "MD", 48)
forest(il6y)
il6y = metalianxuitime(datail6, "MD", 4)
forest(il6y)
il6y = metalianxuitime(datail6, "MD", 24)
forest(il6y)
il6y = metalianxuitime(datail6, "MD", 2)
forest(il6y)

il6y = metalianxuitime(datail6, "MD", 12)
forest(il6y)
il6y = metalianxuitime(datail6, "MD", 52)
forest(il6y)
il6y = metalianxuitime(datail6, "MD", 16)
forest(il6y)
il6y = metalianxuitime(datail6, "MD", 8)
forest(il6y)

########TNF-α####
tnfa = metalianxu(datatnfa, "MD")#"MD", "SM", "SMD", or "ROM"
tnfa
forest(tnfa)
funnel(tnfa)
egger(tnfa, nrow(datatnfa)-1, 'rank', "TNF-α")

########TNF-α亚组####

tnfay = metalianxuitime(datatnfa, "MD", 16)#"MD", "SM", "SMD", or "ROM"
forest(tnfay)
tnfay = metalianxuitime(datatnfa, "MD", 24)#######
forest(tnfay)
tnfay = metalianxuitime(datatnfa, "MD", 12)
forest(tnfay)
tnfay = metalianxuitime(datatnfa, "SMD", 24)
forest(tnfay)
tnfay
forest(tnfay)
funnel(tnfay)
bijiao = unique(datatnfa$time)


timeyazuegger(datatnfa, "MD", 16,  'linreg')###down可

timeyazuegger(datatnfa, "MD", 24,  'rank')#可#####

timeyazuegger(datatnfa, "MD", 12,  'rank')#可

timeyazuegger(datatnfa, "MD", 16,  'rank')#可

timeyazuegger(datatnfa, "SMD", 24,  'linreg')#单可

timeyazuegger(datatnfa, "SMD", 24,  'rank')#单可#


#######IL-10#####
il10 = metalianxu(datail10, "MD")#"MD", "SM", "SMD", or "ROM"
il10
forest(il10)
funnel(il10)
egger(il10, nrow(datail10)-1, 'linreg', "IL-10")

######IL-10亚组均不好####

