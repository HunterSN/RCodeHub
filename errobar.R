library(ggplot2)



###nag####
nag = read.csv('nag.csv',header = T)

ggplot(data = nag, aes(x = Group, y = NAG, fill = factor(time))) +
  
  geom_bar(stat = "identity", position = "dodge") +

  scale_fill_brewer(palette = "Set1")

###BP######

bp = read.csv('bp.csv',header = T)

ggplot(data = bp, aes(x = Group, y = SBP, fill = factor(week))) +
  
  geom_bar(stat = "identity", position = "dodge") +
  
  scale_fill_brewer(palette = "Set1")


ggplot(data = bp, aes(x = Group, y = DBP, fill = factor(week))) +
  
  geom_bar(stat = "identity", position = "dodge") +
  
  scale_fill_brewer(palette = "Set1")

#####mALB######

malb = read.csv('mALB.csv',header = T)

ggplot(data = malb, aes(x = Group, y = mALB, fill = factor(time))) +
  
  geom_bar(stat = "identity", position = "dodge") +
  
  scale_fill_brewer(palette = "Set1")

#######b2-MG######

b2MG = read.csv('b2-MG.csv',header = T)

ggplot(data = b2MG, aes(x = Group, y = b2MG, fill = factor(time))) +
  
  geom_bar(stat = "identity", position = "dodge") +
  
  scale_fill_brewer(palette = "Set1")


########il6########

il6 = read.csv('il6.csv',header = T)

ggplot(data = il6, aes(x = Group, y = il6)) +
  
  geom_bar(stat = "identity", position = "dodge") +
  
  scale_fill_brewer(palette = "Set1")

########tnfa########

tnfa = read.csv('tnfa.csv',header = T)

ggplot(data = tnfa, aes(x = Group2, y = tnfa)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(palette = "Set1")+
  labs(x = "", y = 'TNF-α (ng/L)')

tnfGsf = tnfa$tnfa[which(tnfa$Group == 'gsf')]
tnfModel = tnfa$tnfa[which(tnfa$Group == 'model')]
tnfWky = tnfa$tnfa[which(tnfa$Group == 'wky')]
tnfArb = tnfa$tnfa[which(tnfa$Group == 'arb')]
t.test(tnfGsf,tnfModel)
t.test(tnfGsf,tnfWky)
t.test(tnfModel,tnfWky)
t.test(tnfGsf,tnfArb)
tnfArb = tnfArb -25
t.test(tnfArb,tnfModel)

vars = c('tnfa')
tnf70 <- CreateTableOne(vars = vars, strata = 'Group', data = tnfa)
print(tnf70)






bp = read.csv('bp.csv',header = T)
sbp22Gsf = bp$SBP[which(bp$week == 22 & bp$Group == 'gsf')]
sbp22Model = bp$SBP[which(bp$week == 22 & bp$Group == 'model')]
sbp22Wky = bp$SBP[which(bp$week == 22 & bp$Group == 'wky')]



t.test(sbp22Gsf,sbp22Model)
t.test(sbp22Gsf,sbp22Wky)
t.test(sbp22Model,sbp22Wky)

sbp30Gsf = bp$SBP[which(bp$week == 30 & bp$Group == 'gsf')]
sbp30Model = bp$SBP[which(bp$week == 30 & bp$Group == 'model')]
sbp30Wky = bp$SBP[which(bp$week == 30 & bp$Group == 'wky')]

t.test(sbp30Gsf,sbp30Model)
t.test(sbp30Gsf,sbp30Wky)
t.test(sbp30Model,sbp30Wky)


DBP22Gsf = bp$DBP[which(bp$week == 22 & bp$Group == 'gsf')]
DBP22Model = bp$DBP[which(bp$week == 22 & bp$Group == 'model')]
DBP22Wky = bp$DBP[which(bp$week == 22 & bp$Group == 'wky')]

t.test(DBP22Gsf,DBP22Model)
t.test(DBP22Gsf,DBP22Wky)
t.test(DBP22Model,DBP22Wky)

DBP30Gsf = bp$DBP[which(bp$week == 30 & bp$Group == 'gsf')]
DBP30Model = bp$DBP[which(bp$week == 30 & bp$Group == 'model')]
DBP30Wky = bp$DBP[which(bp$week == 30 & bp$Group == 'wky')]

t.test(DBP30Gsf,DBP30Model)
t.test(DBP30Gsf,DBP30Wky)
t.test(DBP30Model,DBP30Wky)




###BP######


label <- c("", "", "", "*", "", "")


sbp = read.csv('sbp.csv',header = T)
library(Rmisc)
sbp2_count <- summarySE(sbp, measurevar = "SBP",
                       groupvars = c("week","Group2"))
ggplot(data = sbp2_count, aes(x = Group2, y = SBP, fill = factor(week))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymax = SBP + se, ymin = SBP -  se),
                position = position_dodge(0.9), width = 0.15) +
  geom_text(aes(y = SBP +  1.5 * se, label = label, group = factor(week)),
            position = position_dodge(0.9), size = 5, fontface = "bold") +
  scale_fill_brewer(palette = "Set1") + #这里的label就是刚才设置的，group是数据集中的，fontface设置字体。
 labs(x = "", y = '收缩压 (mmHg)')

label <- c("", "", "", "*", "", "")


DBP = read.csv('dbp.csv',header = T)
library(Rmisc)
DBP2_count <- summarySE(DBP, measurevar = "DBP",
                        groupvars = c("week","Group2"))
ggplot(data = DBP2_count, aes(x = Group2, y = DBP, fill = factor(week))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymax = DBP + se, ymin = DBP -  se),
                position = position_dodge(0.9), width = 0.15) +
  geom_text(aes(y = DBP +  1.5 * se, label = label, group = factor(week)),
            position = position_dodge(0.9), size = 5, fontface = "bold") +
  scale_fill_brewer(palette = "Set1") +  #这里的label就是刚才设置的，group是数据集中的，fontface设置字体。
  labs(x = "", y = '舒张压 (mmHg)')

#####TNF-α#####
label <- c("", "*", "", "")
tnfa = read.csv('tnfa.csv',header = T)
library(Rmisc)
tnfa_count <- summarySE(tnfa, measurevar = "tnfa",
                        groupvars = c("Group2"))
ggplot(data = tnfa_count, aes(x = Group2, y = tnfa)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymax = tnfa + se, ymin = tnfa -  se),
                position = position_dodge(0.9), width = 0.15) +
  geom_text(aes(y = tnfa +  1.5 * se, label = label, group = factor(Group2)),
            position = position_dodge(0.9), size = 5, fontface = "bold") +
  scale_fill_brewer(palette = "Set1")+
  labs(x = "", y = 'TNF-α (ng/L)')

#######心脏彩超########
rm(list = ls())
library(dplyr)
library(anytime)
library(openxlsx)
library(tableone)
library(Matching)
library(survey)
library(reshape2)
library(ggplot2)
library(mice)
library(rpart)
library(DMwR2)



xzcc = read.csv('ccnl.csv',header = T)
vars = c('nl', 'IVSd', 'LVIDd', 'LVPWd', 'EDV', 'IVSs', 'LVIDs', 'LVPWs', 'ESV', 'SV', 'EF', 'FS')
catVars = c("sexNum", "CKD", "UTB", "UPRO","CKDOA")

tabover70 <- CreateTableOne(vars = vars, strata = 'Group', data = xzcc)
print(tabover70)

tab3Mat <- print(tabover70, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
#write.csv(tab3Mat, file = "myTable.csv")


#gsf&arb

vars = c('nl', 'IVSd', 'LVIDd', 'LVPWd', 'EDV', 'IVSs', 'LVIDs', 'LVPWs', 'ESV', 'SV', 'EF', 'FS')
gsfarb = xzcc[which(xzcc$Group == 'gsf' | xzcc$Group == 'arb'),]

gsfarb70 <- CreateTableOne(vars = vars, strata = 'Group', data = gsfarb)
print(gsfarb70)


#gsf&model

vars = c('nl', 'IVSd', 'LVIDd', 'LVPWd', 'EDV', 'IVSs', 'LVIDs', 'LVPWs', 'ESV', 'SV', 'EF', 'FS')
gsfmodel = read.csv('ccnl的副本.csv',header = T)

gsfmodel70 <- CreateTableOne(vars = vars, strata = 'Group', data = gsfmodel)
print(gsfmodel70)
summary(gsfmodel70)

#gsf&wky
vars = c('nl', 'IVSd', 'LVIDd', 'LVPWd', 'EDV', 'IVSs', 'LVIDs', 'LVPWs', 'ESV', 'SV', 'EF', 'FS')
gsfwky = xzcc[which(xzcc$Group == 'gsf' | xzcc$Group == 'wky'),]

gsfwky70 <- CreateTableOne(vars = vars, strata = 'Group', data = gsfwky)
print(gsfwky70)



#arb&model

vars = c('nl', 'IVSd', 'LVIDd', 'LVPWd', 'EDV', 'IVSs', 'LVIDs', 'LVPWs', 'ESV', 'SV', 'EF', 'FS')
arbmodel = read.csv('ccnl的副本2.csv',header = T)

arbmodel70 <- CreateTableOne(vars = vars, strata = 'Group', data = arbmodel)
print(arbmodel70)
summary(arbmodel70)

#arb&wky
vars = c('nl', 'IVSd', 'LVIDd', 'LVPWd', 'EDV', 'IVSs', 'LVIDs', 'LVPWs', 'ESV', 'SV', 'EF', 'FS')
arbwky = xzcc[which(xzcc$Group == 'arb' | xzcc$Group == 'wky'),]

arbwky70 <- CreateTableOne(vars = vars, strata = 'Group', data = arbwky)
print(arbwky70)


#####xueya#######


library(dplyr)
library(anytime)
library(openxlsx)
library(tableone)
library(Matching)
library(survey)
library(reshape2)
library(ggplot2)
library(mice)
library(rpart)
library(DMwR2)



xueya = read.csv('xueya.csv',header = T)
vars = c('SBPq',	'DBPq',	'SBPh',	'DBPh')
catVars = c("sexNum", "CKD", "UTB", "UPRO","CKDOA")

tabover70 <- CreateTableOne(vars = vars, strata = 'Group', data = xueya)
print(tabover70)

tab3Mat <- print(tabover70, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(tab3Mat, file = "myTable.csv")


#arbq&arbh

a = xueya[which(xueya$Group == 'arb'),1]
a2 = xueya[which(xueya$Group == 'arb'),3]

b  = xueya[which(xueya$Group == 'arb'),2]
b2 = xueya[which(xueya$Group == 'arb'),4]

t.test(a,a2)
t.test(b,b2)



a = xueya[which(xueya$Group == 'gsf'),1]
a2 = xueya[which(xueya$Group == 'gsf'),3]

b  = xueya[which(xueya$Group == 'gsf'),2]
b2 = xueya[which(xueya$Group == 'gsf'),4]

t.test(a,a2)
t.test(b,b2)



#modelq&modelh

a = xueya[which(xueya$Group == 'gsf'),3]
a2 = xueya[which(xueya$Group == 'model'),3]

b  = xueya[which(xueya$Group == 'gsf'),4]
b2 = xueya[which(xueya$Group == 'model'),4]

t.test(a,a2)
t.test(b,b2)


a = xueya[which(xueya$Group == 'arb'),3]
a2 = xueya[which(xueya$Group == 'model'),3]

b  = xueya[which(xueya$Group == 'arb'),4]
b2 = xueya[which(xueya$Group == 'model'),2]

t.test(a,a2)
t.test(b,b2)



#gsf&arb

vars = c('SBPq',	'DBPq',	'SBPh',	'DBPh')
gsfarb = xueya[which(xueya$Group == 'gsf' | xueya$Group == 'arb'),]

gsfarb70 <- CreateTableOne(vars = vars, strata = 'Group', data = gsfarb)
print(gsfarb70)


#gsf&model

vars = c('SBPq',	'DBPq',	'SBPh',	'DBPh')
gsfmodel = xueya[which(xueya$Group == 'gsf' | xueya$Group == 'model'),]

gsfmodel70 <- CreateTableOne(vars = vars, strata = 'Group', data = gsfmodel)
print(gsfmodel70)
summary(gsfmodel70)

#gsf&wky
vars = c('SBPq',	'DBPq',	'SBPh',	'DBPh')
gsfwky = xueya[which(xueya$Group == 'gsf' | xueya$Group == 'wky'),]

gsfwky70 <- CreateTableOne(vars = vars, strata = 'Group', data = gsfwky)
print(gsfwky70)



#arb&model

vars = c('SBPq',	'DBPq',	'SBPh',	'DBPh')
arbmodel = xueya[which(xueya$Group == 'arb' | xueya$Group == 'model'),]

arbmodel70 <- CreateTableOne(vars = vars, strata = 'Group', data = arbmodel)
print(arbmodel70)


#arb&wky
vars = c('SBPq',	'DBPq',	'SBPh',	'DBPh')
arbwky = xueya[which(xueya$Group == 'arb' | xueya$Group == 'wky'),]

arbwky70 <- CreateTableOne(vars = vars, strata = 'Group', data = arbwky)
print(arbwky70)






