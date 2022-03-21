rm(list = ls())


library(Matrix)
library(arules) #apriori
library(arulesViz)



a = read.csv('m1.csv',header = T,stringsAsFactors = F)#读取数据
b_apriori = split(a$herb,a$num)#按处方分类组成list
c_apriori = apriori(as(b_apriori,'transactions'),parameter = list(supp = 0.1,conf = 0.5,maxlen = 2,minlen = 2)) #apriori计算
summary(c_apriori)
output =inspect(c_apriori)

library(openxlsx)
#write.xlsx(output, file = "output.xlsx", colNames = TRUE)


#A1 = subset(c_apriori, items %in% c('banxia','sanqi','chenpi'))

####图######
plot(c_apriori, method = "grouped",control = list(k = 30), main = "")
plot(c_apriori, method = "graph")
plot(c_apriori, method = "graph", shading = 'lift')
plot(c_apriori)
plot(c_apriori, measure = c("support", "lift"), shading = "confidence")
#plot(c_apriori, interactive = T)

plot(c_apriori,shading = 'order',control = list(main = 'Two-key plot'))




#write.csv(output,"ding_gxb_supp_0.2_conf_0.8_apriori.csv",quote = F,row.names = F,fileEncoding = 'GBK')
#write.csv(output,"ding_gxb_apriori.csv",quote = F,row.names = F,fileEncoding = 'GBK')




#((Coronary Intervention*, Percutaneous[Title/Abstract]) OR (Revascularization*, Percutaneous Coronary [Title/Abstract]) OR (Intervention*, Percutaneous Coronary [Title/Abstract]) OR (Percutaneous Coronary Interventions [Title/Abstract]) OR (Percutaneous Coronary Revascularization* [Title/Abstract]) OR (Coronary Revascularization*, Percutaneous [Title/Abstract]))


#######总#####
c_apriori = apriori(as(b_apriori,'transactions'),parameter = list(supp = 0,conf = 0,maxlen = 2,minlen = 2)) #apriori计算
summary(c_apriori)
output =inspect(c_apriori)
library(openxlsx)
write.xlsx(output, file = "outputzong.xlsx", colNames = TRUE)

####总taps筛选######

rm(list = ls())


library(Matrix)
library(arules) #apriori
library(arulesViz)



a = read.csv('m1.csv',header = T,stringsAsFactors = F)#读取数据
b_apriori = split(a$herb,a$num)#按处方分类组成list
c_apriori = apriori(as(b_apriori,'transactions'),parameter = list(supp = 0,conf = 0,maxlen = 2,minlen = 2)) #apriori计算
summary(c_apriori)
output =inspect(c_apriori)




shai = output
shai$lhs =  gsub('[{]',"",shai$lhs)
shai$lhs =  gsub('[}]',"",shai$lhs)

shai$rhs =  gsub('[{]',"",shai$rhs)
shai$rhs =  gsub('[}]',"",shai$rhs)
library(openxlsx)
write.xlsx(shai, file = "outputzong.xlsx", colNames = TRUE)



taps = read.csv('taps.csv',header = F,stringsAsFactors = F)



aa = subset(shai, (is.element(shai$lhs,taps[,1])) & (is.element(shai$rhs,taps[,1])))
library(openxlsx)
write.xlsx(aa, file = "tapsapriori.xlsx", colNames = TRUE)

####筛taps筛选######



c_apriori = apriori(as(b_apriori,'transactions'),parameter = list(supp = 0.2,conf = 0.5,maxlen = 2,minlen = 2)) #apriori计算
summary(c_apriori)
output =inspect(c_apriori)




shai = output
shai$lhs =  gsub('[{]',"",shai$lhs)
shai$lhs =  gsub('[}]',"",shai$lhs)

shai$rhs =  gsub('[{]',"",shai$rhs)
shai$rhs =  gsub('[}]',"",shai$rhs)
library(openxlsx)
write.xlsx(shai, file = "outputzong.xlsx", colNames = TRUE)

taps = read.csv('taps.csv',header = F,stringsAsFactors = F)



aa = subset(shai, (is.element(shai$lhs,taps[,1])) & (is.element(shai$rhs,taps[,1])))
library(openxlsx)
write.xlsx(aa, file = "shaitapsapriori.xlsx", colNames = TRUE)


