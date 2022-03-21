rm(list = ls())



a = read.csv('gsfsym.csv',header =  T ,stringsAsFactors = F)
b = split(a$targe_num,a$source)  #按处方分类组成list

c = as.data.frame(t(sapply(b, "[", i = 1:max(sapply(b, length)))))
c = cbind(source = row.names(c), c)

library(openxlsx)
write.xlsx(c, file = "lzh.xlsx", colNames = TRUE)
