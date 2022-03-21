rm(list = ls())

a = read.csv('dbn.csv',header = F,stringsAsFactors = F)

lei = as.numeric(ncol(a))
hang = as.numeric(nrow(a))

b= c()
b = as.data.frame(b)


zs = as.numeric(lei*hang)



for (i in 1:hang) {
  for (ii in 1:(lei-1)) {
    xu = as.numeric(as.numeric(nrow(b)) + 1)
    b[xu,1] = a[i,1]
    b[xu,2] = a[i,ii+1]
  }
}



bb = subset(b,b$V2 != "")

write.csv(bb,"dbnhzl.csv",quote = F,row.names = F,fileEncoding = 'GBK')
#write.table(x, file = "pinggu.txt", row.names = FALSE, col.names = FALSE)


#转单列
de = read.csv('qita.csv',header = F,stringsAsFactors = F)
data = as.matrix(de)
nrow(data)
dim(data) = c(nrow(data)*ncol(data),1)
data= subset(data,data != "")
write.csv(data,"qitade.csv",quote = F,row.names = F,fileEncoding = 'GBK')






####总表

rm(list = ls())
library(tidyverse)
separate(data = test, col = pedigree, into = c("m", "trait1"), sep = "_X_")

z = read.csv('zongbiao.csv',header = T,stringsAsFactors = F)

lieshu = ncol(z)
hanghshu = nrow(z)
#分割列
separate(data = z, col = z[,10], into = c(), sep = ",")
