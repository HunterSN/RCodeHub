rm(list = ls())


filename <- dir(pattern=".csv")
hebing = c()
for (i in 1:length(filename)){
  hebing[i] <- gsub('?.csv','',filename[i])
  assign(hebing[i],read.csv(filename[i], header = T,stringsAsFactors = F))
}

clu = 1
TotalOut = read.csv(file = filename[clu], header = T,stringsAsFactors = F)
jishu = nrow(TotalOut)
for (clu in 2:length(filename)) {
  print(clu)
  newData = read.csv(file = filename[clu], header = T,stringsAsFactors = F)
  jishu = jishu + nrow(newData)
  print(jishu)
  TotalOut = rbind(TotalOut,newData)
  print(nrow(TotalOut))
}
print(nrow(TotalOut))
rm(newData)
write.csv(TotalOut,"TotalOut.csv",quote = F,row.names = F,fileEncoding = 'GBK')



