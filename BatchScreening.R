rm(list = ls())


filename <- dir(pattern=".csv")
disnet = c()
for (i in 1:length(filename)){
  disnet[i] <- gsub('?.csv','',filename[i])
  assign(disnet[i],read.csv(filename[i], header = T))
}

jiezhi = 3
bianhao = paste0(paste0("degree>=",jiezhi),".csv")
catch = data.frame()
for (clu in 1:length(filename)) {
  print(clu)
  newData = read.csv(file = filename[clu], header = T)
  outdata = newData[which(newData$Degree >= jiezhi),]
  outname = filename[clu]
  write.csv(outdata,sub(".csv",bianhao,filename[clu]),quote = F,row.names = F,fileEncoding = 'GBK')
  
  outdata$class = sub(".csv",bianhao,filename[clu])
  catch = rbind(catch,outdata)
  
}

  write.csv(catch,'zong.csv',quote = F,row.names = F,fileEncoding = 'GBK')


