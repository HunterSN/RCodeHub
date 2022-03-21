rm(list = ls())
zong = read.csv('rrgs.csv',header = T,stringsAsFactors = F)
yaowuHao = read.csv('yaohao.csv',header = T,stringsAsFactors = F)

sheTuan = as.data.frame(table(zong$分组))
ni = as.numeric(as.character(sheTuan[1,2]))
zhongYao = as.data.frame(table(zong$herb_num))
nj = as.numeric(as.character(zhongYao[1,2]))

nrni =as.numeric(nrow(sheTuan))
nrnj = as.numeric(nrow(zhongYao))
n = as.numeric(nrow(zong))
rrOut = data.frame()

y = 1
while (y <= nrnj) {
  njy = as.numeric(as.character(zhongYao[y,2]))
  
  s = 1
  while (s <= nrni) {
    nis = as.numeric(as.character(sheTuan[s,2]))
    zu2 = zong[which(zong$分组 == s-1),]
    nij = as.data.frame(table(zu2$herb_num))
    nijy = which(nij$Var1 == as.numeric(as.character(zhongYao[y,1])))
    nijyy = nij[nijy,2]
    if (length(nijy) ==0) {
      nijyy = 0
    }
    #nijyy = nij[nijy,2]
    rr = (nijyy/nis)/((njy - nijyy)/(n - nis))
    print(rr)
    hang = as.numeric((y-1)*nrni+s)
    print(hang)
    
    
    ka = nijyy
    kb = njy - nijyy  
    kc = nis - nijyy
    kd = (n - nis) - (njy - nijyy)
    tab = as.table(rbind(c(ka,kb),c(kc,kd)))
    xsq = chisq.test(tab)
    rrOut[hang,1] = s-1 #社团号
    rrOut[hang,2] = yaowuHao[which(yaowuHao$herb_num == as.numeric(as.character(zhongYao[y,1]))),1] #草药 
    rrOut[hang,3] = rr #rr
    rrOut[hang,4] = xsq$p.value #卡方p值
    rrOut[hang,5] = nijyy #nij 该社团内该药物数
    rrOut[hang,6] = nis #ni 该社团内总药物数
    rrOut[hang,7] = njy #nj 该药物总数
    
    s = s + 1
    #print(s)
  }
  y = y + 1
  #print(y)
}

names(rrOut) = c('社团号','草药','RR','卡方p值','该社团内该药物数','该社团内总药物数','该药物总数')
write.csv(rrOut,"RRout.csv",quote = F,row.names = F,fileEncoding = 'GBK')

rr1 = rrOut[which(rrOut$卡方p值 <= 0.05),]
rr2 = rr1[which(rr1$RR >1),]
write.csv(rr2,"RRout_shaixuan.csv",quote = F,row.names = F,fileEncoding = 'GBK')

i = 0
while (i <= 1) {
  rr3 = rr2[which(rr2$社团号 == i),]
  pzz = as.data.frame(quantile(rr3$该社团内该药物数))
  pz = pzz[3,]
  rrzz = as.data.frame(quantile(rr3$RR))
  rrz = rrzz[3,]
  
  rr4 = rr3[which(rr3$该社团内该药物数 >= pz),]
  rr5 = rr4[which(rr4$RR >= rrz),]
  rr5 = rr5[which(rr5$卡方p值 != 0),]
  shuchumingi = paste('RRout',i)
  shuchuming = paste(shuchumingi,".csv")
  write.csv(rr5,shuchuming,quote = F,row.names = F,fileEncoding = 'GBK')  
  i = i+1
}






rr3 = rr2[which(rr2$社团号 == 1),]
pzz = as.data.frame(quantile(rr3$该社团内该药物数))
pz = pzz[3,]
rrzz = as.data.frame(quantile(rr3$RR))
rrz = rrzz[3,]

rr4 = rr3[which(rr3$该社团内该药物数 >= pz),]
rr5 = rr4[which(rr4$RR >= rrz),]
write.csv(rr5,"RRout_1.csv",quote = F,row.names = F,fileEncoding = 'GBK')

