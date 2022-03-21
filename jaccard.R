rm(list = ls())



a = read.csv('zgxjb.csv',header =  T ,stringsAsFactors = F)
b = split(a$targe_num,a$source)  #按处方分类组成list
num = unique(as.vector(a$source)) #提取处方号
quanzhong = t(combn(num,2)) 
r = nrow(quanzhong)
o = 0
aa = c()
while (o <= r - 2000001) {
  o = o + 2000000
  bb = c(o)
  aa = c(aa,bb)
  
}
aa = c(aa,r)

m = 1
for (l in aa) {
  print(l)
  ceshi1 = quanzhong[m:l,]
  jacweight = data.frame()
  zhongjie = c(0,0,0)
  for (i in 1:nrow(ceshi1)) {
    print(i)
    xx = ceshi1[i,1]
    yy = ceshi1[i,2]
    x = as.character(xx)
    y = as.character(yy)
    jiao = length(intersect(as.vector(b[[x]]),as.vector(b[[y]])))
    bing = length(union(as.vector(b[[x]]),as.vector(b[[y]])))
    jaccardnum = jiao/bing 
    
    if (jaccardnum>=0.5) {
      jacweight[i,1] = xx
      jacweight[i,2] = yy
      jacweight[i,3] = jaccardnum
    }
    
    jieguo = c(xx,yy,jaccardnum)
    
    
    hebing = rbind(zhongjie,jieguo)
    zhongjie = hebing
    
  }
  
  shuchujacweight = na.omit(jacweight)
  
  names(shuchujacweight) = c("source","target","weight")
  write.csv(shuchujacweight,paste(as.character(l),"quanzhong.csv") ,quote = F,row.names = F,fileEncoding = 'GBK')
  
  write.csv(zhongjie,paste(as.character(l),'zong.csv'),quote = F,row.names = F,fileEncoding = 'GBK')
  
  rm(hebing)
  
  rm(shuchujacweight)
  rm(jacweight)
  rm(zhongjie)
  
  
  
  
  m = l + 1
  
}

#write.table(shuchujacweight,paste(as.character(l),"quanzhong2.txt") ,sep = "\t",quote = F,row.names = F,fileEncoding = 'GBK')

#write.table(zhongjie,paste(as.character(l),'zong2.txt'),sep = "\t",quote = F,row.names = F,fileEncoding = 'GBK')
