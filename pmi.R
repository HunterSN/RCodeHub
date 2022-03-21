rm(list = ls())

a = read.csv('m1.csv',header = T,stringsAsFactors = F)  #读取数据
b = split(a$herb,a$num)  #按处方分类组成list
herb = unique(as.vector(a$herb)) #提取草药
ch = t(combn(herb,2)) 
hangshu = nrow(ch)
caoyaohuxinxi = data.frame()
nrow(ch)
chufang = length(b)

#正式代码
for (i in 1:nrow(ch)) {
  print(i)
  countherb1 = 0
  countherb2 = 0
  tongshi = 0
  for (m in b) {
    if (is.element(ch[i,1],m)) {
      countherb1 = countherb1 + 1
    }
  }
  for (m in b) {
    if (is.element(ch[i,2],m)) {
      countherb2 = countherb2 + 1
    }
  }
  for (m in b) {
    if (is.element(ch[i,1],m)&&is.element(ch[i,2],m)) {
      tongshi = tongshi + 1
    }
  }
  huxinxi = abs(log((tongshi/chufang)/((countherb1/chufang)*(countherb2/chufang))))
  xx = ch[i,1]
  yy = ch[i,2]
  caoyaohuxinxi[i,1] = xx
  caoyaohuxinxi[i,2] = yy
  caoyaohuxinxi[i,3] = countherb1
  caoyaohuxinxi[i,4] = countherb2
  caoyaohuxinxi[i,5] = tongshi
  caoyaohuxinxi[i,6] = huxinxi
}

names(caoyaohuxinxi) = c("草药A","草药B","A频数","B频数","共用频数","PMI")

write.csv(caoyaohuxinxi,"m1pmi.csv",quote = F,row.names = F,fileEncoding = 'GBK')



