rm(list = ls())



a = read.csv('gsfd.csv',header =  T ,stringsAsFactors = F)
b = split(a$targe_num,a$source)  #按处方分类组成list
num = unique(as.vector(a$source)) #提取处方号

zhongjie = as.data.frame(c())
for (i in 2:length(num)) {

  jiao = length(intersect(as.vector(b[[1]]),as.vector(b[[i]])))
  bing = length(as.vector(b[[1]]))
  #bing = length(union(as.vector(b[[1]]),as.vector(b[[i]])))
  jaccardnum = jiao/bing 
 
  jieguo = c(1,num[i],jaccardnum)
  
  
  hebing = rbind(zhongjie,jieguo)
  zhongjie = hebing
  
}

shuchujacweight = na.omit(zhongjie)

names(shuchujacweight) = c("source","target","weight")

shuchujacweight$fullper = shuchujacweight$weight * 100

cut = table(cut(shuchujacweight$fullper,breaks = c(0,70,100)))
cut[2]/(cut[1] + cut[2])

jishu = table(cut(shuchujacweight$fullper,breaks = c(0,10,20,30,40,50,60,70,80,90,100)))

cc = as.data.frame(jishu)
barplot(jishu)



library(openxlsx)
write.xlsx(shuchujacweight, file = "yGSFquanzhong.xlsx", colNames = TRUE)


