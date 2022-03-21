library(venn)


v = read.csv('targetbitch.csv',header = T)

hnv = v$hn
m2v = v$m2
m7v = v$m7
m0v = v$m0
m6v = v$m6
v = list(hnv,m2v,m7v,m0v,m6v)
venn(v,snames = "hnv,m2v,m7v,m0v,m6v", zcolor = "style", cexil = 1, cexsn = 0.8) #指定集合名称，设置默认的颜色风格，把交集区域数值标签字体改大到1，集合名称标签字体改成0.8



