rm(list = ls())


a = read.csv('gsf_hx.csv',header = T)
attach(a)
table(herb)

par(family='STKaiti')
plot(herb)
plot(target)
