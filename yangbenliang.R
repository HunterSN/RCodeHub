
rm(list = ls())

###队列研究设计 已知p1、p2，两样本量相等时四格表资料样本量估计

p0 = 0.58 #对照组
p1 = 0.86 #暴露组

p = (p0 + p1)/2

a = 0.05
b = 0.1


za = 1.65
zb = 1.28

n = (((za * (2 * p * (1 - p))^0.5) + (zb * (p1 * (1 - p1) + p0 * (1 - p0))^0.5))^2) /((p1-p0)^2)
