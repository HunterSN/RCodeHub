rm(list = ls())


library(pander)
library(dplyr)
library(Epi)
whickham1 <- matrix(c(519, 980, 200, 356), byrow=TRUE, nrow=2)
rownames(whickham1) <- c("TCM", "Non-TCM")
colnames(whickham1) <- c("Fmale", "Male")
pander(addmargins(whickham1))
twoby2(whickham1)
chisq.test(whickham1)

age <- c(rep("18-24", 4), rep("25-34", 4), 
         rep("35-44", 4), rep("45-54", 4), 
         rep("55-64", 4), rep("65-74", 4), 
         rep("75+", 4))
treat <- c(rep(c("TCM", "TCM", "Non-TCM", "Non-TCM"), 7))
sex <- c(rep(c("Fmale", "Male"), 14))
counts <- c(1, 4, 0, 2, 10, 65, 5, 16, 31, 73, 4, 25, 73, 103, 14, 34, 69, 144, 17, 56, 126, 172, 89, 42, 209, 419, 71, 181)
whickham2 <- data.frame(treat, sex, age, counts) %>% tbl_df()
whickham2$treat <- factor(whickham2$treat, levels = c("TCM", "Non-TCM"))
whickham2.tab1 <- xtabs(counts ~ treat + sex + age, data = whickham2)
whickham2.tab1
ftable(whickham2.tab1)
mantelhaen.test(whickham2.tab1)
