rm(list = ls())
library(UpSetR)

#####通用函数#####
tongyongHCT = function(tongyong, nset, sets,xlabname){
  tongyong = tongyong[which(rowSums(tongyong[-1]) > 0),]
  upset(fromList(tongyong), 
        nsets = nset, 
        sets = sets, 
        mainbar.y.label = "intersection size", 
        sets.x.label = xlabname,#"Numbers of **** enrichment gene"
        main.bar.color = "#2a83a2", sets.bar.color = "#3b7960",
        mb.ratio = c(0.5, 0.5),
        order.by = "freq", 
        decreasing = c(TRUE,FALSE),
        queries = list(list(query = intersects, params = list(c(names(kegg)[-1])), color = '#2f5688',active = T))
  )
}


######单个#####
kegg = read.csv('gsfkegg.csv',header = T)
tongyongHCT(kegg, nrow(kegg)-1, names(kegg)[-1],"Numbers of KEGG enrichment gene")

go = read.csv('gsfgo.csv',header = T)
tongyongHCT(go, nrow(go)-1, names(go)[-1],"Numbers of GO enrichment gene")

mol = read.csv('mol.csv',header = T)
tongyongHCT(mol, nrow(mol)-1, names(mol)[-1],"Numbers of MOL enrichment gene")

