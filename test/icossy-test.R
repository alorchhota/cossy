source("./icossy/icossy.R")
source("./icossy/data-preparation.R")



gctfile <- "data/cns.gct"
clsfile <- "data/cns.cls"


icsy <- icossy(gctfile = gctfile,
               chipfile = NA, 
               clsfile = clsfile, 
               network = "pathwayapi", 
               nmis = 15, 
               frank = T, 
               qnorm = F, 
               ztrans = F,
               sig.test = "ttest")
print(icsy$status)

gctfile <- "data/leukemia2.gct"
clsfile <- "data/leukemia2.cls"
chipfile <- "data/leukemia2.chip"
icsy <- icossy(gctfile = gctfile,
               chipfile = chipfile, 
               clsfile = clsfile, 
               network = "pathwayapi", 
               nmis = 15, 
               frank = T, 
               qnorm = F, 
               ztrans = F,
               sig.test = "ttest")
print(icsy$status)