source("./icossy/icossy.R")
source("./icossy/data-preparation.R")



gctfile <- paste0("data/cns.gct")
clsfile <- paste0("data/cns.cls")


icsy <- icossy(gctFile = gctFile,
               chipfile = NA, 
               clsFile = clsfile, 
               network = "pathwayapi", 
               nmis = 5, 
               frank = T, 
               qnorm = F, 
               ztrans = F,
               sig.test = "ttest")
