source("./icossy/icossy.R")
source("./icossy/data-preparation.R")



gctfile <- paste0("data/cns.gct")
clsfile <- paste0("data/cns.cls")


icsy <- icossy(gctfile = gctfile,
               chipfile = NA, 
               clsfile = clsfile, 
               network = "pathwayapi", 
               nmis = 5, 
               frank = T, 
               qnorm = F, 
               ztrans = F,
               sig.test = "ttest")
