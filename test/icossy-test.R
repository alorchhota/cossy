source("./icossy/icossy.R")
source("./icossy/data-preparation.R")



gctfile <- paste0("data/leukemia.gct")
clsfile <- paste0("data/leukemia.cls")

debug(icossy)
icsy <- icossy(gctfile = gctfile,
               chipfile = NA, 
               clsfile = clsfile, 
               network = "pathwayapi", 
               nmis = 15, 
               frank = T, 
               qnorm = F, 
               ztrans = F,
               sig.test = "ttest")