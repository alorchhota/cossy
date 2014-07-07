#source("./icossy/icossy.R")

readGmtFile <- function(gmt_file_path){
  gmt_file <- file(gmt_file_path, "r")
  lines <- readLines(gmt_file)
  close(gmt_file)
  splittedLines <- strsplit(lines, split='\t')
  lapply(splittedLines, function(set) list(id=set[1], name=set[2], genes=set[-c(1,2)]));
}

readJsonGmtFile <- function(json_gmt_file_path){
  txt <- readChar(json_gmt_file_path, file.info(json_gmt_file_path)$size)
  jsonGmt <- fromJSON(txt, simplifyVector = F)
  return(jsonGmt)
}

pathwayapi <- readGmtFile("./data/pathwayapi.gmt")
pathwayapi_json <- readJsonGmtFile("./data/pathwayapi.json")


gctfile <- paste0("data/cns.gct")
clsfile <- paste0("data/cns.cls")


icsy <- icossy(gctFile = gctFile, chipfile = NA, clsFile = clsfile, network = "pathwayapi", nmis = 5, frank = T, qnorm = F, ztrans = F)
