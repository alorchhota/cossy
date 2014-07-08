#################################################################
##### This script prepares data for icossy package  #############
##### The rda file will be binded with the package  #############
#################################################################

library(jsonlite)

dataOutputFile <- "results/icossy.rda"

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


#save(pathwayapi, pathwayapi_json, file = dataOutputFile)
save(pathwayapi, file = "./results/pathwayapi.rda")
save(pathwayapi_json, file = "./results/pathwayapi_json.rda")

cat('Data saved in the directory ', getwd(), '/results/', sep = "")
