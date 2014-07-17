#################################################################
##### This script prepares data for icossy package  #############
##### The rda file will be binded with the package  #############
#################################################################

library(igraph)
library(jsonlite)


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

kegg <- readGmtFile("./data/kegg.gmt")
kegg_json <- readJsonGmtFile("./data/kegg.json")

string <- readGmtFile("./data/string.gmt")
string_json <- readJsonGmtFile("./data/string.json")


save(pathwayapi, file = "./results/pathwayapi.rda", compress="bzip2")
save(pathwayapi_json, file = "./results/pathwayapi_json.rda", compress="bzip2")
save(kegg, file = "./results/kegg.rda", compress="bzip2")
save(kegg_json, file = "./results/kegg_json.rda", compress="bzip2")
save(string, file = "./results/string.rda", compress="bzip2")
save(string_json, file = "./results/string_json.rda", compress="bzip2")

cat('Data saved in the directory ', getwd(), '/results/', sep = "")
