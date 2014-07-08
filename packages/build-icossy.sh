#############################################
##### commands to create icossy package #####
#############################################

cp ./../icossy/icossy.R ./icossy/R/icossy.R
cp ./../data/pathwayapi.rda ./icossy/data/pathwayapi.rda
cp ./../data/pathwayapi_json.rda ./icossy/data/pathwayapi_json.rda

R CMD build icossy
R CMD check icossy