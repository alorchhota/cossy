#############################################
##### commands to create icossy package #####
#############################################

## create data folder
cd icossy
mkdir data
cd ..

## move code files
cp ./../icossy/icossy.R ./icossy/R/icossy.R

## move data files
cp ./../data/pathwayapi.rda ./icossy/data/pathwayapi.rda
cp ./../data/pathwayapi_json.rda ./icossy/data/pathwayapi_json.rda
cp ./../data/kegg.rda ./icossy/data/kegg.rda
cp ./../data/kegg_json.rda ./icossy/data/kegg_json.rda
cp ./../data/string.rda ./icossy/data/string.rda
cp ./../data/string_json.rda ./icossy/data/string_json.rda

## build package
R CMD build icossy --resave-data
R CMD check icossy