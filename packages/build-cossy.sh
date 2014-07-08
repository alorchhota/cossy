############################################
##### commands to create cossy package #####
############################################

cp ./../cossy.R ./cossy/R/cossy.R
R CMD build cossy
R CMD check cossy