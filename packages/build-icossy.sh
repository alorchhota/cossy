#############################################
##### commands to create icossy package #####
#############################################

cp ./../icossy/icossy.R ./icossy/R/icossy.R
R CMD build icossy
R CMD check icossy