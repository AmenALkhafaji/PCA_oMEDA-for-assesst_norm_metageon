rm(list=ls())

install.packages('Rarefy')
install.packages('phyloregion')


#Required packages
require(Rarefy)
require(ade4)
require(adiv)
require(ape)
require(vegan)
require(phyloregion)
require(raster)

###########################################################################



#Importamos los datos
setwd("D:\\")
datos <- read.csv(file="genus.csv", header=TRUE, sep=";")

#remove tags from the dataframe to proccess the recored only
tags <- datos$tags
datos$tags <- NULL

#Establish the rarefaction level: = DL.min that defined as lowest depth of samples.

(raremax <- min(rowSums(datos)))


Rrare <- rrarefy(datos, raremax)


datos=data.frame(Rrare)
datos$tags <-tags


write.csv2(datos,"d:\\rrare.csv", row.names = FALSE)


