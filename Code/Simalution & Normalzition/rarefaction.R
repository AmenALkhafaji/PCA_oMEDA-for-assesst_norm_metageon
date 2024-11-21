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
rm(list=ls())
#Importamos los datos
setwd("D:\\")
datos <- read.csv(file="D:\\Book15021.csv",header = FALSE,sep = ";")
datos <-as.data.frame(datos )


datos <- as.data.frame(apply(datos, 2, as.numeric))
tag<-datos[,37]

if (tag[1]==1) 
{
  datos<-datos[,-37] 

(raremax <- min(rowSums(datos)))

Rraree <- as.data.frame(rrarefy(datos, raremax))
i=1;

Rraree1=NULL
n=51
while ( i<n)
{
  Rraree <- as.data.frame(rrarefy(Rraree, raremax))
  Rraree$itreation <-i
  Rraree$tag <-1
  Rraree1 <- as.data.frame(((rbind(Rraree1,Rraree))))
  i<-i+1 
 }


}


datos2 <- read.csv(file="D:\\Book15022.csv",header = FALSE,sep = ";")
#remove tags from the dataframe to proccess the recored only

#datos2 <-datos2[-1, ] #Establish the rarefaction level: = DL.min that defined as lowest depth of samples.

datos2 <- as.data.frame(apply(datos2, 2, as.numeric))
tag2<-datos2[,37]
if (tag2[1]==2) 
{
  datos2<-datos2[,-37] 
  
  (raremax2 <- min(rowSums(datos)))
  
  Rraree2 <- as.data.frame(rrarefy(datos2, raremax2))
  i2=1;
  
  Rraree4=NULL
  
  while ( i2<n)
  {
    Rraree2 <- as.data.frame(rrarefy(Rraree2, raremax2))
    Rraree2$itreation <-i2
    Rraree2$tag <-2
    Rraree4 <- as.data.frame(((rbind(Rraree4,Rraree2))))
    i2<-i2+1 
  }
  
  
}


Rraree4 <- as.data.frame(((rbind(Rraree1,Rraree4))))
write.csv(Rraree4, "d:\\final_G_50_2.csv", row.names = TRUE)

