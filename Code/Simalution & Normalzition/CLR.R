


#########################################################################
#              Code:  CLR Code                                  #               
#              Author: Amen Adnan Khabeer; Jose Comacho; Crolina Comez  #                              #      
#              Date: 19/02/2024                                         #
#              Email: amen.a.khabeer@uotechnology.edu.iq                #
#                   ;josecamacho@ugr.es;gomezll@ugr.es                  #  
#                                                                       #
#                                                                       #
#                                                                       #
#   Note : please load your Simulation or data                          #
#########################################################################



rm(list=ls())
install.packages("easyCODA")
library(easyCODA)



rm(list=ls())

# Load your data Simalution after Run your Simulation( pylum or Genus)
datos <- read.csv(file="Phylum.csv", header=TRUE, sep=";")

# Remove Tag from the dataset sample
tags <- datos$tags

datos$tags <- NULL
datos[datos==0] <- 0.001
#aplicamos centerd log-ratio transformation
clrlist <- CLR(datos, weight = F)

#a?adimos de nuevo el vector "y"
df <- data.frame(clrlist$LR)
df$tags <- tags 

# Save the CLR values to a CSV file
write.csv2(df,"CLR.csv", row.names = FALSE)

