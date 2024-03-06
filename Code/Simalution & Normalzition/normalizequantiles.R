




#########################################################################
#              Code: Quantile Normalization Code                                  #               
#              Author: Amen Adnan Khabeer; Jose Comacho; Crolina Comez  #                              #      
#              Date: 19/02/2024                                         #
#              Email: amen.a.khabeer@uotechnology.edu.iq                #
#                   ;josecamacho@ugr.es;gomezll@ugr.es                  #  
#                                                                       #
#                                                                       #
#                                                                       #
#   Note : please load your Simulation or data                          #
#########################################################################

install.packages("preprocessCore")
library(preprocessCore)

###########################################################################


# Load your data Simalution after Run your Simulation( pylum or Genus); Take in your account that your data saved in drive D

setwd("D:\\")

data<- read.csv(file="genus.csv", header=TRUE, sep=";")



# Remove Tag from the dataset sample
tags <- data$tags
data$tags <- NULL

# Perform quantile normalization

data_matrix <- as.matrix(data)

# Perform quantile normalization

normalized_data <- normalize.quantiles(data_matrix)

# Convert back to dataframe

normalized_data_df <- as.data.frame(normalized_data)

# Write output CSV file
normalized_data_df$tags<-tags
write.csv(normalized_data_df,"D:\\QuantileNormalization.csv", row.names = FALSE)

