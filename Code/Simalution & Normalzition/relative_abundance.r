
#########################################################################
#              Code: Relative Abundance Code                                  #               
#              Author: Amen Adnan Khabeer; Jose Comacho; Crolina Comez  #                              #      
#              Date: 19/02/2024                                         #
#              Email: amen.a.khabeer@uotechnology.edu.iq                #
#                   ;josecamacho@ugr.es;gomezll@ugr.es                  #  
#                                                                       #
#                                                                       #
#                                                                       #
#   Note : please load your Simulation or data                          #
#########################################################################

# Load your data Simalution after Run your Simulation( pylum or Genus); Take in your account that your data saved in drive D


setwd("D:\\")

data <- read.csv(file="genus.csv", header=TRUE, sep=";")

# Remove Tag from the dataset sample

tags <- data$tags
data$tags <- NULL

# calculate the total read counts for each sample

total_counts <- apply(data, 1, sum)

# calculate the relative abundance of each taxon in each sample
relative_abundance <- as.data.frame(t(t(data) / total_counts))

relative_abundance$tags<-tags

write.csv(relative_abundance, "D:\\relative_abundance.csv")

