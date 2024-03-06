



#########################################################################
#              Code: RPKM  Code                                  #               
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

# Calculate RPKM and length of each gene in base pairs

gene_length <- c(1000, 2000, 1500)

# total number of reads in the experiment
total_reads <- sum(data) 

# Calculate RPKM
rpkm <- as.data.frame(t(t(data) / (gene_length/1000) / (total_reads/1000000)))

# Save the rpkm values to a CSV file
rpkm$tags<-tags
write.csv(rpkm,"F:\\rpkm.csv", row.names = FALSE)
