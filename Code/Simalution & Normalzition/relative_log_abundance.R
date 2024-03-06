



#########################################################################
#              Code: Relative Log Abundance Code                                  #               
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

# Load the gene expression data from the CSV file
gene_expression <- data

# Calculate the total expression levels for each sample
total_expression <- apply(gene_expression, 2, sum)

# Calculate the median expression level across samples
median_expression <- median(total_expression)

# Calculate the Relative Log Expression (RLE) values
rle <- as.data.frame(t(apply(gene_expression, 1, function(x) log2(x/median_expression))))
rle$tags<-tags
# Save the RLE values to a CSV file
write.csv(rle, "D:\\RelativeLog.csv")

