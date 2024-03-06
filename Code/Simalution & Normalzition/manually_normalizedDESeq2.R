


#########################################################################
#              Code: DESeq Code                                         #               
#              Author: Amen Adnan Khabeer; Jose Comacho; Crolina Comez  #                              #      
#              Date: 19/02/2024                                         #
#              Email: amen.a.khabeer@uotechnology.edu.iq                #
#                   ;josecamacho@ugr.es;gomezll@ugr.es                  #  
#                                                                       #
#                                                                       #
#                                                                       #
#   Note : please load your Simulation or data                          #
#########################################################################



# Load your data Simulation after Run your Simulation( pylum or Genus); Take in your account that your data saved in drive D


setwd("D:\\")
data <- read.csv(file="genus.csv", header=TRUE, sep=";")


# Remove Tag from the dataset sample

tags <- data$tags

data$tags <- NULL

data[data == 0] <- 0.001
# compute DESeq2 
log_data = log(data)
log_data = log_data %>% 
rownames_to_column('gene') %>% 
mutate (pseudo_reference = rowMeans(log_data))
filtered_log_data = log_data %>% filter(pseudo_reference != "-Inf")

ratio_data = sweep(filtered_log_data[,2:28], 1, filtered_log_data[,29], "-")
sample_medians = apply(ratio_data, 2, median)
# Find Sacle_Factor 
scaling_factors = exp(sample_medians)
# Apply Sacle factor 
manually_normalized = sweep(data, 2, scaling_factors, "/")

as.data.frame(manually_normalized$tags)<-tags


# Save normaliezed data by DESeq2
write.csv(manually_normalized, "F:\\DESeq.csv")

