#########################################################################
#              Code:  CLR Transformation                               #
#              Author: Amen Adnan Khabeer; Jose Comacho; Carolina Gomez#
#              Date: 12/04/2024                                        #
#              Emails: amen.a.khabeer@uotechnology.edu.iq              #
#                      josecamacho@ugr.es; gomezll@ugr.es              #
#                                                                       #
#   Note : Please load your simulation or data before running          #
#########################################################################

# Clear the current R environment
rm(list = ls())

# Install and load the required package for compositional data analysis
install.packages("easyCODA")
library(easyCODA)

# Clear again just in case (not strictly necessary after the first `rm`)
rm(list = ls())

# Load your dataset (e.g., simulated Phylum or Genus level abundance data)
abundance_data <- read.csv(file = "Phylum.csv", header = TRUE, sep = ";")

# Extract the sample labels or group identifiers
sample_labels <- abundance_data$tags

# Remove the 'tags' column from the abundance matrix
abundance_matrix <- abundance_data[, -6]

# Replace all zero values with a small number to avoid log(0) errors
abundance_matrix[abundance_matrix == 0] <- 1

# Apply Centered Log-Ratio (CLR) transformation to the abundance data
clr_result <- CLR(abundance_matrix, weight = FALSE)

# Create a new data frame containing CLR-transformed values
clr_dataframe <- data.frame(clr_result$LR)

# Add the sample labels back to the CLR-transformed data
clr_dataframe$tags <- sample_labels

# Save the final CLR-transformed data to a CSV file
write.csv2(clr_dataframe, "CLR.csv", row.names = FALSE)
