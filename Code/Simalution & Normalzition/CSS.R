#########################################################################
#              Code: Cumulative Sum Scaling (CSS) Normalization         #
#              Author: Amen Adnan Khabeer; Jose Comacho; Carolina Gomez #
#              Date: 19/02/2024                                          #
#              Emails: amen.a.khabeer@uotechnology.edu.iq               #
#                      josecamacho@ugr.es; gomezll@ugr.es               #
#                                                                       #
#   Note: Load your simulation results (e.g., Phylum or Genus data)     #
#         before running this script. Ensure the file is saved on D:\   #
#########################################################################

# Load required package
library(metagenomeSeq)

# Set working directory where your dataset is saved
setwd("D:\\")

# Load your dataset (e.g., Phylum abundance table with metadata)
abundance_data <- read.csv(file = "Phylum.csv", header = TRUE, sep = ";")

# Extract sample labels (e.g., groups, tags, or IDs)
sample_labels <- abundance_data$tags

# Remove the 'tags' column to retain only abundance values
abundance_matrix <- abundance_data[, -6]

# Replace zero values to avoid issues with log transformations
abundance_matrix[abundance_matrix == 0] <- 0.001

# Create a metagenomeSeq MRexperiment object with the count data
meta_seq_object <- newMRexperiment(abundance_matrix)

# Calculate normalization scaling factors and apply CSS normalization
meta_seq_object_css <- cumNorm(meta_seq_object, p = cumNormStatFast(meta_seq_object))

# Extract normalized, log-transformed abundance values
normalized_log_counts <- data.frame(MRcounts(meta_seq_object_css, norm = TRUE, log = TRUE))

# Add the sample labels back to the normalized dataset
normalized_log_counts$tags <- sample_labels

# Save the final normalized dataset to a CSV file
write.csv(normalized_log_counts, "CSS.csv")
