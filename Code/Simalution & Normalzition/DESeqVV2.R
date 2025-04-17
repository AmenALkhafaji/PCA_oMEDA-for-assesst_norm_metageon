#########################################################################
#              Code: DESeq Manual Normalization                         #
#              Author: Amen Adnan Khabeer; Jose Comacho; Carolina Gomez #
#              Date: 19/02/2024                                          #
#              Emails: amen.a.khabeer@uotechnology.edu.iq               #
#                      josecamacho@ugr.es; gomezll@ugr.es               #
#                                                                       #
#   Note: Load your simulation data (Phylum or Genus level)             #
#         Ensure your data file is saved in drive D:\                   #
#########################################################################

# Load required libraries
library(dplyr)
library(tibble)

# Load dataset (Phylum or Genus level abundance)
raw_data <- read.csv(file = "Phylum.csv", header = TRUE, sep = ";")

# Extract sample tags (e.g., group labels)
sample_tags <- raw_data$tags

# Remove 'tags' column to keep only count data
abundance_matrix <- raw_data[, -6]

# Replace zero counts with 1 to avoid log(0) issues
abundance_matrix[abundance_matrix == 0] <- 1

# Log-transform the data
log_transformed_data <- log(abundance_matrix)

# Add rownames as a column to prepare for normalization
log_transformed_data <- log_transformed_data %>%
  rownames_to_column(var = "feature") %>%
  mutate(pseudo_reference = rowMeans(select(., -feature)))

# Remove any rows where pseudo_reference is -Inf
filtered_log_data <- log_transformed_data %>%
  filter(pseudo_reference != -Inf)

# Compute log ratios to pseudo-reference
log_ratios <- sweep(filtered_log_data[, 2:(ncol(filtered_log_data) - 1)], 1, filtered_log_data$pseudo_reference, "-")

# Compute median of log ratios per sample
sample_medians <- apply(log_ratios, 2, median)

# Calculate scaling factors using the exponent of sample medians
scaling_factors <- exp(sample_medians)

# Apply scaling factors to normalize original (non-log) data
deseq_normalized_data <- sweep(abundance_matrix, 2, scaling_factors, "/")

# Reattach the sample tags
deseq_normalized_data <- as.data.frame(deseq_normalized_data)
deseq_normalized_data$tags <- sample_tags

# Save DESeq-normalized data to CSV
write.csv(deseq_normalized_data, "DESeq.csv", row.names = FALSE)
