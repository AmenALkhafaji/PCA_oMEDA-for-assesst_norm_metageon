#########################################################################
#              Code: Total Sum Scaling (TSS) Normalization              #
#              Author: Amen Adnan Khabeer; Jose Comacho; Carolina Gomez #
#              Date: 19/02/2024                                          #
#              Emails: amen.a.khabeer@uotechnology.edu.iq               #
#                      josecamacho@ugr.es; gomezll@ugr.es               #
#                                                                       #
#   Note: Load your simulation data (Phylum or Genus level)             #
#         Ensure your data file is saved in drive D:\                   #
#########################################################################

# Load dataset (e.g., Phylum abundance table)
raw_abundance_data <- read.csv(file = "Phylum.csv", header = TRUE, sep = ";")

# Extract sample identifiers or group tags
sample_tags <- raw_abundance_data$tags

# Remove the 'tags' column to retain only numeric abundance values
abundance_matrix <- raw_abundance_data[, -6]

# Apply Total Sum Scaling (TSS) normalization:
# Divide each value by the total sum of all values in the matrix
tss_normalized_data <- abundance_matrix / sum(abundance_matrix)

# Optional: Add sample tags back (if needed)
tss_normalized_data$tags <- sample_tags

# Save the TSS-normalized dataset to a CSV file
write.csv(tss_normalized_data, "TSS.csv", row.names = FALSE)
