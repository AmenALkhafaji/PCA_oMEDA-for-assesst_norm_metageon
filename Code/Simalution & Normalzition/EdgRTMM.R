# Install edgeR package if not already installed
if (!requireNamespace("edgeR", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("edgeR")
}

# Load edgeR library
library(edgeR)

# Read the CSV file, assuming the last column is the group info
data <- read.csv("Phylum.csv", header = TRUE, sep=";")

# Extract the group information from the 37th column
group <- as.factor(data[, 6])  # Assuming 1 = h, 2 = S
print(length(group))
# Remove the last column from the data to get only expression values
expression_data <- data[, 1:5]

# Create a DGEList object (Digital Gene Expression List)
dge <- DGEList(counts = t(expression_data), group = t(group))

# TMM normalization
dge <- calcNormFactors(dge, method = "TMM")

# Extract normalized counts (logCPM if needed)
normalized_counts <- cpm(dge, log = FALSE)

# Save normalized counts to CSV file
write.csv(normalized_counts, file = "EdgRTMM.csv", row.names = TRUE)

