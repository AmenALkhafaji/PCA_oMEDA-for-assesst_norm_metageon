# Install required packages (only needs to be done once)
install.packages('Rarefy')
install.packages('phyloregion')

# Load required packages
require(Rarefy)
require(ade4)
require(adiv)
require(ape)
require(vegan)
require(phyloregion)
require(raster)

# Clear the environment
rm(list = ls())

# Set the working directory to where your data file is located
setwd("D:\\")

# Import the CSV file containing species or taxonomic abundance data
raw_data <- read.csv(file = "Phylum.csv", header = FALSE, sep = ";")

# Convert data to a proper data frame
raw_data <- as.data.frame(raw_data)

# Convert all columns to numeric (in case they're read as character)
numeric_data <- as.data.frame(apply(raw_data, 2, as.numeric))

# Extract the grouping or tag column (e.g., sample site, ID, etc.)
sample_tags <- numeric_data[, 6]

# Remove the tag column from the main abundance data
abundance_data <- numeric_data[, -6]

# Find the minimum total count across all rows (samples)
minimum_sample_size <- min(rowSums(abundance_data))

# Perform rarefaction to standardize sampling effort across samples
rarefied_data <- as.data.frame(rrarefy(abundance_data, minimum_sample_size))

# Export the rarefied data to a CSV file
write.csv(rarefied_data, "rarefaction.csv", row.names = TRUE)
