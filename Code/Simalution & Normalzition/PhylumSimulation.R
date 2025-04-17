# Authors: Amen Al Khafaji, Carolina Gómez-Llorente, José Camacho
# Contact: amen.a.khabeer@uotechnology.edu.iq
# Date: 12/4/2024
# Description: Simulation and sampling of gut microbiota profiles for healthy and unhealthy ecosystems
#              with the purpose of classification and comparative analysis.

# Clean the environment
rm(list = ls())

# Load required libraries
library(tidyr)
library(dplyr)

# ---------------------------------------------
# Define healthy microbiota composition
# ---------------------------------------------
healthy_microbiome <- c(rep("Bacillota", 50 * 1000000),
                        rep("Bacteroidetes", 35 * 1000000),
                        rep("Actinobacteria", 0.5 * 1000000),
                        rep("Verrucomicrobiota", 6 * 1000000),
                        rep("Proteobacteria", 2.2 * 1000000))

# ---------------------------------------------
# Define unhealthy microbiota composition
# ---------------------------------------------
unhealthy_microbiome <- c(rep("Bacillota", 45 * 1000000),
                          rep("Bacteroidetes", 40 * 1000000),
                          rep("Actinobacteria", 1.15 * 1000000),
                          rep("Verrucomicrobiota", 3 * 1000000),
                          rep("Proteobacteria", 10 * 1000000))

# Convert microbiome vectors into factors
healthy_microbiome <- type.convert(healthy_microbiome, as.is = FALSE)
unhealthy_microbiome <- type.convert(unhealthy_microbiome, as.is = FALSE)

# Initialize dataframe for storing sample summaries
sample_summary <- NULL
sample_summary$Bacillota <- 0
sample_summary$Bacteroidetes <- 0
sample_summary$Actinobacteria <- 0
sample_summary$Verrucomicrobiota <- 0
sample_summary$Proteobacteria <- 0
sample_summary <- data.frame(sample_summary, stringsAsFactors = TRUE)

# Number of samples to generate per group
num_samples <- 151

# ---------------------------------------------
# Generate samples from healthy microbiome
# ---------------------------------------------
i <- 1
sample_summary <- as.list(sample_summary)

while (i < num_samples) {
  sample_healthy <- sample(x = healthy_microbiome, size = sample(x = 100:300, size = 1), replace = FALSE)
  counts <- summary(sample_healthy)
  sample_summary$Bacillota <- rbind(sample_summary$Bacillota, counts[1])
  sample_summary$Bacteroidetes <- rbind(sample_summary$Bacteroidetes, counts[2])
  sample_summary$Actinobacteria <- rbind(sample_summary$Actinobacteria, counts[3])
  sample_summary$Verrucomicrobiota <- rbind(sample_summary$Verrucomicrobiota, counts[4])
  sample_summary$Proteobacteria <- rbind(sample_summary$Proteobacteria, counts[5])
  i <- i + 1
}

# ---------------------------------------------
# Generate samples from unhealthy microbiome
# ---------------------------------------------
i <- 1
while (i < num_samples) {
  sample_unhealthy <- sample(x = unhealthy_microbiome, size = sample(x = 100:300, size = 1), replace = FALSE)
  counts <- summary(sample_unhealthy)
  sample_summary$Bacillota <- rbind(sample_summary$Bacillota, counts[1])
  sample_summary$Bacteroidetes <- rbind(sample_summary$Bacteroidetes, counts[2])
  sample_summary$Actinobacteria <- rbind(sample_summary$Actinobacteria, counts[3])
  sample_summary$Verrucomicrobiota <- rbind(sample_summary$Verrucomicrobiota, counts[4])
  sample_summary$Proteobacteria <- rbind(sample_summary$Proteobacteria, counts[5])
  i <- i + 1
}

# Combine all collected samples into one dataframe
sample_summary <- data.frame(sample_summary$Bacillota,
                             sample_summary$Bacteroidetes,
                             sample_summary$Actinobacteria,
                             sample_summary$Verrucomicrobiota,
                             sample_summary$Proteobacteria)

# Remove initial zero row
sample_summary <- sample_summary[-1, ]

# Add group labels: 1 for healthy, 2 for unhealthy
sample_summary$group <- cbind(c(rep(1, num_samples - 1), rep(2, num_samples - 1)))
sample_summary$group <- factor(sample_summary$group)

# ---------------------------------------------
# Calculate percentage of zeros in the dataset
# ---------------------------------------------
zero_counts <- nrow(filter(sample_summary, Bacillota == "0")) +
               nrow(filter(sample_summary, Bacteroidetes == "0")) +
               nrow(filter(sample_summary, Actinobacteria == "0")) +
               nrow(filter(sample_summary, Verrucomicrobiota == "0")) +
               nrow(filter(sample_summary, Proteobacteria == "0"))

zero_percentage <- zero_counts / (5 * (num_samples - 1) * 2) * 100

# Rename the columns with accurate phylum names
colnames(sample_summary) <- c("Bacillota", "Bacteroidetes", "Actinobacteria", "Verrucomicrobiota", "Proteobacteria", "group")

# Export the data to a CSV file
write.csv2(sample_summary, "Phylum.csv", row.names = FALSE)
