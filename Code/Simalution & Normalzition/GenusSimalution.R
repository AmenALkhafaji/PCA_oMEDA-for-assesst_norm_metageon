# Authors: Amen Al Khafaji, Carolina Gómez-Llorente, José Camacho
# Contact: amen.a.khabeer@uotechnology.edu.iq
# Date: 12/4/2024

# Clear environment
rm(list=ls())

# Load necessary libraries
library(tidyr)
library(dplyr)

# Simulate a healthy gut microbiome ecosystem with genus-level abundances (scaled to 1,000,000)
healthy_ecosystem <- c(rep("A",0.15*1000000), rep("B",0.6*1000000), rep("C",1.5*1000000), rep("D",1.6*1000000),
                       rep("E",12.2*1000000), rep("FF",4.1*1000000), rep("Gg",30.8*1000000), rep("H",2.3*1000000),
                       rep("I",0.14*1000000), rep("J",1.21*1000000), rep("K",0.18*1000000), rep("L",0.32*1000000),
                       rep("M",3.31*1000000), rep("N",0.63*1000000), rep("O",0.17*1000000), rep("P",0.34*1000000),
                       rep("Q",1.08*1000000), rep("R",0.14*1000000), rep("S",2.15*1000000), rep("TT",2.37*1000000),
                       rep("U",0.033*1000000), rep("V",0.47*1000000), rep("W",0.34*1000000), rep("X",5.46*1000000),
                       rep("Y",3.47*1000000), rep("Z",6.42*1000000), rep("pp",4.25*1000000), rep("ppA",2.81*1000000),
                       rep("ppB",0.36*1000000), rep("ppC",1.02*1000000), rep("ppD",0.83*1000000), rep("ppE",1.80*1000000),
                       rep("ppFF",0.31*1000000), rep("pptr",1), rep("ppH",0.17*1000000), rep("ppI",7.1*1000000))

# Simulate an unhealthy gut microbiome ecosystem with altered abundances
unhealthy_ecosystem <- c(rep("A",0.38*1000000), rep("B",0.6*1000000), rep("C",1.5*1000000), rep("D",1),
                         rep("E",26*1000000), rep("FF",10), rep("Gg",20), rep("H",11*1000000),
                         rep("I",0.7*1000000), rep("J",1*1000000), rep("K",0.14*1000000), rep("L",1*1000000),
                         rep("M",0.15*1000000), rep("N",0.3*1000000), rep("O",1), rep("P",0.3*1000000),
                         rep("Q",3*1000000), rep("R",0.14*1000000), rep("S",20), rep("TT",0.5*1000000),
                         rep("U",1*1000000), rep("V",0.14*1000000), rep("W",2*1000000), rep("X",4.8*1000000),
                         rep("Y",1), rep("Z",0.2*1000000), rep("pp",50), rep("ppA",5*1000000), rep("ppB",3),
                         rep("ppC",2.7*1000000), rep("ppD",0.8*1000000), rep("ppE",8*1000000),
                         rep("ppFF",0.22*1000000), rep("pptr",1), rep("ppH",0.31*1000000), rep("ppI",1.23*1000000))

# Convert character vectors to factors
healthy_ecosystem <- type.convert(healthy_ecosystem, as.is = FALSE)
unhealthy_ecosystem <- type.convert(unhealthy_ecosystem, as.is = FALSE)

# Initialize summary structure for storing genus-level frequencies
sample_summary <- list(
  A=0, B=0, C=0, D=0, E=0, FF=0, G=0, H=0, I=0, Gg=0, K=0, L=0, M=0, N=0, O=0, P=0, Q=0, R=0, S=0, TT=0,
  U=0, V=0, W=0, X=0, Y=0, Z=0, pp=0, ppA=0, ppB=0, ppC=0, ppD=0, ppE=0, ppFF=0, pptr=0, ppH=0, ppI=0
)
sample_summary <- data.frame(sample_summary, stringsAsFactors = TRUE)
num_samples <- 151
sample_summary <- as.list(sample_summary)

# Generate samples for the healthy ecosystem
n <- 1
while (n < num_samples) {
  sample_data <- sample(x = healthy_ecosystem, size = sample(x = 100:300, size = 1), replace = FALSE)
  freq <- summary(sample_data)
  for (i in seq_along(freq)) {
    sample_summary[[names(freq)[i]]] <- rbind(sample_summary[[names(freq)[i]]], freq[i])
  }
  n <- n + 1
}

# Generate samples for the unhealthy ecosystem
n <- 1
while (n < num_samples) {
  sample_data <- sample(x = unhealthy_ecosystem, size = sample(x = 100:300, size = 1), replace = FALSE)
  freq <- summary(sample_data)
  for (i in seq_along(freq)) {
    sample_summary[[names(freq)[i]]] <- rbind(sample_summary[[names(freq)[i]]], freq[i])
  }
  n <- n + 1
}

# Combine all genus frequency data into a single data frame
final_df <- do.call(data.frame, sample_summary)
final_df <- final_df[-1,]  # Remove initial row used for initialization
final_df$Condition <- factor(c(rep("Healthy", num_samples - 1), rep("Unhealthy", num_samples - 1)))

# Rename genus columns to more descriptive names
colnames(final_df) <- c("Bifidobacterium","Butyricimonas","Odoribacter","Paraprevotella",
                        "Bacteroides","Parabacteroides","Prevotella","Unknown1","Unknown2",
                        "Unknown3","Clostridium","Unknown4","Unknown5","Ruminococcus","Anaerostipes",
                        "Blautia","Coprococcus","Dorea","Lachnospira","Roseburia","Unknown6",
                        "Unknown7","Unknown8","Unknown9","Faecalibacterium","Oscillospira",
                        "Ruminococcus_dup","Dialister","Unknown10","Unknown11","Unknown12","Sutterella",
                        "Escherichia","Klebsiella","Haemophilus","Akkermansia","Condition")

# Save the dataset to CSV file
write.csv2(final_df, "Genus.csv", row.names = FALSE)
