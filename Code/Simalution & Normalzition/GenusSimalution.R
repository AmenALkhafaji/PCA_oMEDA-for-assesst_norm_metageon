# Authors: Amen Al Khafaji, Carolina Gómez-Llorente, José Camacho
# Contact: amen.a.khabeer@uotechnology.edu.iq
# Date: 12/4/2024

rm(list=ls())
library(tidyr)
library(dplyr)

# Healthy Ecosystem
healthy_ecosystem <- c( c(rep("A",0.15 *1000000)) ,c(rep("B",0.6*1000000)),
                        c(rep("C",1.5*1000000)), c(rep("D",1.6*1000000)),
                        c(rep("E",12.2*1000000)), c(rep("FF",4.1*1000000)),
                        c(rep("Gg",30.8*1000000)), c(rep("H",2.3*1000000)),
                        c(rep("I",0.14*1000000)), c(rep("J",1.21*1000000)),
                        c(rep("K",0.18*1000000)), c(rep("L",0.32*1000000)),
                        c(rep("M",3.31*1000000)), c(rep("N",0.63*1000000)),
                        c(rep("O",0.17*1000000)), c(rep("P",0.34*1000000)),
                        c(rep("Q",1.08*1000000)), c(rep("R",0.14*1000000)),
                        c(rep("S",2.15*1000000)), c(rep("TT",2.37*1000000)),
                        c(rep("U",0.033*1000000)), c(rep("V",0.47*1000000)),
                        c(rep("W",0.34*1000000)), c(rep("X",5.46*1000000)),
                        c(rep("Y",3.47*1000000)), c(rep("Z",6.42*1000000)),
                        c(rep("pp",4.25*1000000)),
                        c(rep("ppA",2.81*1000000)),c(rep("ppB",0.36*1000000)),
                        c(rep("ppC",1.02*1000000)),c(rep("ppD",0.83*1000000)),
                        c(rep("ppE",1.80*1000000)),
                        c(rep("ppFF",0.31*1000000)),c(rep("pptr",0.001*1000000)),c(rep("ppH",0.17*1000000)),c(rep("ppI",7.1*1000000)))

# Sick Ecosystem
sick_ecosystem <- c( c(rep("A",0.38*1000000)) ,c(rep("B",0.6*1000000)),
                     c(rep("C",1.5*1000000)), c(rep("D",0.001*1000000)),
                     c(rep("E",26*1000000)), c(rep("FF",0.01*1000000)),
                     c(rep("Gg",0.02*1000000)), c(rep("H",11*1000000)),
                     c(rep("I",0.7*1000000)), c(rep("J",1*1000000)),
                     c(rep("K",0.14*1000000)), c(rep("L",1*1000000)),
                     c(rep("M",0.15*1000000)), c(rep("N",0.3*1000000)),
                     c(rep("O",0.001*1000000)), c(rep("P",0.3*1000000)),
                     c(rep("Q",3*1000000)), c(rep("R",0.14*1000000)),
                     c(rep("S",0.02*1000000)), c(rep("TT",0.5*1000000)),
                     c(rep("U",1*1000000)), c(rep("V",0.14*1000000)),
                     c(rep("W",2*1000000)), c(rep("X",4.8*1000000)),
                     c(rep("Y",0.001*1000000)), c(rep("Z",0.2*1000000)),
                     c(rep("pp",0.05*1000000)),
                     c(rep("ppA",5*1000000)),c(rep("ppB",3)),
                     c(rep("ppC",2.7*1000000)),c(rep("ppD",0.8*1000000)),
                     c(rep("ppE",8*1000000)),c(rep("ppFF",0.22*1000000)),c(rep("pptr",0.001*1000000)),
                     c(rep("ppH",0.31*1000000)),c(rep("ppI",1.23*1000000)))

# Convert to factors
healthy_ecosystem <- type.convert(healthy_ecosystem, as.is = FALSE) 
sick_ecosystem <- type.convert(sick_ecosystem, as.is = FALSE) 

# Initialize healthy ecosystem sample data frame
final_result <- data.frame(matrix(0, ncol = 36, nrow = 0))

# Add columns
colnames(final_result) <- c("Bifidobacterium", "Butyricimonas", "Odoribacter", "Paraprevotella",
                            "Bacteroides", "Parabacteroides", "Prevotella", "unknown1", "unknown2",
                            "unknown3", "Clostridium", "unknown4", "unknown5", "Ruminococcus", "Anaerostipes", 
                            "Blautia", "Coprococcus", "Dorea", "Lachnospira", "Roseburia", "unknown6", 
                            "unknown7", "unknown8", "unknown9", "Faecalibacterium", "Oscillospira", 
                            "Ruminococcus", "Dialister", "unknown10", "unknown11", "unknown12", "Sutterella", 
                            "Escherichia", "Klebsiella", "Haemophilus", "Akkermansia", "Tags")

# Healthy ecosystem sample loop
samples_count <- 151
counter <- 1
final_result <- as.list(final_result)

while (counter < samples_count) {
  sample1 <- sample(x = healthy_ecosystem, size = sample(x = 100:300, size = 1), replace = F)
  result <- summary(sample1)
  final_result$Bifidobacterium <- rbind(final_result$Bifidobacterium, result[1])
  final_result$Butyricimonas <- rbind(final_result$Butyricimonas, result[2])
  final_result$Odoribacter <- rbind(final_result$Odoribacter, result[3])
  final_result$Paraprevotella <- rbind(final_result$Paraprevotella, result[4])
  final_result$Bacteroides <- rbind(final_result$Bacteroides, result[5])
  final_result$Parabacteroides <- rbind(final_result$Parabacteroides, result[6])
  final_result$Prevotella <- rbind(final_result$Prevotella, result[7])
  final_result$unknown1 <- rbind(final_result$unknown1, result[8])
  final_result$unknown2 <- rbind(final_result$unknown2, result[9])
  final_result$Clostridium <- rbind(final_result$Clostridium, result[10])
  final_result$unknown4 <- rbind(final_result$unknown4, result[11])
  final_result$unknown5 <- rbind(final_result$unknown5, result[12])
  final_result$Ruminococcus <- rbind(final_result$Ruminococcus, result[13])
  final_result$Anaerostipes <- rbind(final_result$Anaerostipes, result[14])
  final_result$Blautia <- rbind(final_result$Blautia, result[15])
  final_result$Coprococcus <- rbind(final_result$Coprococcus, result[16])
  final_result$Dorea <- rbind(final_result$Dorea, result[17])
  final_result$Lachnospira <- rbind(final_result$Lachnospira, result[18])
  final_result$Roseburia <- rbind(final_result$Roseburia, result[19])
  final_result$unknown6 <- rbind(final_result$unknown6, result[20])
  final_result$unknown7 <- rbind(final_result$unknown7, result[21])
  final_result$unknown8 <- rbind(final_result$unknown8, result[22])
  final_result$unknown9 <- rbind(final_result$unknown9, result[23])
  final_result$Faecalibacterium <- rbind(final_result$Faecalibacterium, result[24])
  final_result$Oscillospira <- rbind(final_result$Oscillospira, result[25])
  final_result$Dialister <- rbind(final_result$Dialister, result[26])
  final_result$unknown10 <- rbind(final_result$unknown10, result[27])
  final_result$unknown11 <- rbind(final_result$unknown11, result[28])
  final_result$unknown12 <- rbind(final_result$unknown12, result[29])
  final_result$Sutterella <- rbind(final_result$Sutterella, result[30])
  final_result$Escherichia <- rbind(final_result$Escherichia, result[31])
  final_result$Klebsiella <- rbind(final_result$Klebsiella, result[32])
  final_result$Haemophilus <- rbind(final_result$Haemophilus, result[33])
  final_result$Akkermansia <- rbind(final_result$Akkermansia, result[34])
  counter <- counter + 1
}

counter <- 1
# Sick ecosystem sample loop
while (counter < samples_count) {
  sample2 <- sample(x = sick_ecosystem, size = sample(x = 100:300, size = 1), replace = F)
  result <- summary(sample2)
  final_result$Bifidobacterium <- rbind(final_result$Bifidobacterium, result[1])
  final_result$Butyricimonas <- rbind(final_result$Butyricimonas, result[2])
  final_result$Odoribacter <- rbind(final_result$Odoribacter, result[3])
  final_result$Paraprevotella <- rbind(final_result$Paraprevotella, result[4])
  final_result$Bacteroides <- rbind(final_result$Bacteroides, result[5])
  final_result$Parabacteroides <- rbind(final_result$Parabacteroides, result[6])
  final_result$Prevotella <- rbind(final_result$Prevotella, result[7])
  final_result$unknown1 <- rbind(final_result$unknown1, result[8])
  final_result$unknown2 <- rbind(final_result$unknown2, result[9])
  final_result$Clostridium <- rbind(final_result$Clostridium, result[10])
  final_result$unknown4 <- rbind(final_result$unknown4, result[11])
  final_result$unknown5 <- rbind(final_result$unknown5, result[12])
  final_result$Ruminococcus <- rbind(final_result$Ruminococcus, result[13])
  final_result$Anaerostipes <- rbind(final_result$Anaerostipes, result[14])
  final_result$Blautia <- rbind(final_result$Blautia, result[15])
  final_result$Coprococcus <- rbind(final_result$Coprococcus, result[16])
  final_result$Dorea <- rbind(final_result$Dorea, result[17])
  final_result$Lachnospira <- rbind(final_result$Lachnospira, result[18])
  final_result$Roseburia <- rbind(final_result$Roseburia, result[19])
  final_result$unknown6 <- rbind(final_result$unknown6, result[20])
  final_result$unknown7 <- rbind(final_result$unknown7, result[21])
  final_result$unknown8 <- rbind(final_result$unknown8, result[22])
  final_result$unknown9 <- rbind(final_result$unknown9, result[23])
  final_result$Faecalibacterium <- rbind(final_result$Faecalibacterium, result[24])
  final_result$Oscillospira <- rbind(final_result$Oscillospira, result[25])
  final_result$Dialister <- rbind(final_result$Dialister, result[26])
  final_result$unknown10 <- rbind(final_result$unknown10, result[27])
  final_result$unknown11 <- rbind(final_result$unknown11, result[28])
  final_result$unknown12 <- rbind(final_result$unknown12, result[29])
  final_result$Sutterella <- rbind(final_result$Sutterella, result[30])
  final_result$Escherichia <- rbind(final_result$Escherichia, result[31])
  final_result$Klebsiella <- rbind(final_result$Klebsiella, result[32])
  final_result$Haemophilus <- rbind(final_result$Haemophilus, result[33])
  final_result$Akkermansia <- rbind(final_result$Akkermansia, result[34])
  counter <- counter + 1
}

final_result <- data.frame(final_result)

# Save to CSV
write.csv2(final_result, "Genus.csv", row.names = FALSE)
