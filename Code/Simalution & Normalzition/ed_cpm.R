BiocManager::install("edgeR")
library(edgeR)
setwd("D:\\")
data <- read.csv(file="D:\\phylum1.csv",header = FALSE,sep = ";")

n=cpm(data, normalized.lib.sizes = TRUE,log = FALSE, prior.count = 2)
write.csv(n,"D:\\egdr.csv", row.names = FALSE)
