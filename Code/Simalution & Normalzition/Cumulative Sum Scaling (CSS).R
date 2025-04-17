library(metagenomeSeq)


#########################################################################
#              Code: Cumialtative Sum Scaling Code                                         #               
#              Author: Amen Adnan Khabeer; Jose Comacho; Crolina Comez  #                              #      
#              Date: 19/02/2024                                         #
#              Email: amen.a.khabeer@uotechnology.edu.iq                #
#                   ;josecamacho@ugr.es;gomezll@ugr.es                  #  
#                                                                       #
#                                                                       #
#                                                                       #
#   Note : please load your Simulation or data                          #
#########################################################################


# Load your data Simulation after Run your Simulation( pylum or Genus); Take in your account that your data saved in drive D


setwd("D:\\")
Rrare <- read.csv(file="Phylum.csv", header=TRUE, sep=";")


# Remove Tag from the dataset sample

tags <- Rrare$tags

Rrare$tags <- NULL

# set replacment of zero reads in samples 

Rrare[Rrare==0] <- 0.001

metaSeqObject      = newMRexperiment(Rrare) 

# find scale factor of CSS

metaSeqObject_CSS  = cumNorm( metaSeqObject, p=cumNormStatFast(metaSeqObject) )

# extract norlamzed data from metaSeqObject_CSS

OTU_read_count_CSS = data.frame(MRcounts(metaSeqObject_CSS, norm=TRUE, log=TRUE))

#retrun Tag to dataset
OTU_read_count_CSS$tags<-tags

# Save the RLE values to a CSV file

write.csv(OTU_read_count_CSS, "CSS.csv")

