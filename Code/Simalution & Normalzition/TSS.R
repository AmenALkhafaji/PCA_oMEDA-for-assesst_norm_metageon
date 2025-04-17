
#########################################################################
#              Code: Total Sum Scaling Code                                         #               
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



df<- read.csv(file="Phylum.csv", header=TRUE, sep=";")

# Remove Tag from the dataset sample

tags <- df$tags

df$tags <- NULL

# Find  total sum scaling
df_scaled <- df / sum(df)

# Save the TTS values to a CSV file
df$tags<-tags
write.csv(df_scaled, "TSS.csv")

