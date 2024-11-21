######Type 2 diabetes######  

rm(list=ls())
library(tidyr)
library(dplyr)

#    ECOSISTEMA SANO

#Firmicutes <- 45.0001%       (A)
#Bacteroidetes <- 45%         (B)
#Actinobacteria <- 3.3333%    (C)
#Verrucomcrobiota <- 3.33333%   (D)
#Proteobacteria <- 3.33333%     (E)
ecosistema1 <- c( c(rep("A",50 *1000)) ,c(rep("B",35*1000)),
                  c(rep("C",0.5*1000)), c(rep("D",6*1000)),
                  c(rep("E",2.2*1000)))

#    ECOSISTEMA ENFERMO

#Firmicutes <- 48%       (A)
#Bacteroidetes <- 22%    (B)
#Actinobacteria <- 12%   (C)
#Proteobacteria <- 12%   (E)
#Fusobacteria <- 6%      (D)

ecosistema2 <- c( c(rep("A",45*1000)) ,c(rep("B",40*1000)),
                  c(rep("C",1.15*1000)), c(rep("D",3*1000)),
                  c(rep("E",10*1000)))

#convertimos en factores

ecosistema1 <- type.convert(ecosistema1, as.is = FALSE) # -> factor
ecosistema2 <- type.convert(ecosistema2, as.is = FALSE) # -> factor

#bucle de muestras sano
fitfinal1 <- NULL
fitfinal1$A <- 0
fitfinal1$B <- 0
fitfinal1$C <- 0
fitfinal1$D <- 0
fitfinal1$E <- 0
fitfinal1 <- data.frame(fitfinal1, stringsAsFactors = T)

nmuestras <- 151   #ej: 21 son 20 muestras 1 y 20 muestras 2 = 40 muestras total 
n <- 1
fitfinal1<- as.list(fitfinal1)
while (n<nmuestras){
  muestra1 <- sample(x=ecosistema1, size=sample(x=100:300,size=1), replace = F)
  fit <- summary(muestra1)
  fitfinal1$A <- rbind(fitfinal1$A,fit[1])
  fitfinal1$B <- rbind(fitfinal1$B,fit[2])
  fitfinal1$C <- rbind(fitfinal1$C,fit[3])
  fitfinal1$D <- rbind(fitfinal1$D,fit[4])
  fitfinal1$E <- rbind(fitfinal1$E,fit[5])
  n <- n+1
}

n <- 1
#bucle de muestras enfermo
while (n<(nmuestras)){
  muestra2 <- sample(x=ecosistema2, size=sample(x=100:300,size=1), replace = F)
  fit <- summary(muestra2)
  fitfinal1$A <- rbind(fitfinal1$A,fit[1])
  fitfinal1$B <- rbind(fitfinal1$B,fit[2])
  fitfinal1$C <- rbind(fitfinal1$C,fit[3])
  fitfinal1$D <- rbind(fitfinal1$D,fit[4])
  fitfinal1$E <- rbind(fitfinal1$E,fit[5])
  n <- n+1
}
fitfinal1 <- data.frame(fitfinal1$A,fitfinal1$B,fitfinal1$C,fitfinal1$D,
                        fitfinal1$E)
fitfinal1 <- fitfinal1[-1,]

fitfinal1$tags <- cbind(c(rep(1, nmuestras-1), rep(2, nmuestras-1)))
fitfinal1[,"tags"]<-factor(fitfinal1[,"tags"]) #tags to FACTORS


#PORCENTAJE DE CEROS DE NUESTRO DATAFRAME
ceros <- (nrow(filter(fitfinal1,A == "0"))) + (nrow(filter(fitfinal1,B == "0"))) + 
  (nrow(filter(fitfinal1,C == "0"))) + (nrow(filter(fitfinal1,D == "0"))) + 
  (nrow(filter(fitfinal1,E == "0")))
porcentajeCeros <- ceros/(5*(nmuestras-1)*2)*100


colnames(fitfinal1) <- c("Firmicutes","Bacteroidetes",
                         "Actinobacteria","Fusobacteria",
                         "Proteobacteria", "tags")

write.csv2(fitfinal1,"D:\\phylum1.csv", row.names = FALSE)


