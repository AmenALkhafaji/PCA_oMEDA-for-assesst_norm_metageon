
rm(list=ls())
library(tidyr)
library(dplyr)

#    echosystems healthy/unhealthy 

ecosistema1 <- c( c(rep("A",0.15 *1000)) ,c(rep("B",0.6*1000)),
                  c(rep("C",1.5*1000)), c(rep("D",1.6*1000)),
                  c(rep("E",12.2*1000)), c(rep("FF",4.1*1000)),
                  c(rep("Gg",30.8*1000)), c(rep("H",2.3*1000)),
                  c(rep("I",0.14*1000)), c(rep("J",1.21*1000)),
                  c(rep("K",0.18*1000)), c(rep("L",0.32*1000)),
                  c(rep("M",3.31*1000)), c(rep("N",0.63*1000)),
                  c(rep("O",0.17*1000)), c(rep("P",0.34*1000)),
                  c(rep("Q",1.08*1000)), c(rep("R",0.14*1000)),
                  c(rep("S",2.15*1000)), c(rep("TT",2.37*1000)),
                  c(rep("U",0.033*1000)), c(rep("V",0.47*1000)),
                  c(rep("W",0.34*1000)), c(rep("X",5.46*1000)),
                  c(rep("Y",3.47*1000)), c(rep("Z",6.42*1000)),
                  c(rep("pp",4.25*1000)),
                  c(rep("ppA",2.81*1000)),c(rep("ppB",0.36*1000)),
                  c(rep("ppC",1.02*1000)),c(rep("ppD",0.83*1000)),
                  c(rep("ppE",1.80*1000)),
                  c(rep("ppFF",0.31*1000)),c(rep("pptr",0.001*1000)),c(rep("ppH",0.17*1000)),c(rep("ppI",7.1*1000)))

#    ECOSISTEMA ENFERMO

#Firmicutes <- 48%       (A)
#Bacteroidetes <- 22%    (B)
#Actinobacteria <- 12%   (C)
#Proteobacteria <- 12%   (E)
#Fusobacteria <- 6%      (D)

ecosistema2 <- c( c(rep("A",0.38*1000)) ,c(rep("B",0.6*1000)),
                  c(rep("C",1.5*1000)), c(rep("D",0.001*1000)),
                  c(rep("E",26*1000)), c(rep("FF",0.01*1000)),
                  c(rep("Gg",0.02*1000)), c(rep("H",11*1000)),
                  c(rep("I",0.7*1000)), c(rep("J",1*1000)),
                  c(rep("K",0.14*1000)), c(rep("L",1*1000)),
                  c(rep("M",0.15*1000)), c(rep("N",0.3*1000)),
                  c(rep("O",0.001*1000)), c(rep("P",0.3*1000)),
                  c(rep("Q",3*1000)), c(rep("R",0.14*1000)),
                  c(rep("S",0.02*1000)), c(rep("TT",0.5*1000)),
                  c(rep("U",1*1000)), c(rep("V",0.14*1000)),
                  c(rep("W",2*1000)), c(rep("X",4.8*1000)),
                  c(rep("Y",0.001*1000)), c(rep("Z",0.2*1000)),
                  c(rep("pp",0.05*1000)),
                  c(rep("ppA",5*1000)),c(rep("ppB",3)),
                  c(rep("ppC",2.7*1000)),c(rep("ppD",0.8*1000)),
                  c(rep("ppE",8*1000)),c(rep("ppFF",0.22*1000)),c(rep("pptr",0.001*1000))
        ,c(rep("ppH",0.31*1000)),c(rep("ppI",1.23*1000)))



#    ECOSISTEMA ENFERMO

#Firmicutes <- 48%       (A)
#Bacteroidetes <- 22%    (B)
#Actinobacteria <- 12%   (C)
#Proteobacteria <- 12%   (E)
#Fusobacteria <- 6%      (D)



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
fitfinal1$FF <- 0
fitfinal1$G <- 0
fitfinal1$H <- 0
fitfinal1$I <- 0
fitfinal1$Gg <- 0
fitfinal1$K <- 0
fitfinal1$L <- 0
fitfinal1$M <- 0
fitfinal1$N <- 0
fitfinal1$O <- 0
fitfinal1$P <- 0
fitfinal1$Q <- 0
fitfinal1$R <- 0
fitfinal1$S <- 0
fitfinal1$TT <- 0
fitfinal1$U <- 0
fitfinal1$V <- 0
fitfinal1$W <- 0
fitfinal1$X <- 0
fitfinal1$Y <- 0
fitfinal1$Z <- 0
fitfinal1$pp <- 0

fitfinal1$ppA<-0
fitfinal1$ppB<-0
fitfinal1$ppC<-0
fitfinal1$ppD<-0
fitfinal1$ppE<-0
fitfinal1$ppFF<-0
fitfinal1$pptr<-0
fitfinal1$ppH<-0
fitfinal1$ppI<-0



fitfinal1 <- data.frame(fitfinal1, stringsAsFactors = T)

nmuestras <- 151  #ej: 21 son 20 muestras 1 y 20 muestras 2 = 40 muestras total 
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
  fitfinal1$FF <- rbind(fitfinal1$FF,fit[6])
  fitfinal1$Gg <- rbind(fitfinal1$Gg,fit[7])
  fitfinal1$H <- rbind(fitfinal1$H,fit[8])
  fitfinal1$I <- rbind(fitfinal1$I,fit[9])
  fitfinal1$G <- rbind(fitfinal1$J,fit[10])
  fitfinal1$K <- rbind(fitfinal1$K,fit[11])
  fitfinal1$L <- rbind(fitfinal1$L,fit[12])
  fitfinal1$M <- rbind(fitfinal1$M,fit[13])
  fitfinal1$N <- rbind(fitfinal1$N,fit[14])
  fitfinal1$O <- rbind(fitfinal1$O,fit[15])
  fitfinal1$P <- rbind(fitfinal1$P,fit[16])
  fitfinal1$Q <- rbind(fitfinal1$Q,fit[17])
  fitfinal1$R <- rbind(fitfinal1$R,fit[18])
  fitfinal1$S <- rbind(fitfinal1$S,fit[19])
  fitfinal1$TT <-rbind(fitfinal1$TT,fit[20])
  fitfinal1$U <- rbind(fitfinal1$U,fit[21])
  fitfinal1$V <- rbind(fitfinal1$V,fit[22])
  fitfinal1$W <- rbind(fitfinal1$W,fit[23])
  fitfinal1$X <- rbind(fitfinal1$X,fit[24])
  fitfinal1$Y <- rbind(fitfinal1$Y,fit[25])
  fitfinal1$Z <- rbind(fitfinal1$Z,fit[26])
  fitfinal1$pp <- rbind(fitfinal1$pp,fit[27])
  
  fitfinal1$ppA<-rbind(fitfinal1$ppA,fit[28])
  fitfinal1$ppB<-rbind(fitfinal1$ppB,fit[29])
  fitfinal1$ppC<-rbind(fitfinal1$ppC,fit[30])
  fitfinal1$ppD<-rbind(fitfinal1$ppD,fit[31])
  fitfinal1$ppE<-rbind(fitfinal1$ppE,fit[32])
  fitfinal1$ppFF<-rbind(fitfinal1$ppFF,fit[33])
  fitfinal1$pptr<-rbind(fitfinal1$pptr,fit[34])
  fitfinal1$ppH<-rbind(fitfinal1$ppH,fit[35])
  fitfinal1$ppI<-rbind(fitfinal1$ppI,fit[36])
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
  fitfinal1$FF <- rbind(fitfinal1$FF,fit[6])
  fitfinal1$Gg <- rbind(fitfinal1$Gg,fit[7])
  fitfinal1$H <- rbind(fitfinal1$H,fit[8])
  fitfinal1$I <- rbind(fitfinal1$I,fit[9])
  fitfinal1$G <- rbind(fitfinal1$J,fit[10])
  fitfinal1$K <- rbind(fitfinal1$K,fit[11])
  fitfinal1$L <- rbind(fitfinal1$L,fit[12])
  fitfinal1$M <- rbind(fitfinal1$M,fit[13])
  fitfinal1$N <- rbind(fitfinal1$N,fit[14])
  fitfinal1$O <- rbind(fitfinal1$O,fit[15])
  fitfinal1$P <- rbind(fitfinal1$P,fit[16])
  fitfinal1$Q <- rbind(fitfinal1$Q,fit[17])
  fitfinal1$R <- rbind(fitfinal1$R,fit[18])
  fitfinal1$S <- rbind(fitfinal1$S,fit[19])
  fitfinal1$TT <-rbind(fitfinal1$TT,fit[20])
  fitfinal1$U <- rbind(fitfinal1$U,fit[21])
  fitfinal1$V <- rbind(fitfinal1$V,fit[22])
  fitfinal1$W <- rbind(fitfinal1$W,fit[23])
  fitfinal1$X <- rbind(fitfinal1$X,fit[24])
  fitfinal1$Y <- rbind(fitfinal1$Y,fit[25])
  fitfinal1$Z <- rbind(fitfinal1$Z,fit[26])
  fitfinal1$pp <- rbind(fitfinal1$pp,fit[27])
  
  fitfinal1$ppA<-rbind(fitfinal1$ppA,fit[28])
  fitfinal1$ppB<-rbind(fitfinal1$ppB,fit[29])
  fitfinal1$ppC<-rbind(fitfinal1$ppC,fit[30])
  fitfinal1$ppD<-rbind(fitfinal1$ppD,fit[31])
  fitfinal1$ppE<-rbind(fitfinal1$ppE,fit[32])
  fitfinal1$ppFF<-rbind(fitfinal1$ppFF,fit[33])
  fitfinal1$pptr<-rbind(fitfinal1$pptr,fit[34])
  fitfinal1$ppH<-rbind(fitfinal1$ppH,fit[35])
  fitfinal1$ppI<-rbind(fitfinal1$ppI,fit[36])
  n <- n+1
}
fitfinal1 <- data.frame(fitfinal1$A,fitfinal1$B,fitfinal1$C,fitfinal1$D, fitfinal1$E,fitfinal1$FF ,fitfinal1$Gg ,fitfinal1$H ,fitfinal1$I ,fitfinal1$G ,fitfinal1$K ,
                        fitfinal1$L ,fitfinal1$M ,fitfinal1$N ,fitfinal1$O ,fitfinal1$P ,fitfinal1$Q ,fitfinal1$R ,fitfinal1$S ,fitfinal1$TT ,fitfinal1$U ,fitfinal1$V ,fitfinal1$W ,fitfinal1$X ,fitfinal1$Y ,fitfinal1$Z ,fitfinal1$pp,fitfinal1$ppA,fitfinal1$ppB,fitfinal1$ppC,fitfinal1$ppD,fitfinal1$ppE,fitfinal1$ppFF,fitfinal1$pptr,fitfinal1$ppH,fitfinal1$ppI )
               

fitfinal1 <- fitfinal1[-1,]

fitfinal1$tags <- cbind(c(rep("1", nmuestras-1), rep("2", nmuestras-1)))
fitfinal1[,"tags"]<-factor(fitfinal1[,"tags"]) #tags to FACTORS


colnames(fitfinal1) <- c("Bifidobacterium","Butyricimonas","Odoribacter","Paraprevotella",
                          "Bacteroides","Parabacteroides","Prevotella","unkwon1","unkwon2",
                          "unkwon3", "Clostridium","unkwon4","unkwon5","Ruminococcus","Anaerostipes","Blautia",
                           "Coprococcus","Dorea","Lachnospira","Roseburia","unkwon6",
                           "unkwon7","unkwon8","unkwon9","Faecalibacterium","Oscillospira",
                             "Ruminococcus","Dialister","unkwon10","unkwon11","unkwon12","Sutterella",
                             "Escherichia","Klebsiella","Haemophilus","Akkermansia","Tags")

write.csv2(fitfinal1,"F:\\genus.csv", row.names = FALSE)

