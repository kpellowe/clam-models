##################################################
# LIFE HISTORY PARAMETERIZATION AND YPR AND SSB/R FUNCTION USED IN PELLOWE & LESLIE ECOSPHERE MANUSCRIPT
library("MASS")
mu1 <- 85.7504
mu2 <- 0.1215
sigma1 <- ((2.302)^2)*311
sigma2 <- ((0.016)^2)*311
X1 <- -0.6 #estimate from Hart 2003
Linfk <- mvrnorm(n=10000, mu=c(mu1, mu2), Sigma = matrix(c(sigma1,X1, 
                                                           X1, sigma2), ncol = 2, byrow = TRUE), empirical = TRUE)
library("truncnorm")
distM = rtruncnorm(10000, a=0, b=1, mean = 0.3538, sd = (0.068*sqrt(18)))
dista = rnorm(10000, mean = 0.0021, sd = (0.09*sqrt(2485)))
distb = rnorm(10000, mean = 2.5445, sd = (0.021*sqrt(2485)))

# YPR FUNCTION
yprfunc <- function(SL50) #YPR is a function of length at 50% selectivity. This is how we modeled different slot limit treatments.
{
  #Create matrix to store results
  Datos2 <- matrix(1:440000, nrow = 110000) 
  
  for (uo in 1:10000) #10,000 Monte Carlo Simulations
  {
    #set variables that with change with each iteration via monte carlo simulation
    #select Linf and k
    Array <- Linfk[sample.int(nrow(Linfk), size=1), 1:2]
    while (Array[1]<60 | Array[1]>120 | Array[2]<0 | Array[2]> 0.2) 
    {
      Array <- Linfk[sample.int(nrow(Linfk), size=1), 1:2]
    }
    Linf <- Array[1]
    #Linf
    k <- Array[2]
    #k
    #select M
    M <- sample(distM, size = 1, replace = TRUE)
    #M
    #select a and b. Some code to try to control for crazy a and b parameters, which can result in unrealistic weights
    Xa = matrix(ncol = 1, nrow=1000)
    for(year in 0:1) {
      #for each Fm value sample the ypr 10000 times
      for(i in 1:1000){
        Xa[i,year] = sample(dista, 1) }}
    a <- apply(Xa,2,mean)
    Xb = matrix(ncol = 1, nrow=1000)
    for(year in 0:1) {
      #for each Fm value sample the ypr 10000 times
      for(i in 1:1000){
        Xb[i,year] = sample(distb, 1) }}
    b <- apply(Xb,2,mean)
    
    #Now set variables that will not change with each iteration
    Sp <- 0.4782
    #SL50 is set in function call
    SL50two <- 90+8.5273 #set according to upper limit. 8.5273mm above upper limit
    m <- 0.5 #estimate from Villalejo-Fuerte et al. study
    L50 <- 42 #estimate from Villalejo-Fuerte et al. study
    
    #set up steps of function
    for(Fm in c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))
    {
      #At beginning of forloop (t=0), set these parameters to zero
      yprprev <- 0
      ypr <- 0
      SSBR <- 0
      SSBRprev <- 0
      Selectiv <- 0
      AccMort <- 0
      Length <- 0
      SumS <- 0
      DYDF <- 0
      DYDFprev <- 0
      
      
      for(t in 0:11)  #Approx maximum lifespan of a Mexican chocolate clam
      {
        #YPR & SSB/R calculations step-by-step
        Length <- Length+(t+((Linf-t)*(1-exp(-k*1))))
        
        #Stop loop if Length grows larger than 110mm, the max size observed from fishers' catch
        if (Length > 110) {
          break
        }
        
        Weight <- a*(Length^b)
        if(t<1){SelectivPrev <- 0} else {SelectivPrev <- Selectiv }
        #Selectivity for only min size. Must choose one depending on slot limit scenario
        #Selectiv <- 1/(1+exp(-Sp*(Length-SL50))) #for only min size
        
        # Selectivity for both min and max size
        Selectiv <- (1-1/(1+exp(-Sp*(Length-SL50two))))/(1+exp(-Sp*(Length-SL50))) #for min & max size
        if(t < 1){AccMort = 0} else {AccMort <- AccMort +(Selectiv*Fm)+M}
        #AccMort
        yprprev <- ypr + yprprev
        # yprprev
        ypr <- Weight*Fm*Selectiv/(Fm*Selectiv+M)*(1-exp(-(Fm*Selectiv+M)))*exp(-AccMort)
        #set YPRsum to previous timestep YPR + current timestep YPR
        yprsum <- ypr + yprprev

        ## SSB/R calculations
        Pmature <- 1/(1+exp(-m*(Length-L50)))
        PropFem = 0.5
        SSBRprev <- SSBR + SSBRprev
        SSBR <- Pmature*PropFem*exp(-AccMort)*Weight
        SSBRsum <- SSBR + SSBRprev
        #SSBRsum
        
      }
      ###########################################################
      Datos2[uo+((uo-1)*10)+(Fm*10)] <- uo
      Datos2[uo+((uo-1)*10)+(Fm*10)+110000] <- Fm
      Datos2[uo+((uo-1)*10)+(Fm*10)+220000] <- yprsum
      Datos2[uo+((uo-1)*10)+(Fm*10)+330000] <- SSBRsum 
    }}
  #Print dataframe as function output
  return(Datos2)
}


####### MIN SIZE 44MM ##################
DataS44.df <- as.data.frame(yprfunc(35.4727))
write.csv(DataS44.df,'DataS44.csv')
setwd("~/Dropbox/2019 UMaine Postdoc/1_Fall 2019/Ecosphere 2nd Revision")
DataS44 <- read.table("~/Dropbox/2019 UMaine Postdoc/1_Fall 2019/Ecosphere 2nd Revision/DataS44.csv", header=TRUE, sep=",")
DataData44 <- data.frame("F" = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))
DataData44$YPRmean <- tapply(DataS44$V3,DataS44$V2,mean)
DataData44$YPRerr <- tapply(DataS44$V3,DataS44$V2,sd)/sqrt(10000)
DataData44$SSBmean <- tapply(DataS44$V4,DataS44$V2,mean)
DataData44$SSBerr <- tapply(DataS44$V4,DataS44$V2,sd)/sqrt(10000)
DataData44$YPRCIup <- DataData44$YPRmean + (1.96*DataData44$YPRerr)
DataData44$YPRCIlo <- DataData44$YPRmean - (1.96*DataData44$YPRerr)
DataData44$SSBCIup <- DataData44$SSBmean + (1.96*DataData44$SSBerr)
DataData44$SSBCIlo <- DataData44$SSBmean - (1.96*DataData44$SSBerr)
DataData44$Treatment <- "Scenario 1: Min 44"
View(DataData44)

####### MIN SIZE 44MM MAX SIZE 80 ##################
DataS4480.df <- as.data.frame(yprfunc(35.4727))
write.csv(DataS4480.df,'DataS4480.csv')
setwd("~/Dropbox/2019 UMaine Postdoc/1_Fall 2019/Ecosphere 2nd Revision")
DataS4480 <- read.table("~/Dropbox/2019 UMaine Postdoc/1_Fall 2019/Ecosphere 2nd Revision/DataS4480.csv", header=TRUE, sep=",")
DataData4480 <- data.frame("F" = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))
DataData4480$YPRmean <- tapply(DataS4480$V3,DataS4480$V2,mean)
DataData4480$YPRerr <- tapply(DataS4480$V3,DataS4480$V2,sd)/sqrt(10000)
DataData4480$SSBmean <- tapply(DataS4480$V4,DataS4480$V2,mean)
DataData4480$SSBerr <- tapply(DataS4480$V4,DataS4480$V2,sd)/sqrt(10000)
DataData4480$YPRCIup <- DataData4480$YPRmean + (1.96*DataData4480$YPRerr)
DataData4480$YPRCIlo <- DataData4480$YPRmean - (1.96*DataData4480$YPRerr)
DataData4480$SSBCIup <- DataData4480$SSBmean + (1.96*DataData4480$SSBerr)
DataData4480$SSBCIlo <- DataData4480$SSBmean - (1.96*DataData4480$SSBerr)
DataData4480$Treatment <- "Scenario 2: Min 44, Max 80"
View(DataData4480)

####### MIN SIZE 44MM MAX SIZE 90 ##################
DataS4490.df <- as.data.frame(yprfunc(35.4727))
write.csv(DataS4490.df,'DataS4490.csv')
setwd("~/Dropbox/2019 UMaine Postdoc/1_Fall 2019/Ecosphere 2nd Revision")
DataS4490 <- read.table("~/Dropbox/2019 UMaine Postdoc/1_Fall 2019/Ecosphere 2nd Revision/DataS4490.csv", header=TRUE, sep=",")
DataData4490 <- data.frame("F" = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))
DataData4490$YPRmean <- tapply(DataS4490$V3,DataS4490$V2,mean)
DataData4490$YPRerr <- tapply(DataS4490$V3,DataS4490$V2,sd)/sqrt(10000)
DataData4490$SSBmean <- tapply(DataS4490$V4,DataS4490$V2,mean)
DataData4490$SSBerr <- tapply(DataS4490$V4,DataS4490$V2,sd)/sqrt(10000)
DataData4490$YPRCIup <- DataData4490$YPRmean + (1.96*DataData4490$YPRerr)
DataData4490$YPRCIlo <- DataData4490$YPRmean - (1.96*DataData4490$YPRerr)
DataData4490$SSBCIup <- DataData4490$SSBmean + (1.96*DataData4490$SSBerr)
DataData4490$SSBCIlo <- DataData4490$SSBmean - (1.96*DataData4490$SSBerr)
DataData4490$Treatment <- "Scenario 3: Min 44, Max 90"
View(DataData4490)


####### MIN SIZE 54MM ##################
DataS54.df <- as.data.frame(yprfunc(45.4727))
write.csv(DataS54.df,'DataS54.csv')
setwd("~/Dropbox/2019 UMaine Postdoc/1_Fall 2019/Ecosphere 2nd Revision")
DataS54 <- read.table("~/Dropbox/2019 UMaine Postdoc/1_Fall 2019/Ecosphere 2nd Revision/DataS54.csv", header=TRUE, sep=",")
DataData54 <- data.frame("F" = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))
DataData54$YPRmean <- tapply(DataS54$V3,DataS54$V2,mean)
DataData54$YPRerr <- tapply(DataS54$V3,DataS54$V2,sd)/sqrt(10000)
DataData54$SSBmean <- tapply(DataS54$V4,DataS54$V2,mean)
DataData54$SSBerr <- tapply(DataS54$V4,DataS54$V2,sd)/sqrt(10000)
DataData54$YPRCIup <- DataData54$YPRmean + (1.96*DataData54$YPRerr)
DataData54$YPRCIlo <- DataData54$YPRmean - (1.96*DataData54$YPRerr)
DataData54$SSBCIup <- DataData54$SSBmean + (1.96*DataData54$SSBerr)
DataData54$SSBCIlo <- DataData54$SSBmean - (1.96*DataData54$SSBerr)
DataData54$Treatment <- "Scenario 4: Min 54"
View(DataData54)

####### MIN SIZE 54MM MAX SIZE 80 ##################
DataS5480.df <- as.data.frame(yprfunc(45.4727))
write.csv(DataS5480.df,'DataS5480.csv')
setwd("~/Dropbox/2019 UMaine Postdoc/1_Fall 2019/Ecosphere 2nd Revision")
DataS5480 <- read.table("~/Dropbox/2019 UMaine Postdoc/1_Fall 2019/Ecosphere 2nd Revision/DataS5480.csv", header=TRUE, sep=",")
DataData5480 <- data.frame("F" = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))
DataData5480$YPRmean <- tapply(DataS5480$V3,DataS5480$V2,mean)
DataData5480$YPRerr <- tapply(DataS5480$V3,DataS5480$V2,sd)/sqrt(10000)
DataData5480$SSBmean <- tapply(DataS5480$V4,DataS5480$V2,mean)
DataData5480$SSBerr <- tapply(DataS5480$V4,DataS5480$V2,sd)/sqrt(10000)
DataData5480$YPRCIup <- DataData5480$YPRmean + (1.96*DataData5480$YPRerr)
DataData5480$YPRCIlo <- DataData5480$YPRmean - (1.96*DataData5480$YPRerr)
DataData5480$SSBCIup <- DataData5480$SSBmean + (1.96*DataData5480$SSBerr)
DataData5480$SSBCIlo <- DataData5480$SSBmean - (1.96*DataData5480$SSBerr)
DataData5480$Treatment <- "Scenario 5: Min 54, Max 80"
View(DataData5480)

####### MIN SIZE 54MM MAX SIZE 90 ##################
DataS5490.df <- as.data.frame(yprfunc(45.4727))
write.csv(DataS5490.df,'DataS5490.csv')
setwd("~/Dropbox/2019 UMaine Postdoc/1_Fall 2019/Ecosphere 2nd Revision")
DataS5490 <- read.table("~/Dropbox/2019 UMaine Postdoc/1_Fall 2019/Ecosphere 2nd Revision/DataS5490.csv", header=TRUE, sep=",")
DataData5490 <- data.frame("F" = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))
DataData5490$YPRmean <- tapply(DataS5490$V3,DataS5490$V2,mean)
DataData5490$YPRerr <- tapply(DataS5490$V3,DataS5490$V2,sd)/sqrt(10000)
DataData5490$SSBmean <- tapply(DataS5490$V4,DataS5490$V2,mean)
DataData5490$SSBerr <- tapply(DataS5490$V4,DataS5490$V2,sd)/sqrt(10000)
DataData5490$YPRCIup <- DataData5490$YPRmean + (1.96*DataData5490$YPRerr)
DataData5490$YPRCIlo <- DataData5490$YPRmean - (1.96*DataData5490$YPRerr)
DataData5490$SSBCIup <- DataData5490$SSBmean + (1.96*DataData5490$SSBerr)
DataData5490$SSBCIlo <- DataData5490$SSBmean - (1.96*DataData5490$SSBerr)
DataData5490$Treatment <- "Scenario 6: Min 54, Max 90"
View(DataData5490)

####### MIN SIZE 64MM ##################
DataS64.df <- as.data.frame(yprfunc(55.4727))
write.csv(DataS64.df,'DataS64.csv')
setwd("~/Dropbox/2019 UMaine Postdoc/1_Fall 2019/Ecosphere 2nd Revision")
DataS64 <- read.table("~/Dropbox/2019 UMaine Postdoc/1_Fall 2019/Ecosphere 2nd Revision/DataS64.csv", header=TRUE, sep=",")
View(DataS64)
DataData64 <- data.frame("F" = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))
DataData64$YPRmean <- tapply(DataS64$V3,DataS64$V2,mean)
DataData64$YPRerr <- tapply(DataS64$V3,DataS64$V2,sd)/sqrt(10000)
DataData64$SSBmean <- tapply(DataS64$V4,DataS64$V2,mean)
DataData64$SSBerr <- tapply(DataS64$V4,DataS64$V2,sd)/sqrt(10000)
DataData64$YPRCIup <- DataData64$YPRmean + (1.96*DataData64$YPRerr)
DataData64$YPRCIlo <- DataData64$YPRmean - (1.96*DataData64$YPRerr)
DataData64$SSBCIup <- DataData64$SSBmean + (1.96*DataData64$SSBerr)
DataData64$SSBCIlo <- DataData64$SSBmean - (1.96*DataData64$SSBerr)
DataData64$Treatment <- "Scenario 7: Min 64"
View(DataData64)

####### MIN SIZE 64MM MAX SIZE 80 ##################
DataS6480.df <- as.data.frame(yprfunc(55.4727))
write.csv(DataS6480.df,'DataS6480.csv')
setwd("~/Dropbox/2019 UMaine Postdoc/1_Fall 2019/Ecosphere 2nd Revision")
DataS6480 <- read.table("~/Dropbox/2019 UMaine Postdoc/1_Fall 2019/Ecosphere 2nd Revision/DataS6480.csv", header=TRUE, sep=",")
#View(DataS6480)
#DataS6480 <- DataS6480[ which(DataS6480$V3 >= 0 & DataS6480$V4 > 0), ]
DataData6480 <- data.frame("F" = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))
DataData6480$YPRmean <- tapply(DataS6480$V3,DataS6480$V2,mean)
DataData6480$YPRerr <- tapply(DataS6480$V3,DataS6480$V2,sd)/sqrt(10000)
DataData6480$SSBmean <- tapply(DataS6480$V4,DataS6480$V2,mean)
DataData6480$SSBerr <- tapply(DataS6480$V4,DataS6480$V2,sd)/sqrt(10000)
DataData6480$YPRCIup <- DataData6480$YPRmean + (1.96*DataData6480$YPRerr)
DataData6480$YPRCIlo <- DataData6480$YPRmean - (1.96*DataData6480$YPRerr)
DataData6480$SSBCIup <- DataData6480$SSBmean + (1.96*DataData6480$SSBerr)
DataData6480$SSBCIlo <- DataData6480$SSBmean - (1.96*DataData6480$SSBerr)
DataData6480$Treatment <- "Scenario 8: Min 64, Max 80"
View(DataData6480)

####### MIN SIZE 64MM MAX SIZE 90 ##################
DataS6490.df <- as.data.frame(yprfunc(55.4727))
write.csv(DataS6490.df,'DataS6490.csv')
setwd("~/Dropbox/2019 UMaine Postdoc/1_Fall 2019/Ecosphere 2nd Revision")
DataS6490 <- read.table("~/Dropbox/2019 UMaine Postdoc/1_Fall 2019/Ecosphere 2nd Revision/DataS6490.csv", header=TRUE, sep=",")
DataData6490 <- data.frame("F" = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))
DataData6490$YPRmean <- tapply(DataS6490$V3,DataS6490$V2,mean)
DataData6490$YPRerr <- tapply(DataS6490$V3,DataS6490$V2,sd)/sqrt(10000)
DataData6490$SSBmean <- tapply(DataS6490$V4,DataS6490$V2,mean)
DataData6490$SSBerr <- tapply(DataS6490$V4,DataS6490$V2,sd)/sqrt(10000)
DataData6490$YPRCIup <- DataData6490$YPRmean + (1.96*DataData6490$YPRerr)
DataData6490$YPRCIlo <- DataData6490$YPRmean - (1.96*DataData6490$YPRerr)
DataData6490$SSBCIup <- DataData6490$SSBmean + (1.96*DataData6490$SSBerr)
DataData6490$SSBCIlo <- DataData6490$SSBmean - (1.96*DataData6490$SSBerr)
DataData6490$Treatment <- "Scenario 9: Min 64, Max 90"
View(DataData6490)

####### MIN SIZE 74MM ##################
DataS74.df <- as.data.frame(yprfunc(65.4727))
write.csv(DataS74.df,'DataS74.csv')
setwd("~/Dropbox/2019 UMaine Postdoc/1_Fall 2019/Ecosphere 2nd Revision")
DataS74 <- read.table("~/Dropbox/2019 UMaine Postdoc/1_Fall 2019/Ecosphere 2nd Revision/DataS74.csv", header=TRUE, sep=",")
DataData74 <- data.frame("F" = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))
DataData74$YPRmean <- tapply(DataS74$V3,DataS74$V2,mean)
DataData74$YPRerr <- tapply(DataS74$V3,DataS74$V2,sd)/sqrt(10000)
DataData74$SSBmean <- tapply(DataS74$V4,DataS74$V2,mean)
DataData74$SSBerr <- tapply(DataS74$V4,DataS74$V2,sd)/sqrt(10000)
DataData74$YPRCIup <- DataData74$YPRmean + (1.96*DataData74$YPRerr)
DataData74$YPRCIlo <- DataData74$YPRmean - (1.96*DataData74$YPRerr)
DataData74$SSBCIup <- DataData74$SSBmean + (1.96*DataData74$SSBerr)
DataData74$SSBCIlo <- DataData74$SSBmean - (1.96*DataData74$SSBerr)
DataData74$Treatment <- "Scenario 10: Min 74"
View(DataData74)

####### MIN SIZE 74MM MAX SIZE 80 ##################
DataS7480.df <- as.data.frame(yprfunc(65.4727))
write.csv(DataS7480.df,'DataS7480.csv')
setwd("~/Dropbox/2019 UMaine Postdoc/1_Fall 2019/Ecosphere 2nd Revision")
DataS7480 <- read.table("~/Dropbox/2019 UMaine Postdoc/1_Fall 2019/Ecosphere 2nd Revision/DataS7480.csv", header=TRUE, sep=",")
DataData7480 <- data.frame("F" = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))
DataData7480$YPRmean <- tapply(DataS7480$V3,DataS7480$V2,mean)
DataData7480$YPRerr <- tapply(DataS7480$V3,DataS7480$V2,sd)/sqrt(10000)
DataData7480$SSBmean <- tapply(DataS7480$V4,DataS7480$V2,mean)
DataData7480$SSBerr <- tapply(DataS7480$V4,DataS7480$V2,sd)/sqrt(10000)
DataData7480$YPRCIup <- DataData7480$YPRmean + (1.96*DataData7480$YPRerr)
DataData7480$YPRCIlo <- DataData7480$YPRmean - (1.96*DataData7480$YPRerr)
DataData7480$SSBCIup <- DataData7480$SSBmean + (1.96*DataData7480$SSBerr)
DataData7480$SSBCIlo <- DataData7480$SSBmean - (1.96*DataData7480$SSBerr)
DataData7480$Treatment <- "Scenario 11: Min 74, Max 80"
View(DataData7480)

####### MIN SIZE 74MM MAX SIZE 90 ##################
DataS7490.df <- as.data.frame(yprfunc(65.4727))
write.csv(DataS7490.df,'DataS7490.csv')
setwd("~/Dropbox/2019 UMaine Postdoc/1_Fall 2019/Ecosphere 2nd Revision")
DataS7490 <- read.table("~/Dropbox/2019 UMaine Postdoc/1_Fall 2019/Ecosphere 2nd Revision/DataS7490.csv", header=TRUE, sep=",")
DataData7490 <- data.frame("F" = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))
DataData7490$YPRmean <- tapply(DataS7490$V3,DataS7490$V2,mean)
DataData7490$YPRerr <- tapply(DataS7490$V3,DataS7490$V2,sd)/sqrt(10000)
DataData7490$SSBmean <- tapply(DataS7490$V4,DataS7490$V2,mean)
DataData7490$SSBerr <- tapply(DataS7490$V4,DataS7490$V2,sd)/sqrt(10000)
DataData7490$YPRCIup <- DataData7490$YPRmean + (1.96*DataData7490$YPRerr)
DataData7490$YPRCIlo <- DataData7490$YPRmean - (1.96*DataData7490$YPRerr)
DataData7490$SSBCIup <- DataData7490$SSBmean + (1.96*DataData7490$SSBerr)
DataData7490$SSBCIlo <- DataData7490$SSBmean - (1.96*DataData7490$SSBerr)
DataData7490$Treatment <- "Scenario 12: Min 74, Max 90"
View(DataData7490)


#Combine all dataframes by row using dplyr
library(dplyr)
alldata <- dplyr::bind_rows(DataData44, DataData4480, DataData4490, DataData54, DataData5480, DataData5490, DataData64, DataData6480, DataData6490, DataData74, DataData7480, DataData7490)
#specify order levels so that treatments go from 1 to 12
alldata$Treatment <- factor(alldata$Treatment, levels = c("Scenario 1: Min 44", "Scenario 2: Min 44, Max 80", "Scenario 3: Min 44, Max 90", "Scenario 4: Min 54", "Scenario 5: Min 54, Max 80", "Scenario 6: Min 54, Max 90", "Scenario 7: Min 64", "Scenario 8: Min 64, Max 80", "Scenario 9: Min 64, Max 90", "Scenario 10: Min 74", "Scenario 11: Min 74, Max 80", "Scenario 12: Min 74, Max 90"))
#View data
View(alldata)

