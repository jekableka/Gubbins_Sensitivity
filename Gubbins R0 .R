##Recreating Gubbins et al 2008
##Blue Tongue R0 Model with Two Hosts and One Vector
##Jessica Rowland April 2013

##Load Required Packages

##Package for Latin Hypercube Sampling
library("lhs")

## Code Outline
## 1. LHS Sampling and Temperature Calibration
## 2. Parameters with LHS Output
## 3. R0 Model Components
## 4. R0 Models
## 5. Visualization R0 Parameters and R0 Estimates

##Latin Hypercube Sampling of Parameters

LHStest <-randomLHS(100,8)
LHStable <-matrix(0,nrow=100,ncol=11)
LHStable[,1] <-qunif(LHStest[,1],1,100)
LHStable[,2] <-qunif(LHStest[,2],.8,1)
LHStable[,3] <-qunif(LHStest[,3],.001,.15)
LHStable[,4] <-qunif(LHStest[,4],0,5000)
LHStable[,5] <-qunif(LHStest[,5],0,5000)
LHStable[,6] <-qgamma(LHStest[,6],16.4)
LHStable[,7] <-qgamma(LHStest[,7],20.6)
LHStable[,8] <-qunif(LHStest[,8],0,35)

##Temperature calibration of parameters

bitingrate.temp <- .0002*(LHStable[,8])*(LHStable[,8]-3.7)*(41.9-LHStable[,8])^(1/27)
EIP.temp <- (1/(.0003*(LHStable[,8])*(LHStable[,8]-10.4)))
d.vector.temp <- .009*exp(.16*LHStable[,8])

## Adding Temperature Calibrated Parameters to the LHS Matrix

LHStable[,9] <-bitingrate.temp
LHStable[,10] <-EIP.temp
LHStable[,11] <- d.vector.temp

## Range of Temperature Calibration Parameters 
temperature <-c(0:35)

bitingrate.temp.range <- .0002*(temperature)*(temperature-3.7)*(41.9-temperature)^(1/27)
EIP.temp.range <- (1/(.0003*(temperature)*(temperature-10.4)))
d.vector.temp.range <- .009*exp(.16*temperature)

plot(temperature,bitingrate.temp.range)
plot(temperature, EIP.temp.range)
plot(temperature, d.vector.temp.range)


##Parameters


beta.vector <- LHStable[,2] ## probabilty of transmission from vector to host
d.sheep <- .01 ## infection induced mortality in sheep ranges from 0.001 to 0.01
d.cattle <-0 ##infection induced mortaltiy it cattle
beta.host <- LHStable[,3] ## probability of transmission from host to vector
d.vector <- LHStable[,11] ## vector death rate
bite <- LHStable[,9] ## recipricol of days between bites
EIP.recip <-  1/LHStable[,10]##  reciprical of days in EIP
v.cattle <- LHStable[,4] ## ratio of vectors to cattle
v.sheep <- LHStable[,5] ## ratio of vectors to sheep
pb.cattle <- LHStable[,5]/(LHStable[,5]+.1*LHStable[,4]) ## proportion of bites on cattle
pb.sheep <- 1-pb.cattle  ## proportion of bites on sheep
viremia.cattle.recip <- 1/LHStable[,7]  ## recipricol of days of viremia in cattle
viremia.sheep.recip <- 1/LHStable[,6] ## recipricol of days of viremia in sheep

##Bluetongue R0 Model Components

R0.contact <- ((beta.vector*beta.host*(bite^2))/d.vector)  ## contact component

R0.vector  <- ((EIP.recip)/(d.vector+(EIP.recip)))  ## vector component

R0.host.c <- ((v.cattle*(pb.cattle^2))/(viremia.cattle.recip+d.cattle)) ##  catlle component

R0.host.s <- ((v.sheep*(pb.sheep^2))/(viremia.sheep.recip+d.sheep))  ## sheep component

## Bluetongue R0 models for cattle, sheep, and both
R0.cattle <- sqrt(R0.contact*R0.vector*R0.host.c)

R0.sheep <- sqrt(R0.contact*R0.vector*R0.host.s)

R0.both <- sqrt(R0.contact*R0.vector*(R0.host.c + R0.host.s))



## All of the Data in a Dataframe
R0data <- data <-as.data.frame(LHStable)
names(R0data)[1] <-paste("EIPStages")
names(R0data)[2] <-paste("VtoH")
names(R0data)[3] <-paste("HtoV")
names(R0data)[4] <-paste("VtoCattle")
names(R0data)[5] <-paste("VtoSheep")
names(R0data)[6] <-paste("SheepViremia")
names(R0data)[7] <-paste("CattleViremia")
names(R0data)[8] <-paste("Temperature")
names(R0data)[9] <-paste("BitingRate")
names(R0data)[10] <-paste("EIP")
names(R0data)[11] <-paste("VDeathRate")
R0data$Cattle <- R0.cattle
R0data$Sheep <- R0.sheep
R0data$Both <-R0.both

## R0 Estimate Summaries

summary(R0data[,12:14])

## Visualization of Parameters and R0 for Cattle

plot(R0data$EIPStages,R0data$Cattle)
lines(lowess(R0data$EIPStages,R0data$Cattle))


plot(R0data$VtoH,R0data$Cattle)
lines(lowess(R0data$VtoH,R0data$Cattle))

plot(R0data$HtoV,R0data$Cattle)
lines(lowess(R0data$HtoV,R0data$Cattle))

plot(R0data$SheepViremia,R0data$Cattle)
lines(lowess(R0data$SheepViremia,R0data$Cattle))

plot(R0data$CattleViremia,R0data$Cattle)
lines(lowess(R0data$CattleViremia,R0data$Cattle))


plot(R0data$Temperature,R0data$Cattle) 
lines(lowess(R0data$Temperature,R0data$Cattle))

plot(R0data$BitingRate,R0data$Cattle) 
lines(lowess(R0data$BitingRate,R0data$Cattle))

plot(R0data$EIP,R0data$Cattle) 
lines(lowess(R0data$EIP,R0data$Cattle))

plot(R0data$VDeathRate,R0data$Cattle) 
lines(lowess(R0data$VDeathRate,R0data$Cattle))







