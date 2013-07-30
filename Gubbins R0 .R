##Recreating Gubbins et al 2008
##Bluetongue R0 Model with Two Hosts and One Vector
##Jessica Rowland April 2013

##Testing Commit on 18July2013 JR

##Load Required Packages

##Package for Latin Hypercube Sampling
library("lhs")

## Code Outline
## 1. LHS Sampling and Temperature Calibration
## 2. Parameters with LHS Output
## 3. R0 Model Components
## 4. R0 Models
## 5. Visualization of R0 Parameters and R0 Estimates
## 6. Partial Rank Correlation Coefficient

##Latin Hypercube Sampling of Parameters

##Setting the Seed for the Random Number Generator

set.seed(1)

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

## Setting Temperature Value for Temperature functions

temp.LHS <- LHStable[,8]

## Biting Rate Temperature function

temp.func.br <- function(temp,a=.0002,b=3.7,c=41.9,d=1/27){
  bitingrate.temp <- a*temp*(temp-b)*(c-temp)^(d)
  return(bitingrate.temp=bitingrate.temp)
}

##Add values from function to table of values
LHStable[,9] <-temp.func.br(temp=temp.LHS)

## EIP Temperature function


temp.func.EIP <- function(temp,e=1,f=.0003,g=10.4){
  EIP.temp <- (e/(f*(temp)*(temp-g)))
  EIP.temp[EIP.temp<0] <- Inf
  return(EIP.temp=EIP.temp)
}

##Add values from function to table of values
LHStable[,10] <-temp.func.EIP(temp=temp.LHS)

## Vector Death Rate Temperature function

temp.func.d.vec <- function(temp,h=.009,i=.16){
  d.vector.temp <- h*exp(i*temp)
  return(d.vector.temp=d.vector.temp)
}

##Add values from function to table of values
LHStable[,11] <-temp.func.d.vec(temp=temp.LHS)


## Range of Temperature Calibration Parameters 
temperature <-c(0:35)

bitingrate.temp.range <- temp.func.br(temp=temperature)
EIP.temp.range <- (e/(f*(temperature)*(temperature-g)))
d.vector.temp.range <- h*exp(i*temperature)

plot(temperature, bitingrate.temp.range)
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


## Partial Rank Correlation Coefficient

## Rank Transform Parameters

rank.beta.vector <- rank(beta.vector)
rank.d.sheep <- rank(d.sheep)
rank.d.cattle <-rank(d.cattle)
rank.beta.host <- rank(beta.host)
rank.d.vector <- rank(d.vector)
rank.bite <- rank(bite)
rank.EIP.recip <-  rank(EIP.recip)
rank.v.cattle <- rank(v.cattle)
rank.v.sheep <- rank(v.sheep)
rank.pb.cattle <- rank(pb.cattle)
rank.pb.sheep <- rank(pb.sheep)
rank.viremia.cattle.recip <- rank(viremia.cattle.recip)
rank.viremia.sheep.recip <- rank(viremia.sheep.recip)
rank.R0.cattle <- rank(R0.cattle)

## Linear Regression of each input parameter
## Cannot include d.sheep and d.cattle because they were not part of LHS
## Excluded sheep parameters because only looking at R0.cattle

## Beta.vector Input
reg.beta.vector <- lm(rank.beta.vector ~ 
                      rank.beta.host +
                      rank.d.vector +
                      rank.bite +
                      rank.EIP.recip +
                      rank.v.cattle +
                      rank.pb.cattle +
                      rank.viremia.cattle.recip)

## Residuals Beta.vetor Input
resid.beta.vector <- residuals(reg.beta.vector)

## Beta.vector Output

reg.R0.cattle.beta.vector <- lm(rank.R0.cattle ~ 
                                  rank.beta.host +
                                  rank.d.vector +
                                  rank.bite +
                                  rank.EIP.recip +
                                  rank.v.cattle +
                                  rank.pb.cattle +
                                  rank.viremia.cattle.recip)

## Residuals Beta.vector Output
resid.R0.cattle.beta.vector <- residuals(reg.R0.cattle.beta.vector)

## Plot Residuals Beta.vector

plot(resid.beta.vector,resid.R0.cattle.beta.vector)

## PRCC for Beta.vector

Cor.beta.vector <- cor(resid.beta.vector,resid.R0.cattle.beta.vector)


## Beta.host Input
reg.beta.host <- lm( rank.beta.host~  rank.beta.vector+
                        rank.d.vector +
                        rank.bite +
                        rank.EIP.recip +
                        rank.v.cattle +
                        rank.pb.cattle +
                        rank.viremia.cattle.recip)

## Residuals Beta.host Input
resid.beta.host <- residuals(reg.beta.host)

## Beta.host Output

reg.R0.cattle.beta.host <- lm(rank.R0.cattle ~ 
                                  rank.beta.vector +
                                  rank.d.vector +
                                  rank.bite +
                                  rank.EIP.recip +
                                  rank.v.cattle +
                                  rank.pb.cattle +
                                  rank.viremia.cattle.recip)

## Residuals Beta.host Output
resid.R0.cattle.beta.host <- residuals(reg.R0.cattle.beta.host)

## Plot Residuals Beta.host

plot(resid.beta.host,resid.R0.cattle.beta.host)

## PRCC for Beta.host

Cor.beta.host <- cor(resid.beta.host,resid.R0.cattle.beta.host)



## D.vector Input
reg.d.vector <- lm(rank.d.vector ~  rank.beta.host+
                      rank.beta.vector+
                       rank.bite +
                       rank.EIP.recip +
                       rank.v.cattle +
                       rank.pb.cattle +
                       rank.viremia.cattle.recip)

## Residuals d.vector Input
resid.d.vector <- residuals(reg.d.vector)

## d.vector Output

reg.R0.cattle.d.vector<- lm(rank.R0.cattle ~ 
                                rank.beta.vector +
                                rank.beta.host+
                                rank.bite +
                                rank.EIP.recip +
                                rank.v.cattle +
                                rank.pb.cattle +
                                rank.viremia.cattle.recip)

## Residuals d.vector Output
resid.R0.cattle.d.vector <- residuals(reg.R0.cattle.d.vector)

## Plot Residuals d.vector

plot(resid.d.vector,resid.R0.cattle.d.vector)

## PRCC for d.vector

Cor.d.vector <- cor(resid.d.vector,resid.R0.cattle.d.vector)

## bite Input
reg.bite <- lm(rank.bite ~  rank.beta.host+
                     rank.beta.vector+
                     rank.d.vector+
                     rank.EIP.recip +
                     rank.v.cattle +
                     rank.pb.cattle +
                     rank.viremia.cattle.recip)

## Residuals bite Input
resid.bite <- residuals(reg.bite)

## bite Output

reg.R0.cattle.bite<- lm(rank.R0.cattle ~ 
                              rank.beta.vector +
                              rank.beta.host+
                              rank.d.vector
                              rank.EIP.recip +
                              rank.v.cattle +
                              rank.pb.cattle +
                              rank.viremia.cattle.recip)

## Residuals bite Output
resid.R0.cattle.bite <- residuals(reg.R0.cattle.bite)

## Plot Residuals bite

plot(resid.bite,resid.R0.cattle.bite)

## PRCC for bite

Cor.bite <- cor(resid.bite,resid.R0.cattle.bite)


## EIP.recip Input
reg.EIP.recip <- lm(rank.EIP.recip ~  rank.beta.host+
                 rank.beta.vector+
                 rank.d.vector+
                 rank.bite  +
                 rank.v.cattle +
                 rank.pb.cattle +
                 rank.viremia.cattle.recip)

## Residuals EIP.recip Input
resid.EIP.recip <- residuals(reg.EIP.recip)

## EIP.recip Output

reg.R0.cattle.EIP.recip<- lm(rank.R0.cattle ~ 
                          rank.beta.vector +
                          rank.beta.host+
                          rank.d.vector +
                          rank.bite  +
                          rank.v.cattle +
                          rank.pb.cattle +
                          rank.viremia.cattle.recip)

## Residuals EIP.recip Output
resid.R0.cattle.EIP.recip <- residuals(reg.R0.cattle.EIP.recip)

## Plot Residuals EIP.recip

plot(resid.EIP.recip,resid.R0.cattle.EIP.recip)

## PRCC for EIP.recip

Cor.EIP.recip <- cor(resid.EIP.recip,resid.R0.cattle.EIP.recip)


## v.cattle Input
reg.v.cattle <- lm(rank.v.cattle ~  rank.beta.host+
                      rank.beta.vector+
                      rank.d.vector+
                      rank.bite  +
                      rank.EIP.recip +
                      rank.pb.cattle +
                      rank.viremia.cattle.recip)

## Residuals v.cattle Input
resid.v.cattle <- residuals(reg.v.cattle)

## v.cattle Output

reg.R0.cattle.v.cattle <- lm(rank.R0.cattle ~ 
                               rank.beta.vector +
                               rank.beta.host+
                               rank.d.vector +
                               rank.bite  +
                               rank.EIP.recip +
                               rank.pb.cattle +
                               rank.viremia.cattle.recip)

## Residuals v.cattle Output
resid.R0.cattle.v.cattle <- residuals(reg.R0.cattle.v.cattle)

## Plot Residuals v.cattle

plot(resid.v.cattle,resid.R0.cattle.v.cattle)

## PRCC for v.cattle

Cor.v.cattle <- cor(resid.v.cattle,resid.R0.cattle.v.cattle)


## pb.cattle Input
reg.pb.cattle <- lm( rank.pb.cattle ~  rank.beta.host+
                     rank.beta.vector+
                     rank.d.vector+
                     rank.bite  +
                     rank.EIP.recip +
                     rank.v.cattle  +
                     rank.viremia.cattle.recip)

## Residuals pb.cattle Input
resid.pb.cattle <- residuals(reg.pb.cattle)

## pb.cattle Output

reg.R0.cattle.pb.cattle <- lm(rank.R0.cattle ~ 
                               rank.beta.vector +
                               rank.beta.host+
                               rank.d.vector +
                               rank.bite  +
                               rank.EIP.recip +
                               rank.v.cattle +
                               rank.viremia.cattle.recip)

## Residuals pb.cattle Output
resid.R0.cattle.pb.cattle <- residuals(reg.R0.cattle.pb.cattle)

## Plot Residuals pb.cattle

plot(resid.pb.cattle,resid.R0.cattle.pb.cattle)

## PRCC for pb.cattle

Cor.pb.cattle <- cor(resid.pb.cattle,resid.R0.cattle.pb.cattle)

## viremia.cattle.recip Input
reg.viremia.cattle.recip <- lm( rank.viremia.cattle.recip~  rank.beta.host+
                      rank.beta.vector+
                      rank.d.vector+
                      rank.bite  +
                      rank.EIP.recip +
                      rank.v.cattle +
                      rank.pb.cattle)

## Residuals viremia.cattle.recip Input
resid.viremia.cattle.recip <- residuals(reg.viremia.cattle.recip)

## viremia.cattle.recip Output

reg.R0.cattle.viremia.cattle.recip <- lm(rank.R0.cattle ~ 
                                rank.beta.vector +
                                rank.beta.host+
                                rank.d.vector +
                                rank.bite  +
                                rank.EIP.recip +
                                rank.pb.cattle +
                                rank.v.cattle)

## Residuals viremia.cattle.recip Output
resid.R0.cattle.viremia.cattle.recip <- residuals(reg.R0.cattle.viremia.cattle.recip)

## Plot Residuals viremia.cattle.recip

plot(resid.viremia.cattle.recip,resid.R0.cattle.viremia.cattle.recip)

## PRCC for viremia.cattle.recip

Cor.viremia.cattle.recip <- cor(resid.viremia.cattle.recip,resid.R0.cattle.viremia.cattle.recip)