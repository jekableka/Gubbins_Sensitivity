## Set seed for random number generator

set.seed(1)

##Load Hmisc package for correlation coeffecient with p values
library("Hmisc")

##Generate random input variables
x <- rnorm(10)
y <- rnorm(10)
z <- rnorm(10)

##Generate random output variable
w <-rnorm(10)

##Rank transform variables
a <- rank(x)
b <- rank(y)
c <- rank(z)
out <- rank(w)

##PRCC Function
prcc.func <- function(a,b,c,out){
	model.input <- lm(a ~ b + c)
	model.output <- lm(out ~ b + c)
	model.corr <- rcorr((residuals(model.input)),(residuals(model.output)))
	prcc <- model.corr$r[1,2]
	pvalue <- model.corr$P[1,2]
	return(list(prcc=prcc,pvalue=pvalue))
}

prcc.func(a,b,c,out)