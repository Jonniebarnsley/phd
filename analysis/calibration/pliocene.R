setwd("~/Code/Calibration")
load("Grant2020_PAGES_workspace.RData")

library(astrochron)
library(ggplot2)
library(dplyr)

data <- read.csv('Pliocene data.csv')
slc <- na.omit( data$slc_10k_minus_control)

iso=iso(PlioSeaNZ,3032,3310) # isolate record to mPWP
strats(iso) # a strat assessment to determiene sampling variability
lin_PlioSeaNZ=linterp(iso,2) # linear interpolation at mean or median sampling

site<-lin_PlioSeaNZ
outmain<-lin_PlioSeaNZ[1]
n=2 #dt sampling

#### This calculates the minima and maxima at various window durations (Milankovitch periodicities)

for(i in 1:nrow(site)) {
  #period = 20 
  celln <- (20/n) 
  
  low <- min(site[i:(i+celln),2])  # lowest and highest levels for target section/group combination
  high <- max(site[i:(i+celln),2])
  outmain[i,2]<-high-low 
  
}

colnames(outmain)<-c("time.ka","t=20")
PlioSeaNZ_amp<-outmain

par(mfrow=c(1,1))

x.values <- seq(0, 1, by=0.02)

a=4
b=2
beta.dist <- dbeta(x.values, shape1=a, shape2=b)
plot(x.values*7.2, beta.dist, type='l', )

beta.sample <- rbeta(100000, shape1=a, shape2=b)
greenland.sample <- beta.sample*7.2
plot(density(greenland.sample))

t20 <- na.omit(PlioSeaNZ_amp$`t=20`)
quantile(t20)

slc.antarctica <- c(sapply(t20, function(x) x - greenland.sample))

prior.global <- density(t20)
prior.antarctica <- density(slc.antarctica)

plot(prior.global)
plot(prior.antarctica)

mean <- mean(slc.antarctica)
mean
sigma <- sd(slc.antarctica)

sample <- rnorm(10000, mean(slc.antarctica), sd(slc.antarctica))

plot(density(sample))
lines(prior.antarctica)

distances = slc - mean

weightsraw = sapply(distances, function(x) exp( -0.5 * sum( x * x , na.rm=TRUE) / (sigma * sigma)))
weights = weightsraw / sum(weightsraw)

prior <- density(slc)
posterior <- density(slc, weights = weights)


par(mfrow=c(1,1)) 

plot(prior, ylim=c(0, 0.1))
lines( posterior, col = "red" )
lines(prior.antarctica, col='blue')
legend('topright', legend=c('Prior', 'Likelihood', 'Posterior'), col = c('black','blue', 'red'), pch-16, lty=1)


do_calibration <- function(model, obs) {
  
  
  
}
