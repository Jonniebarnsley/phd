library(dplyr)
library(RobustGaSP)
library(lhs)

############################### PART 1: Setting up the Emulator ##############################

# import data
setwd("~/code/phd/Emulation/Pliocene")
data <- read.csv('pliocene_minus_control.csv')
ensemble <- na.omit(dplyr::select(
  data,
  gamma0, 
  UMV, 
  LRP,
  PDDi,
  WeertC,
  model,
  Pliocene,
  Control,
  Plio_minus_ctrl
  ))

# emulator design
inputs <- dplyr::select(
  ensemble,
  gamma0,
  UMV,
  LRP,
  PDDi,
  WeertC
  )

# Create the pairs plot
panel.corr <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y), digits=3)
  txt <- paste0("Corr: ", r)
  text(0.5, 0.5, txt, cex = 1)
}

panel.density <- function(x, pch=NULL, cex=NULL, cex.axis=NULL){
  #usr <- par('usr'); on.exit(par(usr))
  dens <- density(x)
  par(usr = c(min(x), max(x), 0, 1.5*max(dens$y)))
  lines(dens)
}

install.packages('ggplot')
library(ggplot)

panel.scat <- function(x, y){
  par(las = 0, cex.axis = 1.2)
  axis(1, at=0, labels='0')
  points(x, y, pch = 1, cex = 0.5, col = rgb(0, 0, 0, alpha=0.5))
}

nearest_order <- round(10^round(log10(mean(inputs$gamma0))))
plot(inputs$gamma0/nearest_order, inputs$UMV, xlab=substr(as.character(nearest_order), 2, 10))

options(scipen = -1)
pairs(
  inputs,
  diag.panel = panel.density,
  upper.panel = NULL,
  lower.panel = panel.scat,
  cex.labels=1,
  gap=1.5
  )

plot(density(inputs$gamma0))
log10(nearest_order)


# emulator response
output <- select(ensemble, Control)

# set a linear mean basis function - tells the emulator that we expect some kind of 
# relationship between the inputs and outputs. Will default to linear in the absence of data.
trend <- as.matrix(cbind(1, inputs))

# create model
model <- rgasp(
  design = inputs, 
  response = output,
  #nugget.est = TRUE, # nugget=TRUE accounts for factors not included in inputs (in this case, GCM)
  trend=trend,
  kernel_type = 'matern_3_2',
  )

########################## PART 2: Evaluating the Emulator ###############################

### INERT INPUTS ###

# check for inert inputs (inputs that do not evoke a statistically significant response
# in the outputs)
P <- findInertInputs(model)

### LEAVE-ONE-OUT ANALYSIS ###

# We want to test the accuracy of the emulator. We do this by removing a point from the
# inputs and asking the emulator to predict it. My "LeaveOneOut.R" script loads some
# custom functions for plotting these, including error bars and colouring according to
# whether the actual value is within the emulator uncertainty.

source("../LeaveOneOut.R")

loo <- leave_one_out(model)

pdf('./plots/leave_one_out.pdf')
par(mfrow=c(1, 1)); plot.leave_one_out(loo)
dev.off()

summary.leave_one_out(loo)

### SENSITIVITY ANALYSIS ###

# We want to test how sensitive the outputs are to each input and what relationship they
# have with sea level contribution.

# Start by setting out limits of the parameter space
LB <- c(9620, 6e17, 0, 0.008, 7600)
UB <- c(471000, 1e21, 8e-4, 0.02, 62000)
range <- UB-LB
midpoints <- (LB+UB)/2

num.testpoints = 1000
# iterate over parameters
pdf('./plots/maineffects.pdf')
par(mfrow=c(2,3)); for (i in 1:5) {
  
  # generate a sensitivity.inputs matrix which takes a range of points for our
  # parameter and the midpoint value for every other parameter on every line
  plot_param <- as.matrix(seq(LB[i], UB[i], length.out=num.testpoints))
  sensitivity.inputs <- matrix(rep(midpoints, each = num.testpoints), nrow = num.testpoints)
  sensitivity.inputs[, i] <- plot_param
  
  # define testing trend in the same way as the design trend
  sensitivity.trend <- as.matrix(cbind(1, sensitivity.inputs))
  
  # emulate the output at each point
  sensitivity.predict <- predict(
    model,
    sensitivity.inputs,
    testing_trend=sensitivity.trend
  )
  
  # plot args
  xmin <- plot_param[1]
  xmax <- plot_param[num.testpoints]
  ymin <- min(sensitivity.predict$lower95)
  ymax <- max(sensitivity.predict$upper95)
  param <- colnames(inputs)[i]
  
  # make empty plot
  plot(1, 
       1,
       type='l',
       xlim=c(xmin, xmax),
       ylim=c(ymin, ymax),
       xlab=param,
       ylab='sea level contribution (m)'
  )
  # shade 5-95% uncertainty region
  polygon(
    c(plot_param, rev(plot_param)), 
    c(sensitivity.predict$lower95, rev(sensitivity.predict$upper95)), 
    col='grey80', 
    border=F
    )
  # plot emulator mean
  lines(plot_param, sensitivity.predict$mean, type='l')
}
dev.off()

### SEA LEVEL CONTRIBUTION ###

# We want to re-sample the parameter space many times and create a histogram / probability
# distribution for sea level contribution.

# Latin Hypercube sample with 2000 samples from 5-D space
LHS <- maximinLHS(n=2000, k=5)

# re-scale the LHS to match the parameter space
for(i in 1:5) {
  LHS[,i] = LB[i] + range[i]*LHS[,i]
}

# define testing trend in the same way as design trend
testing_trend <- cbind(1, LHS)

# emulate the output at each point
model.predict <- predict(
  model, 
  LHS,
  testing_trend=testing_trend
  )

# kernel density estimate
kde <- density(model.predict$mean)
kd <- density(ensemble$Control)

# plot histogram with kernel density estimate overlay

pdf('./plots/simulatorvsemulator.pdf')
par(mfrow=c(1, 2))
hist(
  ensemble$Control, 
  breaks=seq(-26, 26, 2), 
  freq=F,
  ylim=c(0, 0.12),
  xlab='Simulator sea \n level contribution (m)',
  main=''
)
lines(kd, lwd=2, col='black')

hist(
  model.predict$mean, 
  xlim=c(-25, 25),
  ylim=c(0, 0.12),
  breaks=seq(-26, 26, 2),
  freq=F,
  xlab='Emulator sea \n level contribution (m)',
  main=''
)
lines(kde, lwd=2, col='black')
dev.off()

test <- leave_one_out_rgasp(model)
plot(ensemble$Control, test$mean)
leave_one_out_rgasp(model)
