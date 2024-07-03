library(dplyr)
library(RobustGaSP)
library(lhs)

############################### PART 1: Setting up the Emulator ##############################

# set working directory
setwd("~/code/phd")

# import data
data <- read.csv('data/pliocene_summary.csv')
ensemble <- na.omit(dplyr::select(
  data,
  gamma0,         # ice shelf basal melt sensitivity
  UMV,            # upper mantle viscosity
  LRP,            # precipitation lapse rate
  PDDi,           # positive degree day factor ice
  WeertC,         # Weertman friction coefficient
  model,          # forcing GCM
  Pliocene,       # Sea level contribution in the Pliocene
  Control,        # Sea level contribution in the control
  Plio_minus_ctrl # Difference between previous two
  ))

# spatial and depth-averaged ocean thermal forcing for each GCM
tf <- c("cesm"=5.63, "ccsm4uoft"=4.57, "hadcm3"=4.13, "cosmos"=3.86, "control"=2.3)

# spatially-averaged precipitation for each GCM
pr <- c('cesm'=2.1, 'ccsm4uoft'=1.96, 'hadcm3'=1.77, 'cosmos'=1.56, 'control'=0.54)

ensemble$ocean_forcing <- tf[ensemble$model]
ensemble$precip <- pr[ensemble$model]

# old attempts to incorporate thermal forcing, precipitation, and temp into the emulator
# --- ignore this ---
#tf <- c("cesm" = 3.66, "ccsm4uoft" = 2.22, "hadcm3" = 1.75, "cosmos" = 1.57, 'control'=0.72)
#tf_trimmed <- c("cesm" = 1.5, "ccsm4uoft" = 0.75, "hadcm3" = 1.17, "cosmos" = 0.76, 'control'=0.3)
#tf_all_depths <- c("cesm" = 5.63, "ccsm4uoft" = 4.57, "hadcm3" = 4.13, "cosmos" = 3.86, "control"=2.3)
#tf_all_depths_trimmed <- c("cesm" = 3.06, "ccsm4uoft" = 2.4, "hadcm3" = 2.1, "cosmos" = 1.92, 'control'=0.97)
#tf_all_depths_basins <- c("cesm" = 4.65, "ccsm4uoft" = 3.8, "hadcm3" = 2.41, "cosmos" = 2.85)
#pr <- c('cesm' = 2.1, 'ccsm4uoft' = 1.96, 'hadcm3' = 1.77, 'cosmos' = 1.56, 'control'=0.54)
#pr_trimmed <- c('cesm'=1, 'ccsm4uoft'=0.85, 'hadcm3'=0.73, 'cosmos'=0.78, 'control'=0.93)
#tas <- c('cesm' = -6.73, 'ccsm4uoft' = -10.25, 'hadcm3' = -10.72, 'cosmos' = -11.6)
#tas_trimmed <- c('cesm'=-25.32, 'ccsm4uoft'=-29.35, 'hadcm3'=-29.68, 'cosmos'=-29.22)
# --- ignore this ---


# emulator design
inputs <- dplyr::select(
  ensemble,
  gamma0,
  UMV,
  LRP,
  PDDi,
  WeertC,
  ocean_forcing,
  precip
  )

# Create the pairs plot
pairs(inputs)

# normalize inputs
normalized <- as.data.frame(scale(inputs))

# emulator response
output <- ensemble$Plio_minus_ctrl

# set a linear mean basis function - tells the emulator that we expect some kind of 
# relationship between the inputs and outputs. Will default to linear in the absence of data.
trend <- as.matrix(cbind(1, normalized))

# create model
model <- rgasp(
  design = normalized,
  response = output,
  trend=trend,
  kernel_type = 'matern_5_2',  # more smooth emulator predictions
  #kernel_type = 'matern_3_2',  # more spikey emulator predictions
  lower_bound = T,
  nugget.est = TRUE, # nugget=TRUE accounts for factors not included in inputs (e.g. full pr and tf fields)
  #zero.mean='No',
  #optimization = 'lbfgs',
  #isotropic=F,
  #alpha=rep(1.5, 5)
  )

########################## PART 2: Evaluating the Emulator ###############################

### INERT INPUTS ###

# check for inert inputs (inputs that do not evoke a statistically significant response
# in the outputs)
P <- findInertInputs(model)

### LEAVE-ONE-OUT ANALYSIS ###

# We want to test the accuracy of the emulator. We do this by removing a point from
# the inputs and asking the emulator to predict it. See "LeaveOneOut.R" for custom
# functions to plot these, including error bars and colouring according to whether
# the actual value is within the emulator uncertainty.

source("analysis/emulation/LeaveOneOut.R"); loo <- leave_one_out(model)

# plot the simulator outputs against emulator predictions
#pdf('./plots/pdf/leave_one_out.pdf')
par(mfrow=c(1, 1)); plot(loo)
#dev.off()

# Normalized Euclidean Distance and RMSE are used to quantitatively compare emulators
# Ideally, 95% of simulations will lie within the emulator uncertainty (Pass)
summary(loo)


### SENSITIVITY ANALYSIS ###

# We want to test how sensitive the outputs are to each input and what relationship they
# have with sea level contribution. See main_effects.R for the code to produce this plot.

source("analysis/emulation/main_effects.R")

#pdf('./plots/pdf/main_effects_plio_minus_ctrl.pdf', height=2000, width=2000, res=250)
main_effects(model, normalized)
#dev.off()


######################## Part 3: Sea level contribution ########################

# We want to re-sample the parameter space many times and create a histogram / probability
# distribution for sea level contribution.

# Latin Hypercube sample with 1000 samples from 7-D space
LHS <- maximinLHS(n=1000, k=7)

# define testing trend in the same way as design trend
testing_trend <- as.matrix(cbind(1, LHS))

# emulate the output at each point
model.predict <- predict(
  model, 
  LHS,
  testing_trend=testing_trend
  )

# kernel density estimate
kde <- density(model.predict$mean)
kd <- density(ensemble$Plio_minus_ctrl)

# plot histogram with kernel density estimate overlay
#pdf('./plots/pdf/simulatorvsemulator.pdf')
par(mfrow=c(1, 2))
hist(
  ensemble$Plio_minus_ctrl, 
  xlim=c(-10, 15),
  ylim=c(0, 0.3),
  breaks=seq(-10, 15, 0.5),
  freq=F,
  xlab='Simulator sea \n level contribution (m)',
  main=''
)
lines(kd, lwd=2, col='black')

hist(
  model.predict$mean, 
  xlim=c(-10, 15),
  ylim=c(0, 0.3),
  breaks=seq(0, 20, 0.5),
  freq=F,
  xlab='Emulator sea \n level contribution (m)',
  main=''
)
lines(kde, lwd=2, col='black')
#dev.off()

# okay, these two plots obviously look way different, probably because the LHS above is 
# sampling evenly from parameter space, but gamma0, thermal forcing, and precipitation are 
# not evenly distributed in the ensemble. gamma0 in particular is skewed towards lower values
# in the ensemble (see pairs plot), and has a strong positive relationship with sea level
# contribution (see main effects plot). When the LHS samples evenly, it picks up a lot more 
# higher values of gamma0 and so is yielding a much higher sea level distribution curve 
# (5-15m instead of -5-10m). I'm sure this could be fixed with some scaling, but that's
# work for another day.
