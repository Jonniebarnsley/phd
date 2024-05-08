library(dplyr)
library(RobustGaSP)
library(lhs)

############################### PART 1: Setting up the Emulator ##############################

# import data
setwd("~/code/phd")
data <- read.csv('data/pliocene_summary.csv')
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

tf <- c("cesm" = 3.66, "ccsm4uoft" = 2.22, "hadcm3" = 1.75, "cosmos" = 1.57, 'control'=0.72)
tf_trimmed <- c("cesm" = 1.5, "ccsm4uoft" = 0.75, "hadcm3" = 1.17, "cosmos" = 0.76, 'control'=0.3)
tf_all_depths <- c("cesm" = 5.63, "ccsm4uoft" = 4.57, "hadcm3" = 4.13, "cosmos" = 3.86, "control"=2.3)
tf_all_depths_trimmed <- c("cesm" = 3.06, "ccsm4uoft" = 2.4, "hadcm3" = 2.1, "cosmos" = 1.92, 'control'=0.97)
tf_all_depths_basins <- c("cesm" = 4.65, "ccsm4uoft" = 3.8, "hadcm3" = 2.41, "cosmos" = 2.85)

pr <- c('cesm' = 2.1, 'ccsm4uoft' = 1.96, 'hadcm3' = 1.77, 'cosmos' = 1.56, 'control'=0.54)
pr_trimmed <- c('cesm'=1, 'ccsm4uoft'=0.85, 'hadcm3'=0.73, 'cosmos'=0.78, 'control'=0.93)

tas <- c('cesm' = -6.73, 'ccsm4uoft' = -10.25, 'hadcm3' = -10.72, 'cosmos' = -11.6)
tas_trimmed <- c('cesm'=-25.32, 'ccsm4uoft'=-29.35, 'hadcm3'=-29.68, 'cosmos'=-29.22)

ensemble$ocean_forcing <- tf_all_depths[ensemble$model]
ensemble$precip <- pr[ensemble$model]
ensemble$temp <- tas_trimmed[ensemble$model]

# emulator design
inputs <- dplyr::select(
  ensemble,
  gamma0,
  UMV,
  LRP,
  PDDi,
  WeertC,
  ocean_forcing,
  precip,
  #temp,
  #interaction_g0_of2,
  #interaction_pr_lrp,
)

#inputs <- dplyr::select(
#  data,
#  gamma0,
#  UMV,
#  LRP,
#  PDDi,
#  WeertC,
#  ocean_forcing,
#  precip
#)

# Create the pairs plot
pairs(inputs)

# normalize inputs
normalized <- as.data.frame(scale(inputs))

# emulator response
output <- as.data.frame(scale(ensemble$Pliocene))


# set a linear mean basis function - tells the emulator that we expect some kind of 
# relationship between the inputs and outputs. Will default to linear in the absence of data.
trend <- as.matrix(cbind(1, normalized))

# create model for Pliocene outputs
pliocene <- rgasp(
  design = normalized,
  response = ensemble$Plio_minus_ctrl,
  trend=trend,
  kernel_type = 'matern_5_2',
  lower_bound = T,
  nugget.est = T, # nugget=TRUE accounts for factors not included in inputs (in this case, GCM)
  #zero.mean='No',
  #optimization = 'lbfgs',
  #isotropic=F,
  #alpha=rep(2, 5)
)

########################## PART 2: Evaluating the Emulator ###############################

### INERT INPUTS ###

# check for inert inputs (inputs that do not evoke a statistically significant response
# in the outputs)
P <- findInertInputs(pliocene)

### LEAVE-ONE-OUT ANALYSIS ###

# We want to test the accuracy of the emulator. We do this by removing a point from the
# inputs and asking the emulator to predict it. My "LeaveOneOut.R" script loads some
# custom functions for plotting these, including error bars and colouring according to
# whether the actual value is within the emulator uncertainty.

source("analysis/emulation/LeaveOneOut.R"); loo <- leave_one_out(pliocene)

# plot the simulator outputs against emulator predictions
#pdf('./plots/pdf/loo_control.pdf')
png('./plots/png/loo_plio_minus_ctrl.png', height=2000, width=2000, res=250)
par(mfrow=c(1, 1)); plot(loo)
dev.off()

# Normalized Euclidean Distance and RMSE can be used to quantitatively compare emulators
summary(loo)

######################### SENSITIVITY ANALYSIS ##########################

# We want to test how sensitive the outputs are to each input and what relationship they
# have with sea level contribution.

source("analysis/emulation/main_effects.R")

#png('./plots/png/main_effects_plio_minus_ctrl.png', height=2000, width=2000, res=250)
main_effects(pliocene, normalized)
#dev.off()

######################### CONTROL #############################

inputs <- dplyr::select(
  ensemble,
  gamma0,
  UMV,
  LRP,
  PDDi,
  WeertC,
)

normalized <- as.data.frame(scale(inputs))
trend <- as.matrix(cbind(1, normalized))

control <- rgasp(
  design = normalized,
  response = ensemble$Control,
  trend=trend,
  kernel_type = 'matern_5_2',
  lower_bound = F,
  nugget.est = T, # nugget=TRUE accounts for factors not included in inputs (in this case, GCM)
  #zero.mean='No',
  #optimization = 'lbfgs',
  #isotropic=F,
  #alpha=rep(1.7, 5)
)
P <- findInertInputs(control)


loo <- leave_one_out(control)
par(mfrow=c(1, 1)); plot(loo)
summary(loo)

main_effects(control, normalized)

plio_minus_ctrl <- rgasp(
  design = normalized,
  response = ensemble$Plio_minus_ctrl,
  trend=trend,
  kernel_type = 'matern_5_2',
  lower_bound = T,
  nugget.est = TRUE, # nugget=TRUE accounts for factors not included in inputs (in this case, GCM)
  zero.mean='No',
  #optimization = 'lbfgs',
  #isotropic=F,
  alpha=rep(2, 5)
)

P <- findInertInputs(plio_minus_ctrl)

loo <- leave_one_out(plio_minus_ctrl)
par(mfrow=c(1, 1)); plot(loo)
summary(loo)

main_effects(plio_minus_ctrl, inputs)

######################## SEA LEVEL CONTRIBUTION #########################

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
pliocene.predict <- predict(
  pliocene, 
  LHS,
  testing_trend=testing_trend
)

# kernel density estimate
kde <- density(pliocene.predict$mean)
kd <- density(ensemble$Control)

# plot histogram with kernel density estimate overlay

#pdf('./plots/simulatorvsemulator.pdf')
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
#dev.off()

test <- leave_one_out_rgasp(model)
plot(ensemble$Control, test$mean)
leave_one_out_rgasp(model)
