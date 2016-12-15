###################################################
### STEP 1: PRIOR ELICITATION                   ###
###                                             ###
### This script, once configured, produces      ###
### the prior parameters a, b, c, and d         ###
### corresponding to elicited prior information ###
### about trial arm AESI rate quantiles.        ###
###                                             ###
### Patrick Schnell, 2016                       ###
###################################################

source("0-functions.R")

###################################################
### CONFIGURATION                               ###
###                                             ###
### Please replace all NAs with decimal values. ###
###################################################

# Quantiles of the control arm AESI rate
# Ex: 2% control arm AESI rate -> 0.02
control.rate.quantiles <- c(
  NA, # lower quantile (decimal)
  NA  # upper quantile (decimal)
)

# Prior probabilities of control arm
# AESI rate being at most the above quantiles
# Ex: 90% sure control AESI rate at most 2% -> 0.90
control.rate.probabilities <- c(
  NA, # lower quantile probability (decimal)
  NA  # upper quantile probability (decimal)
)

# Quantiles of the treatment arm AESI rate
# Ex: 2% treatment arm AESI rate -> 0.02
treatment.rate.quantiles <- c(
  NA, # lower quantile (decimal)
  NA  # upper quantile (decimal)
)

# Prior probabilities of treatment arm
# AESI rate being at most the above quantiles
# Ex: 90% sure treatment AESI rate at most 2% -> 0.90
treatment.rate.probabilities <- c(
  NA, # lower quantile probability (decimal)
  NA  # upper quantile probability (decimal)
)

############################
### END OF CONFIGURATION ###
############################

# Compute prior parameters
control.prior <- find.prior.params(control.rate.quantiles,
                                   control.rate.probabilities)
treatment.prior <- find.prior.params(treatment.rate.quantiles,
                                     treatment.rate.probabilities)

# Final output
print("Control AESI rate prior parameters:")
print(paste("alpha[1] =", control.prior$par['shape']))
print(paste("beta[1] =", control.prior$par['rate']))
print("Treatment AESI rate prior parameters:")
print(paste("alpha[2] =", treatment.prior$par['shape']))
print(paste("beta[2] =", treatment.prior$par['rate']))
