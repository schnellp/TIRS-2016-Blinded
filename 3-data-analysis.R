###################################################
### STEP 3: DATA ANALYSIS                       ###
###                                             ###
### This script, once configured, samples from  ###
### the full posterior distribution using a     ###
### Gibbs sampler and summarizes the posterior  ###
### probability that the treatment AESI rate    ###
### exceeds the specified threshold.            ###
###                                             ###
### Patrick Schnell, 2016                       ###
###################################################

###################################################
### CONFIGURATION                               ###
###                                             ###
### Please replace all NAs with decimal values. ###
###################################################

# Pseudo-random number generator seed
seed <- NA

# Gibbs sampler parameters
retained.iters <- NA
burn.in <- NA

# Prior parameters from "1-prior-elicitation.R"
# or other source (control, treatment)
alpha <- c(NA, NA)
beta <- c(NA, NA)

# Critical treatment arm AESI rate
critical <- NA

# Probability threshold from "2-simulation-study.R"
p <- NA 

# Randomization ratio
r <- NA

# Load data
# Must be a data frame with rows corresponding to enrolled patients and
# containing at least the following named columns:
#   t: the patient's follow-up time
#   y: the patient's reported number of AESI
data <- load(NA)

############################
### END OF CONFIGURATION ###
############################

set.seed(seed)

# Sample from joint posterior
post <- gibbs.sample(retained.iters, burn.in,
                     alpha, beta, r,
                     data$t, data$y)

# Summarize posterior
post.ctl <- mean(post$control)
post.trt <- mean(post$treatment)
follow.up <- sum(data$t)
events <- sum(data$y)
prob <- mean(post$treatment > critical)

print(paste("Control rate posterior means:", post.ctl))
print(paste("Treatment rate posterior means:", post.trt))
print(paste("Total follow-up:", follow.up))
print(paste("Total events:", events))
print(paste("Probability of treatment rate exceeding critical value:", prob))
print(paste("Raise safety signal?", ifelse(prob > p, "YES", "NO")))
plot(post$treatment ~ post$control,
     main="Joint Posterior at Interim 1",
     xlab="Control SAE Rate",
     ylab="Treatment SAE Rate")
