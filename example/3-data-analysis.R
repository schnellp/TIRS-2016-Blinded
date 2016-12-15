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
seed <- 1

# Prior parameters from "1-prior-elicitation.R"
lambda.0.hat <- 0.013 # so that the expected number of >0's
                      # roughly matches what was seen
                      # in the review paper
# (1 - dpois(0, 18*0.013))
alpha <- c(339, 339 / 100 + 4) # expected number of events in review
                               # paper given lambda.0.hat
# mean(rpois(10000, n.tot*18*0.013))
beta <- c(26064, 26064 / 100 + 167) # number of person-months in review

# Critical treatment arm AESI rate
critical <- 1.2 * 0.013

# Probability threshold
p <- 0.62

# Randomization ratio
r <- 2 / 3

### FOR EXAMPLE ONLY: SIMULATED DATA ###

# Trial design
accrual.length <- 6
max.follow.up <- 18
N <- 1500
interim.times <- c(4, 8, 12, 16)

# mock info
control.rate <- 0.011
# (1 - dpois(0, (12 * 2/3 + 18 * 1/3)*0.011)) ~= 0.143
treatment.rate <- 0.02
# (1 - dpois(0, (12 * 2/3 + 18 * 1/3)*0.02)) ~= 0.245


data <- generate.data(seed, N, r, accrual.length, max.follow.up,
                      control.rate, treatment.rate)

### END OF SIMULATED DATA ###

prob <- follow.up <- events <- post.trt <- post.ctl <-
  enrolled <- numeric(length(interim.times))

for (i in 1:length(interim.times)) {
  # Modify simulated data to reflect data accumulated by interim analysis
  data.now <- data.at(data, interim.times[i])
  
  # Sample from the joint posterior
  post <- gibbs.sample(1000, 100, alpha, beta, r,
                       data.now$t, data.now$y)
  if (i == 1) {
    plot(post[, 2] ~ post[, 1],
         main="Joint Posterior at Interim 1",
         xlab="Control SAE Rate",
         ylab="Treatment SAE Rate")
  }
  enrolled[i] <- nrow(data.now)
  post.trt[i] <- mean(post[, 2])
  post.ctl[i] <- mean(post[, 1])
  prob[i] <- mean(post[, 2] > critical)
  follow.up[i] <- sum(data.now$t)
  events[i] <- sum(data.now$y)
}

# Summaries at each interim analysis
print("Control rate posterior means:")
print(post.ctl)
print("Treatment rate posterior means:")
print(post.trt)
print("Total follow-up:")
print(follow.up)
print("Total events:")
print(events)
print("Probability of treatment rate exceeding critical value:")
print(prob)
