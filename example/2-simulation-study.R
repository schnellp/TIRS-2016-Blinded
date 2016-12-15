###################################################
### STEP 2: SIMULATION STUDY                    ###
###                                             ###
### This script, once configured, finds a       ###
### critical probability threshold for the      ###
### safety signal and performs a simulation     ###
### study to determine the operating            ###
### characteristics of the procedure.           ###
###                                             ###
### Patrick Schnell, 2016                       ###
###################################################

source("0-functions.R")

###################################################
### CONFIGURATION                               ###
###                                             ###
### Please replace all NAs with decimal values. ###
###################################################

# Prior parameters from "1-prior-elicitation.R"
alpha <- c(340, 340 / 100 + 4)
beta <- c(26064, 26064 / 100 + 167)

# Estimated AESI rates
# Ex: prior medians
set.control.rate <- 0.013
set.treatment.rate <- set.control.rate * 1.2
# set.treatment.rate <- 0.016

# Maximum allowed treatment AESI rate
critical.rate <- set.control.rate * 1.2

# Desired signal rate (probability of raising a safety signal)
# under set values of AESI rates
signal.rate <- 0.5

# Simulate at these control arm AESI rates
control.rate <- c(0.9, 1, 1.1) * set.control.rate

# Simulate at these treatment arm AESI rates
treatment.rate <- c(0.9, 1, 1.1, 1.2, 1.3) * set.control.rate

# Number of trials per parameter configuration
n.sims <- 1000

N <- 1500
accrual.length <- 6
max.follow.up <- 18

# Perform interim analyses at these times
# Must be strictly increasing
# Vector may be of any length greater than zero
interim.times <- c(4, 8, 12, 16)

# Probability of a patient being randomized to treatment arm
# Must be strictly between 0 and 1
r <- 2 / 3

############################
### END OF CONFIGURATION ###
############################

set.seed(1)

# Find safety signal probability threshold p
print("Finding probability threshold...")
p <- find.probability.threshold(n.sims, N, r, accrual.length, max.follow.up,
                                set.control.rate, set.treatment.rate,
                                interim.times,
                                critical.rate, alpha, beta,
                                signal.rate)
print(p)

# Simulate
sensitivity <- q1 <- q2 <- q3 <-
  matrix(NA, nrow=length(control.rate), ncol=length(treatment.rate))
rownames(sensitivity) <- rownames(q1) <- rownames(q2) <-
  rownames(q3) <- control.rate
colnames(sensitivity) <- colnames(q1) <- colnames(q2) <-
  colnames(q3) <- treatment.rate

for (i in 1:length(control.rate)) {
  print(paste("control:", control.rate[i]))
  for (j in 1:length(treatment.rate)) {
    print(paste("treatment:", treatment.rate[j]))
    probs <- matrix(NA, nrow=n.sims, ncol=length(interim.times))
    maxprobs <- numeric(length(interim.times))
    stopat <- numeric(length(interim.times))
    for (k in 1:n.sims) {
      probs[k, ] <- simulate(k, N, r, accrual.length, max.follow.up,
                             control.rate[i], treatment.rate[j],
                             interim.times,
                             critical.rate, alpha, beta)
      maxprobs[k] <- max(probs[k, ])
      stopat[k] <- min(c(which(probs[k, ] > p), length(interim.times)+1))
    }
    q1[i, j] <- quantile(stopat, 0.25)
    q2[i, j] <- quantile(stopat, 0.50)
    q3[i, j] <- quantile(stopat, 0.75)
    sensitivity[i, j] <- mean(maxprobs > p)
  }
}

# Print and plot results
result <- list(signal.rate=sensitivity,
               segments.completed.q1=q1,
               segments.completed.q2=q2,
               segments.completed.q3=q3)
print(result)
save(result, file="result.RData")

matplot(treatment.rate, t(sensitivity), type="l",
        xlab="Treatment AESI Rate", ylab="Probability",
        ylim=c(0, 1))
title("Probability of Producing a Signal")
abline(v=critical.rate, col="gray")
abline(h=signal.rate, col="gray")
legend("topleft", legend=rev(control.rate), lty=3:1, col=3:1, title="Control AESI Rate")
