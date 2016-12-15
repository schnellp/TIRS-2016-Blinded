###################################################
### STEP 0: DEFINING NECESSARY FUNCTIONS        ###
###                                             ###
### This script defines functions used in       ###
### carrying out blinded safety analyses.       ###
###                                             ###
### Patrick Schnell, 2016                       ###
###################################################

###################################################
### CONFIGURATION                               ###
###                                             ###
### None needed.                                ###
###################################################

params.from.quantiles <- function(q, p, family="gamma", start=c(0.5, 0.5)) {
  # Finds parameter values of distributions that closely match
  # the specified quantiles.
  #
  # Args:
  #   q: A vector of quantiles (length 2, strictly increasing)
  #   p: A vector of percentiles (length 2, strictly increasing)
  #   family: One of "normal", "beta", or "gamma" specifying
  #           which family of distributions to use
  #   start: starting values for parameter search
  #
  # Returns:
  #    A list containing a named vector of parameter values (par),
  #    a named vector containing the fitted quantiles (cdf),
  #    and the family used (family)
  
  # basic input validation
  stopifnot(length(p) == 2,
            length(q) == 2,
            p[1] < p[2],
            q[1] < q[2],
            0 < p[1],
            p[2] < 1)
  
  # select distribution function
  g <- switch(family,
              normal=pnorm,
              beta=pbeta,
              gamma=pgamma)
  
  # a wrapper of the cdf for optim
  f <-function(theta) {
    # if parameter values are invalid, return NA
    tryCatch({
      sum((g(q, theta[1], theta[2]) - p)^2)
    }, warning=function(w) {
      NA
    })
  }
  
  # find best parameter values
  sol <- optim(start, f)
  
  # prepare output
  out <- list()
  out$par <- sol$par
  names(out$par) <- switch(family,
                           normal=c("mean", "sd"),
                           beta=c("a", "b"),
                           gamma=c("shape", "rate"))
  
  out$cdf <- g(q, out$par[1], out$par[2])
  names(out$cdf) <- q
  
  out$family <- family
  
  out
}

generate.data <- function(seed, N, r, accrual.length, max.follow.up,
                          control.rate, treatment.rate) {
  # Generates exposure-time data.
  
  set.seed(seed)
  
  time.enrolled <- runif(N, 0, accrual.length)
  arm <- rbinom(N, 1, r)
  events <- rpois(N, ((1-arm) * control.rate + arm * treatment.rate) *
                  max.follow.up)
  patient.data <- data.frame(cbind(1:N, time.enrolled, events))
  colnames(patient.data) <- c("id", "time.enrolled", "events")
  
  event.data <- data.frame()
  
  for (i in 1:N) {
    if (events[i] > 0) {
      ts <- runif(events[i], 0, max.follow.up)
      trial.t <- time.enrolled[i] + ts
      for (j in 1:events[i]) {
        event.data <- rbind(event.data, c(i, ts[j], trial.t[j]))
      }
    }
  }
  
  colnames(event.data) <- c("patient.id", "patient.time", "trial.time")
  
  data <- list(patient.data=patient.data,
               event.data=event.data,
               follow.up=max.follow.up,
               accrual.length=accrual.length)
}

data.at <- function(data, trial.time) {
  # Produces a snapshot of the data provided as it would have
  # been available at the specified trial time.
  
  past.events <- which(data$event.data$trial.time <= trial.time)
  enrolled.patients <- which(data$patient.data$time.enrolled <= trial.time)
  N.now <- length(enrolled.patients)
  
  t <- trial.time - data$patient.data$time.enrolled[enrolled.patients]
  y <- numeric(N.now)
  for (i in 1:N.now) {
    y[i] <- sum(data$event.data$patient.id[past.events] == enrolled.patients[i])
  }
  
  return(data.frame(cbind(t, y)))
}

gibbs.sample <- function(iters, burn.in=100,
                         alpha, beta, r,
                         t, y) {
  # Samples from the joint posterior of the AESI rates
  
  N <- length(y)
  
  # inits
  lambda <- c(1, 1)
  a <- rbinom(N, 1, r)
  
  sample <- matrix(NA, ncol=4, nrow=(burn.in + iters))
  colnames(sample) <- c("control", "treatment")
  sample[1, ] <- c(lambda)
  
  for (iter in 2:(burn.in + iters)) {
    # print(iter)
    
    # lambdas
    lambda <- rgamma(2,
                     shape=c(
                       alpha[1] + sum((1 - a) * y),
                       alpha[2] + sum(a * y)
                     ),
                     rate=c(
                       beta[1] + sum((1 - a) * t),
                       beta[2] + sum(a * t)
                     ))
    
    # assignments
    cur.log.dens <- dpois(sum((1 - a) * y),
                          lambda[1] * sum((1 - a) * t),
                          log=TRUE) +
      dpois(sum(a * y),
            lambda[2] * sum(a * t),
            log=TRUE)
    
    pro.a <- rbinom(N, 1, r)
    
    pro.log.dens <- dpois(sum((1 - pro.a) * y),
                          lambda[1] * sum((1 - pro.a) * t),
                          log=TRUE) +
      dpois(sum(pro.a * y),
            lambda[2] * sum(pro.a * t),
            log=TRUE)
    
    log.ratio <- pro.log.dens - cur.log.dens
    if (log.ratio > 0 || rbinom(1, 1, exp(log.ratio))) {
      a <- pro.a
    }
    
    # record
    sample[iter, ] <- lambda
  }
  sample[burn.in + 1:iters, ]
}

simulate <- function(seed, N, r, accrual.length, max.follow.up,
                     control.rate, treatment.rate,
                     interim.times,
                     critical, alpha, beta) {
  # Simulates a trial by generating data and performing
  # interim analyses.
  
  data <- generate.data(seed, N, r, accrual.length, max.follow.up,
                        control.rate, treatment.rate)
  
  
  prob <- numeric(length(interim.times)) 
  
  for (i in 1:length(interim.times)) {
    data.now <- data.at(data, interim.times[i])
    post <- gibbs.sample(1000, 100, alpha, beta, r,
                         data.now$t, data.now$y)
    prob[i] <- mean(post[, 2] > critical)
  }
  prob
}

find.probability.threshold <- function(n.sims, N, r, accrual.length, max.follow.up,
                                       control.rate, treatment.rate,
                                       interim.times,
                                       critical, alpha, beta,
                                       signal.rate) {
  # Determines the probability threshold needed to achieve the
  # specified probability of producing a signal given
  # the specified trial and analysis parameters.
  
  maxprobs <- numeric(length(interim.times))
  for (i in 1:n.sims) {
    print(i)
    probs <- simulate(i, N, r, accrual.length, max.follow.up,
                      control.rate, treatment.rate,
                      interim.times,
                      critical, alpha, beta)
    maxprobs[i] <- max(probs)
  }
  
  as.numeric(quantile(maxprobs, 1-signal.rate, type=1))
}

power.calc <- function(n.sims, N, r, accrual.length, max.follow.up,
                       control.rate, treatment.rate,
                       interim.times,
                       critical, alpha, beta,
                       max.pp) {
  # Computes the probability of producing a signal under
  # the specified trial and analysis parameters.
  
  maxprobs <- numeric(length(interim.times))
  for (i in 1:n.sims) {
    print(i)
    probs <- simulate(i, N, r, accrual.length, max.follow.up,
                      control.rate, treatment.rate,
                      interim.times,
                      critical, alpha, beta)
    maxprobs[i] <- max(probs)
  }
  
  mean(maxprobs > max.pp)
}
