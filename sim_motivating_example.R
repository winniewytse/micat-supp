# R script for the Motivating Example in the manuscript, 
# "Does Unique Factor Invariance Matter? Valid Group Mean 
# Comparisons with Ordered-Categorical Items". 

# Setup ----------------------------------------------------------------

# load libraries
library(SimDesign)
library(tidyverse)
library(lavaan)

# set design factors 
design_factor <- createDesign(
  n = c(500), # group size
  nninv = c(3), # number of noninvariant items
  dninv = c(1.5^2),  # degree of noninvariance (theta2 / theta1)
  alpha2 = c(.2),  # latent mean of group 2
  rc = c(2, 5) # number of response categories
)

# the test of measurement invariance on the data of Sharman et al. (2019)
# can be found in tutorial_mi_test.Rmd
# below are the parameter estimates in the partial scalar model 
par_val <- list(
  rc2 = list(
    lambda = c(2.6773, 2.1622, 2.2052, 1.8245, 1.4301, 1.3843, 1.3138), 
    tau = matrix(c(-2.9202, -2.6639, -2.1266, -2.097, -2.3735, -1.4055, -1.4051), 
                 nrow = 7), 
    theta = rep(1, 7)
  ), 
  rc5 = list(
    lambda = c(1.7135, 1.6416, 2.1017, 1.7575, 1.1337, 1.085, 0.9797), 
    tau = matrix(c(-3.3364, -1.8428, -0.2168, 2.2663, -3.2859, -1.9857, -0.691, 
                   0.9879, -3.5786, -1.9017, -0.3326, 1.3201, -3.3091, -1.9407, 
                   -0.4069, 1.3702, -3.4745, -1.8771, -0.7426, 0.573, -2.0902, 
                   -1.1452, -0.3404, 1.2491, -2.3554, -1.1025, -0.5331, 0.7743), 
                 nrow = 7), 
    theta = rep(1, 7)
  )
)

# CFA model for the simulated data
mod <- '
factor =~ X1 + X2 + X3 + X4 + X5 + X6 + X7
'

# Data generating function ---------------------------------------------

datgen <- function(condition, fixed_objects = list(par_val)){
  
  rc <- condition$rc
  
  if (rc == 2) {
    lambda <- par_val$rc2$lambda
    tau <- par_val$rc2$tau
    theta <- par_val$rc2$theta
  } else if (rc == 5) {
    lambda <- par_val$rc5$lambda
    tau <- par_val$rc5$tau
    theta <- par_val$rc5$theta
  }
  
  n <- condition$n
  nitem <- length(lambda)
  ntau <- rc - 1
  nninv <- condition$nninv
  dninv <- condition$dninv
  theta1 <- theta
  theta2 <- theta * c(rep(1, nitem - nninv), rep(dninv, nninv))
  
  alpha2 <- condition$alpha2
  
  eta <- c(rnorm(n), 
           rnorm(n, mean = alpha2))
  
  # generate latent responses
  ystar <- matrix(nrow = n*2, ncol = nitem)
  for(i in 1:nitem){
    ystar[1:n, i] <- eta[1:n]*lambda[i] + rnorm(n, sd = sqrt(theta1[i]))
    ystar[(n+1):(n*2), i] <- eta[(n+1):(n*2)]*lambda[i] + rnorm(n, sd = sqrt(theta2[i]))
  }
  
  # generate observed responses
  y <- matrix(nrow = n*2, ncol = nitem)
  if (ntau == 1) {
    for(i in 1:nitem) {
      y[, i] <- if_else(ystar[, i] < tau[i], 0, 1)
    }
  } else if (ntau == 4) {
    for(i in 1:nitem) {
      y[, i] <- case_when(
        ystar[, i] < tau[i, 1] ~ 0, 
        ystar[, i] < tau[i, 2] ~ 1, 
        ystar[, i] < tau[i, 3] ~ 2, 
        ystar[, i] < tau[i, 4] ~ 3, 
        TRUE ~ 4
      )
    }
  }
  
  dat <- data.frame(y, group = rep(c("g1", "g2"), each = n))
  return(dat)
}

# Analysis function ---------------------------------------------

analyze <- function(condition, dat, fixed_objects = NULL){
  rc <- condition$rc
  if (rc == 2) {
    lambda <- par_val$rc2$lambda
    tau <- par_val$rc2$tau
  } else if (rc == 5) {
    lambda <- par_val$rc5$lambda
    tau <- par_val$rc5$tau
  }
  n <- condition$n
  nitem <- length(lambda)
  nninv <- condition$nninv
  
  ## observed mean comparison
  dat$hatmu <- rowSums(dat[, 1:nitem])/nitem
  # t-test
  t <- t.test(dat$hatmu[(n+1):(n*2)], dat$hatmu[1:n])
  
  ## latent mean comparison
  # scalar invariance
  scal_fit <- cfa(model = fixed_objects$mod, 
                  data = dat, 
                  estimator = "WLSMV", 
                  ordered = TRUE, 
                  group = "group", 
                  group.equal = c("loadings", "thresholds"), 
                  parameterization = "theta", 
                  std.lv = TRUE)
  
  # partial strict invariance
  g_partial <- c("X5 ~~ X5", "X6 ~~ X6", "X7 ~~ X7")
  pstrict_fit <- cfa(model = fixed_objects$mod, 
                     data = dat, 
                     estimator = "WLSMV", 
                     ordered = TRUE, 
                     group = "group", 
                     group.equal = c("loadings", "thresholds", "residuals"), 
                     group.partial = g_partial, 
                     parameterization = "theta", 
                     std.lv = TRUE)
  
  # results
  scal_res <- subset(parameterestimates(scal_fit, standardized = TRUE), 
                     (lhs == "factor" & op == "~1" & group == 2)
  )[, c("std.all", "z", "pvalue", "ci.lower", "ci.upper")]
  pstrict_res <- subset(parameterestimates(pstrict_fit, standardized = TRUE), 
                        (lhs == "factor" & op == "~1" & group == 2)
  )[, c("std.all", "z", "pvalue", "ci.lower", "ci.upper")]
  
  results <- c(-diff(t$estimate), t$statistic, t$p.value, t$conf.int, round(t$parameter), 
               unlist(scal_res), unlist(pstrict_res))
  names(results) <- c("md", "t", "p", "ci_lo", "ci_up", "df", 
                      paste0(rep(c("sc_", "ps_"), each = 5), 
                             paste0(c("est", "z", "p", "ci_lo", "ci_up"))))
  return(results)
}

# Simulation ---------------------------------------------

set.seed(1063)
dat_rc2 <- datgen(design_factor[1, ])
res_rc2 <- analyze(design_factor[1, ], dat_rc2, fixed_objects = list(mod = mod))

set.seed(1063)
dat_rc5 <- datgen(design_factor[2, ])
res_rc5 <- analyze(design_factor[2, ], dat_rc5, fixed_objects = list(mod = mod))
round(res_rc5, 3)

