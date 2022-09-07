# R script for the Simulation Study in the manuscript, 
# "Does Unique Factor Invariance Matter? Valid Group Mean 
# Comparisons with Ordered-Categorical Items". 

# Setup ----------------------------------------------------------------

library(tidyverse)
library(lavaan)
library(SimDesign)
library(RPushbullet)
library(here)

# Helper functions

source("helper_fun.R")

# Design factors

design_factor <- createDesign(
  n = c(100, 200, 500), # group size
  nninv = c(0, 1, 3), # number of noninvariant items
  dninv = c(1.25^2, 1.5^2),  # degree of noninvariance (theta2 / theta1)
  alpha2 = c(0, .2, .5),  # latent mean of group 2
  rc = c(2, 5),  # number of response categories
  subset = !(nninv == 0 & dninv == 1.5^2)
)

# Single group model for the simulated data

mod <- '
factor =~ X1 + X2 + X3 + X4 + X5 + X6 + X7
'

# Parameter values
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


# Population observed mean difference
pop <- cbind(unique(design_factor[, c("alpha2", "rc")]),
             pop_ez = NA)
for (i in 1:nrow(pop)) {
  if (pop[i, ]$rc == 2) {
    lambda <- par_val$rc2$lambda
    tau <- par_val$rc2$tau
    theta1 <- theta2 <- par_val$rc2$theta
  } else {
    lambda <- par_val$rc5$lambda
    tau <- par_val$rc5$tau
    theta1 <- theta2 <- par_val$rc2$theta
  }
  nitem <- length(lambda)
  
  alpha2 <- pop[i, ]$alpha2
  pop[i, "pop_ez"] <- effsize(lambda, tau, theta1, theta2, 0,
                              alpha2, 1)
  print(i)
}


# Data Generation --------------------------------------------------

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
  # bal <- condition$bal
  
  # if (isTRUE(bal)) {
  #   theta1 <- theta * c(rep(1, nitem - nninv/2), rep(dninv, nninv/2))
  #   theta2 <- theta * c(rep(dninv, nninv/2), rep(1, nitem - nninv/2))
  # } else {
  theta1 <- theta
  theta2 <- theta * c(rep(1, nitem - nninv), rep(dninv, nninv))
  # }
  
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

# dat <- datgen(design_factor[1, ])

# Data Analysis ----------------------------------------------------------------

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
  t <- t.test(dat$hatmu[1:n], dat$hatmu[(n+1):(n*2)])
  var1 <- var(dat$hatmu[1:n]) # variance of group 1
  d <- diff(t$estimate)/sqrt(var1) # cohen's d
  g <- d*(1 - 3/(4*(2*n-2)-1)) # hedge's g
  p <- t$p.value
  
  ## latent mean comparison
  # strong invariance & delta approach
  sd_fit <- cfa(model = fixed_objects$mod, 
                data = dat, 
                estimator = "WLSMV", 
                ordered = c(paste0("X", 1:7)), 
                group = "group", 
                group.equal = c("loadings", "thresholds"), 
                parameterization = "delta", 
                std.lv = TRUE)
  
  # strong invariance & theta approach
  st_fit <- cfa(model = fixed_objects$mod, 
                data = dat, 
                estimator = "WLSMV", 
                ordered = c(paste0("X", 1:7)), 
                group = "group", 
                group.equal = c("loadings", "thresholds"), 
                parameterization = "theta", 
                std.lv = TRUE)
  
  # partial strict invariance
  if (nninv == 0) {
    pstrict_fit <- cfa(model = fixed_objects$mod, 
                       data = dat, 
                       estimator = "WLSMV", 
                       ordered = c(paste0("X", 1:7)), 
                       group = "group", 
                       group.equal = c("loadings", "thresholds", "residuals"), 
                       parameterization = "theta", 
                       std.lv = TRUE)
  } else {
    if (nninv == 1) {
      g_partial <- "X7 ~~ X7"
    } else if (nninv == 3) {
      g_partial <- c("X5 ~~ X5", "X6 ~~ X6", "X7 ~~ X7")
    }
    pstrict_fit <- cfa(model = fixed_objects$mod, 
                       data = dat, 
                       estimator = "WLSMV", 
                       ordered = c(paste0("X", 1:7)), 
                       group = "group", 
                       group.equal = c("loadings", "thresholds", "residuals"), 
                       group.partial = g_partial, 
                       parameterization = "theta", 
                       std.lv = TRUE)
  }
  
  # error handling
  fit_names <- c("sd_fit", "st_fit", "pstrict_fit")
  errors <- vapply(fit_names, 
                   function(x) {
                     inherits(get(x), "try-error")
                   }, FUN.VALUE = logical(1L))
  if (any(errors)) {
    stop(paste(fit_names[errors], collapse = ", "), "did not coverge", 
         call. = FALSE)
  }
  
  # results
  sd_res <- subset(parameterestimates(sd_fit), 
                   (lhs == "factor" & op == "~1" & group == 2) | 
                     (lhs == "factor" & rhs == "factor")
  )[, c("est", "se", "pvalue")]
  st_res <- subset(parameterestimates(st_fit), 
                   (lhs == "factor" & op == "~1" & group == 2) | 
                     (lhs == "factor" & rhs == "factor")
  )[, c("est", "se", "pvalue")]
  pstrict_res <- subset(parameterestimates(pstrict_fit), 
                        (lhs == "factor" & op == "~1" & group == 2) | 
                          (lhs == "factor" & rhs == "factor")
  )[, c("est", "se", "pvalue")]
  
  results <- c(unname(d), unname(g), var1, p, 
               unlist(sd_res)[c(paste0("est", 1:3), "se3", "pvalue3")], 
               unlist(st_res)[c(paste0("est", 1:3), "se3", "pvalue3")], 
               unlist(pstrict_res)[c(paste0("est", 1:3), "se3", "pvalue3")])
  names(results) <- c("d", "g", "var1", "p", 
                      paste0(rep(c("sd_", "st_", "pstrict_"), each = 5), 
                             c(paste0(c("fv1_", "fv2_", "fm2_"), "est"), 
                               "fm2_se", "fm2_pval")))
  return(results)
}

# results <- analyze(design_factor[1, ], dat, list(mod = mod))

# Evaluation ----------------------------------------------------------------

evaluate <- function(condition, results, fixed_objects = NULL){
  # observed mean comparison
  # true effect size of observed scores
  pop_ez <- subset(pop, alpha2 == condition$alpha2 & rc == condition$rc)$pop_ez
  
  omc_summary <- c(pop_ez = pop_ez, 
                   obs_d_est = mean(results$d),
                   obs_d_bias = bias(results$d, parameter = pop_ez), 
                   obs_d_std_bias = bias(results$d, parameter = pop_ez,
                                         type = "standardized"),
                   obs_v1_est = mean(results$var1), 
                   obs_d_rej_rate = mean(results$p < .05))
  
  # latent mean comparison
  est <- lapply(results %>% select(contains("_est")), 
                function(x) {
                  c(mean(x), 
                    bias(x, parameter = condition$alpha2), 
                    bias(x, parameter = condition$alpha2, type = "standardized"))
                }) %>% unlist()
  
  se <- lapply(results %>% select(contains("_se")), mean) %>% unlist()
  
  rej_rate <- lapply(results %>% select(contains(c("_pval"))), 
                     function(x){ mean(x < .05) }) %>% unlist()
  
  lmc_summary <- c(est, se, rej_rate)
  names(lmc_summary) <- c(paste0(rep(c("sd_", "st_", "pstrict_"), each = 9), 
                                 rep(c("fv1_", "fv2_", "fm2_"), each = 3), 
                                 rep(c("est", "bias", "std_bias"))), 
                          paste0(c("sd_", "st_", "pstrict_"), 
                                 rep(c("fm2_se", "fm2_rej_rate"), each = 3))
  )
  summary <- c(omc_summary, lmc_summary)
  return(summary)
}

# Run Simulation ----------------------------------------------------------------

runSimulation(
  design_factor, 
  replications = 2500, 
  generate = datgen, 
  analyse = analyze, 
  summarise = evaluate, 
  fixed_objects = list(model = mod), 
  seed = rep(123, nrow(design_factor)), 
  packages = c("tidyverse", "lavaan"), 
  parallel = TRUE, 
  ncores = 8, 
  filename = "analysis/simulation/mean_comp/simres", 
  save_results = TRUE, 
  save_details = list(save_results_dirname = 
                        "analysis/simulation/mean_comp/sim_details"), 
  notification = "condition"
)

