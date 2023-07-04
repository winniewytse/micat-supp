# R script for the Simulation Study in the manuscript, 
# "Does Unique Factor Invariance Matter? Valid Group Mean 
# Comparisons with Ordered-Categorical Items". 

# Observed Mean Comparisons

# Setup ------------------------------------------------------------------------

library(SimDesign)

# Helper functions

source("helper_fun.R")

# Design factors

design_factor <- createDesign(
  n = c(100, 200, 500), # group size
  nninv = c(0, 1, 3), # number of noninvariant items
  dninv = c(1.25^2, 1.5^2),  # degree of noninvariance (theta2 / theta1)
  alpha2 = c(0, .2, .5),  # latent mean of group 2
  rc = c("2", "5", "7"),  # number of response categories
  subset = !(nninv == 0 & dninv == 1.5^2)
)

# To generate negatively skewed observed data
par_val <- list(
  `2` = list(
    lambda = rep(.6, 10), 
    tau = matrix(rep(-0.594, 10), nrow = 10), 
    theta = rep(.64, 10)
  ), 
  `5` = list(
    lambda = rep(.6, 10), 
    tau = matrix(c(rep(c(-1.555, -1.08, -0.553, 0.151), 10)), 
                 nrow = 10, byrow = TRUE), 
    theta = rep(.64, 10)
  ), 
  `7` = list(
    lambda = rep(.6, 10), 
    tau = matrix(c(rep(c(-1.645, -1.227, -0.915, -0.613, -0.279, 0.305), 10)), 
                 nrow = 10, byrow = TRUE),
    theta = rep(.64, 10)
  )
)

# Population observed mean difference
pop <- unique(design_factor[, c("alpha2", "rc")])
pop$pop_md <- lapply(1:nrow(pop), function(i) {
  rc <- pop[i, ]$rc
  lambda <- par_val[[rc]][["lambda"]]
  tau <- par_val[[rc]][["tau"]]
  theta1 <- theta2 <- par_val[[rc]][["theta"]]
  nitem <- length(lambda)
  alpha2 <- pop[i, ]$alpha2
  md(lambda, tau, theta1, theta2, 0, alpha2, 1)
}) |> 
  unlist()

# Fixed objects
fixed_objs <- list(
  par_val = par_val, 
  nitem = length(par_val[["2"]][["lambda"]]), 
  pop = pop
)


# Data Generation --------------------------------------------------------------

Generate <- function(condition, fixed_objects = NULL){
  rc <- condition$rc
  lambda <- fixed_objects$par_val[[rc]][["lambda"]]
  tau <- fixed_objects$par_val[[rc]][["tau"]]
  theta <- fixed_objects$par_val[[rc]][["theta"]]
  
  n <- condition$n
  nitem <- fixed_objs$nitem
  nninv <- condition$nninv
  dninv <- condition$dninv
  theta1 <- theta
  theta2 <- theta * c(rep(1, nitem - nninv), rep(dninv, nninv))
  
  alpha2 <- condition$alpha2
  
  eta1 <- rnorm(n)
  eta2 <- rnorm(n, mean = alpha2)
  
  # generate errors for latent responses
  err1 <- matrix(
    unlist(lapply(theta1, function(x) { rnorm(n, sd = sqrt(x)) })), 
    ncol = nitem
  )
  err2 <- matrix(
    unlist(lapply(theta2, function(x) { rnorm(n, sd = sqrt(x)) })), 
    ncol = nitem
  )
  
  # generate latent responses
  ystar <- rbind(eta1 %*% t(lambda) + err1, 
                 eta2 %*% t(lambda) + err2)
  
  # generate observed responses
  y <- matrix(
    unlist(lapply(1:nitem, function(i) {
      cut(ystar[, i], breaks = c(-Inf, tau[i, ], Inf), labels = FALSE)
    })), 
    ncol = nitem
  ) - 1
  
  dat <- data.frame(y, group = rep(c("g1", "g2"), each = n))
  return(dat)
}

# Data Analysis ----------------------------------------------------------------

# Observed mean comparison
Analyse.omc <- function(condition, dat, fixed_objects = NULL) {
  nitem <- fixed_objects$nitem
  n <- condition$n
  dat$mean_score <- rowSums(dat[, 1:nitem]) / nitem
  # t-test
  t <- t.test(dat$mean_score[dat$group == "g1"], 
              dat$mean_score[dat$group == "g2"])
  md <- unname(diff(t$estimate)) # mean difference
  c(md = md, se = t$stderr, pval = t$p.value)
}

# Evaluation -------------------------------------------------------------------

Summarise <- function(condition, results, fixed_objects = NULL) {
  pop <- fixed_objects$pop
  pop_md <- pop$pop_md[pop$alpha2 == condition$alpha2 & pop$rc == condition$rc]
  omc_summary <- c(
    est.md = mean(results[, "omc.md"]), 
    se.md = mean(results[, "omc.se"]), 
    bias.md = bias(results[, "omc.md"], parameter = pop_md), 
    stdbias.md = bias(results[, "omc.md"], parameter = pop_md, 
                      type = "standardized"), 
    rejrate.md = mean(results[, "omc.pval"] < .05)
  )
  omc_summary
}

# Run Simulation ---------------------------------------------------------------

res1 <- runSimulation(
  design_factor, 
  replications = 2500, 
  generate = Generate, 
  analyse = Analyse.omc, 
  summarise = Summarise, 
  fixed_objects = fixed_objs, 
  seed = rep(123, nrow(design_factor)), 
  parallel = TRUE, 
  ncores = 7, 
  save = TRUE, 
  save_seeds = TRUE,
  save_results = TRUE
)

