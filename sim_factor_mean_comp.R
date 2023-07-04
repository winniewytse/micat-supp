# R script for the Simulation Study in the manuscript, 
# "Does Unique Factor Invariance Matter? Valid Group Mean 
# Comparisons with Ordered-Categorical Items". 

# Factor Mean Comparisons

# Setup ------------------------------------------------------------------------

library(SimDesign)
library(lavaan)

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

# Fixed objects
fixed_objs <- list(
  par_val = par_val, 
  nitem = length(par_val[["2"]][["lambda"]]), 
  mod = '
factor =~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10
'
)

# Data Generation --------------------------------------------------

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

ExtractOut <- function(fit) {
  out_names <- paste0(c("fv_", "fm_"), rep(c("est", "se", "pval"), each = 2))
  out <- subset(
    parameterestimates(fit), 
    lhs == "factor" & op %in% c("~1", "~~") & group == 2
  )[, c("est", "se", "pvalue")] |> 
    as.vector() |> 
    unlist()
  names(out) <- out_names
  out
}

# Factor mean comparison (scalar invariance with delta approach)
Analyse.sd <- function(condition, dat, fixed_objects = NULL) {
  nitem <- fixed_objects$nitem
  
  zero_cat1 <- unlist(apply(dat[dat$group == "g1", ], 2, 
                            function(x) { table(x) == 0 }))
  zero_cat2 <- unlist(apply(dat[dat$group == "g2", ], 2, 
                            function(x) { table(x) == 0 }))
  if (any(zero_cat1) | any(zero_cat2)) {
    stop(paste0("One or more items contain 0 in a response category, ", 
                "which will lead to model nonconvergence"))
  }
  
  sd_fit <- cfa(model = fixed_objects$mod, 
                data = dat, 
                estimator = "WLSMV", 
                ordered = c(paste0("X", 1:nitem)), 
                group = "group", 
                group.equal = c("loadings", "thresholds"), 
                parameterization = "delta", 
                std.lv = TRUE)
  out <- ExtractOut(sd_fit)
  if (anyNA(out)) {
    stop("Scalar model (delta) did not obtain SE")
  }
  out
}

# Factor mean comparison (scalar invariance with theta approach)
Analyse.st <- function(condition, dat, fixed_objects = NULL) {
  nitem <- fixed_objects$nitem
  
  zero_cat1 <- unlist(apply(dat[dat$group == "g1", ], 2, 
                            function(x) { table(x) == 0 }))
  zero_cat2 <- unlist(apply(dat[dat$group == "g2", ], 2, 
                            function(x) { table(x) == 0 }))
  if (any(zero_cat1) | any(zero_cat2)) {
    stop(paste0("One or more items contain 0 in a response category, ", 
                "which will lead to model nonconvergence"))
  }
  
  st_fit <- cfa(model = fixed_objects$mod, 
                data = dat, 
                estimator = "WLSMV", 
                ordered = c(paste0("X", 1:nitem)), 
                group = "group", 
                group.equal = c("loadings", "thresholds"), 
                parameterization = "theta", 
                std.lv = TRUE)
  out <- ExtractOut(st_fit)
  if (anyNA(out)) {
    stop("Scalar model (theta) did not obtain SE")
  }
  out
}

# Factor mean comparison (strict/partial strict invariance)
Analyse.str <- function(condition, dat, fixed_objects = NULL) {
  nitem <- fixed_objects$nitem
  
  zero_cat1 <- unlist(apply(dat[dat$group == "g1", ], 2, 
                            function(x) { table(x) == 0 }))
  zero_cat2 <- unlist(apply(dat[dat$group == "g2", ], 2, 
                            function(x) { table(x) == 0 }))
  if (any(zero_cat1) | any(zero_cat2)) {
    stop(paste0("One or more items contain 0 in a response category, ", 
                "which will lead to model nonconvergence"))
  }
  
  nninv <- condition$nninv
  nninv_items <- paste0("X", seq(nitem + 1 - nninv, nitem))
  if (nninv == 0) {
    g_partial <- ""
  } else {
    g_partial <- paste(nninv_items, "~~", nninv_items)
  }
  
  pstr_fit <- cfa(model = fixed_objects$mod, 
                  data = dat, 
                  estimator = "WLSMV", 
                  ordered = c(paste0("X", 1:nitem)), 
                  group = "group", 
                  group.equal = c("loadings", "thresholds", "residuals"), 
                  group.partial = g_partial, 
                  parameterization = "theta", 
                  std.lv = TRUE)
  out <- ExtractOut(pstr_fit)
  if (anyNA(out)) {
    stop("Strict/Partial strict model did not obtain SE")
  }
  out
}

# Evaluation ----------------------------------------------------------------

Summarise <- function(condition, results, fixed_objects = NULL) {
  models <- c("sd", "st", "str")
  fmc_summary <- lapply(
    results[, paste0(rep(models, 4), 
                     c(".fv_est", ".fm_est", ".fv_se", ".fm_se"))], 
    mean
  ) |> unlist()
  names(fmc_summary) <- paste0(
    "mean.", 
    paste0(rep(models, 4), c(".fv_est", ".fm_est", ".fv_se", ".fm_se"))
  )
  fm_est <- results[, paste0(models, ".fm_est")]
  alpha2 <- condition$alpha2
  fm_summary <- c(
    bias = bias(fm_est, alpha2), 
    stdbias = bias(fm_est, alpha2, type = "standardized")
  )
  fm_rejrate <- lapply(
    results[, paste0(models, ".fm_pval")], 
    function(x){ mean(x < .05) }
  ) |> unlist()
  names(fm_rejrate) <- paste0("rejrate.", models)
  c(fmc_summary, fm_summary, fm_rejrate)
}

# Run Simulation ----------------------------------------------------------------

res <- runSimulation(
  design_factor, 
  replications = 2500, 
  generate = Generate, 
  analyse = list(sd = Analyse.sd, 
                 st = Analyse.st, 
                 str = Analyse.str), 
  summarise = Summarise, 
  fixed_objects = fixed_objs, 
  seed = rep(123, nrow(design_factor)), 
  packages = "lavaan", 
  parallel = TRUE, 
  ncores = 7, 
  filename = "simres_sass_fmc_neg", 
  save = TRUE, 
  save_seeds = TRUE,
  save_results = TRUE
)
