
# pbnorm() ------------------------------------------------------------

# integrate bivariate normal density

pbnorm <- function(lo1, up1, lo2, up2, mu1, mu2, sigma1, sigma2, rho) {
  cubature::cuhre(
    function(arg) {
      y1 <- arg[1]
      y2 <- arg[2]
      # bivariate normal density function
      1 / (2 * pi * sigma1 * sigma2 * sqrt(1 - rho^2)) * 
        exp(
          -((y1 - mu1)^2 / sigma1^2 - 
              2 * rho * (y1 - mu1) * (y2 - mu2) / (sigma1 * sigma2) + 
              (y2 - mu2)^2 / sigma2^2) / 
            (2*(1 - rho^2))
        )
    },
    lowerLimit = c(lo1, lo2), upperLimit = c(up1, up2)
  )$integral
}

# obscov() ------------------------------------------------------------

# calculates the population variance-covariance matrix of observed scores

obscov <- function(lambda, tau, theta, alpha, psi) {
  
  nitem <- length(lambda)
  ntau <- ncol(tau)
  
  # means/SEs of latent responses
  mu <- lambda*alpha
  sigma <- sqrt(lambda^2 * psi + theta)
  
  # cov matrix of latent responses
  Sigma <- matrix(nrow = nitem, ncol = nitem) 
  for (i in 1:nitem) {
    Sigma[i, ] <- lambda[i] * lambda * psi
  }
  diag(Sigma) <- sigma^2
  
  if (ntau == 1) {
    muhat <- 1 - pnorm(tau, mu, sigma)
    varhat <- muhat * (1 - muhat)
  } else {
    # probability of each response category per item (e.g., P(X1 = 1))
    p <- matrix(nrow = nitem, ncol = ntau)
    for (i in 1:nitem) {
      for (j in 1:(ntau - 1)) {
        p[i, j] <- pnorm(tau[i, j + 1], mu[i], sigma[i]) -
          pnorm(tau[i, j], mu[i], sigma[i])
      }
      p[i, ntau] <- 1 - pnorm(tau[i, ntau], mu[i], sigma[i])
    }
    
    # means/vars of observed responses
    muhat <- as.vector(p %*% (1:ntau))
    varhat <- - muhat^2
    for (i in 1:ntau) {
      varhat <- varhat + p[, i]*i^2
    }
  }
  
  # cov matrix of observed responses
  covhat <- matrix(NA, nrow = nitem, ncol = nitem)
  for (i in 1:(nitem-1)) {
    for (j in 1:(nitem-i)) {
      rho <- Sigma[i + j, i] / (sigma[i] * sigma[i + j])
      mu1 <- mu[i]
      mu2 <- mu[i + j]
      sigma1 <- sigma[i]
      sigma2 <- sigma[i + j]
      pj <- matrix(ncol = ntau, nrow = ntau)
      if (ntau == 1) {
        pj <- pbnorm(tau[i], Inf, tau[i + j], Inf, mu1, mu2, sigma1, sigma2, rho)
      } else {
        for (a in 1:(ntau-1)) {
          for (b in 1:(ntau-1)) {
            pj[a, b] <- pbnorm(tau[i, a], tau[i, a + 1],
                               tau[i + j, b], tau[i + j, b + 1],
                               mu1, mu2, sigma1, sigma2, rho)*a*b
            pj[ntau, b] <- pbnorm(tau[i, ntau], Inf,
                                  tau[i + j, b], tau[i + j, b + 1],
                                  mu1, mu2, sigma1, sigma2, rho)*ntau*b
          }
          pj[a, ntau] <- pbnorm(tau[i, a], tau[i, a + 1],
                                tau[i + j, ntau], Inf,
                                mu1, mu2, sigma1, sigma2, rho)*a*ntau
          pj[ntau, ntau] <- pbnorm(tau[i, ntau], Inf,
                                   tau[i + j, ntau], Inf,
                                   mu1, mu2, sigma1, sigma2, rho)*ntau*ntau
        }
      }
      covhat[i + j, i] <- sum(pj) - muhat[i] * muhat[i + j]
    }
  }
  covhat[upper.tri(covhat)] <- t(covhat)[upper.tri(t(covhat))]
  diag(covhat) <- varhat
  
  return(covhat)
}

# effsize() ------------------------------------------------------------

# calculates the population effect size of observed scores

effsize <- function(lambda, tau, theta1, theta2, alpha1, alpha2, psi) {
  
  nitem <- length(lambda)
  ntau <- ncol(tau)
  
  # means/SEs of latent responses
  mu1 <- lambda * alpha1
  mu2 <- lambda * alpha2
  sigma1 <- sqrt(lambda^2 * psi + theta1)
  sigma2 <- sqrt(lambda^2 * psi + theta2)
  
  # means/vars of observed responses
  if (ntau == 1) {
    muhat1 <- 1 - pnorm(tau, mu1, sigma1)
    muhat2 <- 1 - pnorm(tau, mu2, sigma2)
  } else {
    # probability of each response category per item (e.g., P(X1 = 1))
    p1 <- p2 <- matrix(nrow = nitem, ncol = ntau)
    for (i in 1:nitem) {
      for (j in 1:(ntau-1)) {
        p1[i, j] <- pnorm(tau[i, j + 1], mu1[i], sigma1[i]) -
          pnorm(tau[i, j], mu1[i], sigma1[i])
        p2[i, j] <- pnorm(tau[i, j + 1], mu2[i], sigma2[i]) -
          pnorm(tau[i, j], mu2[i], sigma2[i])
      }
      p1[i, ntau] <- 1 - pnorm(tau[i, ntau], mu1[i], sigma1[i])
      p2[i, ntau] <- 1 - pnorm(tau[i, ntau], mu2[i], sigma2[i])
    }
    muhat1 <- as.vector(p1 %*% (1:ntau))
    muhat2 <- as.vector(p2 %*% (1:ntau))
  }
  
  var1 <- sum(obscov(lambda, tau, theta1, alpha1, psi))
  diff <- mean(muhat2) - mean(muhat1)
  return(diff / sqrt(var1 / nitem ^ 2))
}