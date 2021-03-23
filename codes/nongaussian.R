library(dplyr)
library(rstan)
rstan_options(auto_write = TRUE)
library(foreach)
library(doParallel)

# Create a synthetic data set based on the the assumed graph
create_data <- function(n) {
  
  sratio_logit_lpmf <- function(x, mu, thres) {
    y <- x + 1
    ncat <- length(thres) + 1
    p <- numeric(ncat);
    q <- numeric(ncat - 1)
    k <- 1
    while (k <= min(y, ncat - 1)) {
      q[k] <- 1 - plogis(thres[k] - mu)
      p[k] <- 1 - q[k]
      if(k>1) {
        for (kk in 1:(k - 1)) p[k] <- p[k] * q[kk]
      }
      k <- k + 1;
    }
    if (y == ncat) {
      p[ncat] <- prod(q)
    }
    log(p[y])
  }
  
  sratio_logit_rng <- Vectorize(function(mu, thres) {
    
    ncat <- length(thres) + 1
    u <- runif(1)
    k <- 0
    p <- exp(sratio_logit_lpmf(k, mu, thres))
    while (p < u && k < (ncat - 1)) {
      k = k + 1
      p = p + exp(sratio_logit_lpmf(k, mu, thres))
    }
    k
  }, "mu")
  
  # confounder between W and Y
  u1 <- rnorm(n)
  # confounder between W and X
  u2 <- rnorm(n)
  # confounder between W and S
  u3 <- rnorm(n)
  
  g <- rbinom(n, 1, 0.5)
  s <- rnorm(n, 36 + 3 * u3, 5)
  
  w <- as.integer(cut(u1 + u2 + u3, c(-Inf, -1.1, 1.9, Inf)))
  
  mu_z <- plogis(-1.2 + 0.4 * g + 0.05 * s + 0.1 * (w == 2) + 0.3 * (w == 3))
  phi <- exp(2.2 + 0.2 * g)
  z <- rbeta(n, mu_z * phi, (1 - mu_z) * phi)
  
  x <- sratio_logit_rng(-0.5 * g + 0.04 * s + 13.5 * z + 2*u2, c(12.5, 14))
  
  mu_y <- exp(9.3 + 0.02 * s - 0.5 * g + 0.2 * (x == 1) + 0.5 * (x == 2) + 0.4 * u1)
  shape_y <- 10000
  y <- rgamma(n, shape_y, shape_y / mu_y)  
  
  s <- s - 36 # center for numerical efficiency 
  
  list(income = y, 
    education = as.integer(x), 
    grade = z, 
    ses = as.integer(w), 
    gender = g, 
    itpa = s,
    u1 = u1,
    u2 = u2,
    u3 = u3,
    n = n,
    N = 250,
    M = 250)
}

estimate_once <- function(n, i) {
  
  d_stan <- create_data(n)
  
  parnames <-  c("mean_ydox", "median_ydox")
  
  do_x <- structure(1:3, .Label = c("Secondary", "Lowest/lower tertiary", 
    "Higher tertiary"), class = c("ordered", "factor"))
  
  z <- c("Z ~ P(Z | X = x, S = s, G = g)", 
    "Z ~ P(Z | X = x)", 
    "Z ~ P(Z | S = s, G = g)",
    "Z ~ P(Z)", "DGP")
  z <- factor(z, levels = z)
  
  fit_x <- sampling(model_x, data = d_stan, chains = 1, cores = 1,
    iter = 6000, warmup = 1000,
    refresh = 0, init = 0, pars = parnames,
    save_warmup = FALSE)
  
  fit_onlyx <- sampling(model_onlyx, data = d_stan, chains = 1, cores = 1,
    iter = 6000, warmup = 1000, pars = parnames,
    refresh = 0, init = 0,
    save_warmup = FALSE)
  
  fit_nox <- sampling(model_nox, data = d_stan, chains = 1, cores = 1,
    iter = 6000, warmup = 1000, pars = parnames,
    refresh = 0, init = 0,
    save_warmup = FALSE)
  
  fit_z <- sampling(model_z, data = d_stan, chains = 1, cores = 1,
    iter = 6000, warmup = 1000, pars = parnames,
    refresh = 0, init = 0,
    save_warmup = FALSE)
  
  fit_dgp <- sampling(model_dgp, data = d_stan, chains = 1, cores = 1,
    iter = 6000, warmup = 1000, pars = parnames,
    refresh = 0, init = 0,
    save_warmup = FALSE)
  
  data.frame(replication = i,
    value = c(get_posterior_mean(fit_x, pars = parnames), 
      get_posterior_mean(fit_onlyx, pars = parnames), 
      get_posterior_mean(fit_nox, pars = parnames),
      get_posterior_mean(fit_z, pars = parnames),
      get_posterior_mean(fit_dgp, pars = parnames)),
    variable = rep(c("mean", "median"), each = 3),
    x = do_x,
    method = rep(z, each = 6))
}


# z ~ P(Z | X, S, G)
model_x <- stan_model("lifecourse_withx.stan")
# z ~ P(Z | X)
model_onlyx <- stan_model("lifecourse_onlyx.stan")
# z ~ P(Z | S, G)
model_nox <- stan_model("lifecourse_nox.stan")
# z ~ P(Z)
model_z <- stan_model("lifecourse_justz.stan")
# DGP
model_dgp <- stan_model("lifecourse_dgp.stan")

set.seed(123)
n <- 500
nsim <- 1000

print("start sampling")
a <- proc.time()

cl <- makeCluster(40)
registerDoParallel(cl)
estdf <- 
  foreach(i = 1:nsim,
    .packages = "rstan", 
    .combine = "rbind", .inorder = FALSE) %dopar% {
      estimate_once(n, i)
    }
stopCluster(cl)

saveRDS(estdf, file = "nongaussian_experiment.rds")

print("ready")
proc.time() - a

