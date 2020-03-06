library(dplyr)
library(rstan)
rstan_options(auto_write = TRUE)
library(lavaan)
library(future)


model_fixed <- stan_model("models/gaussian_fixed_z.stan")
model_pz <- stan_model("models/gaussian_marginalize_z.stan")
model_pzx <- stan_model("models/gaussian_conditional_z.stan")
model_constr <- stan_model("models/gaussian_constrained.stan")
model_true <- stan_model("models/gaussian_with_u.stan")

sem_model <- "y ~ x
              x ~ z
              z ~ w
              x ~~ w
              y ~~ w"

estimate_once <- function(n) {
  

  # simulate new data
  v <- rnorm(n, 1, 1)
  u <- rnorm(n, 1, 1)
  w <- rnorm(n, 1 + u + v, 1)
  z <- rnorm(n, 1 + w, 1)
  x <- rnorm(n, 1 + z + v, 1)
  y <- rnorm(n, 1 + x + u, 0.1)
  
  do_x <- c(0, 3, 6, 9)
  fixed_z <- 0
  
  dim(fixed_z) <- 1L
  
  iter <- 2e4
  warmup <- 1e4
  
  d <- tibble::lst(n, y, x, z, w, do_x, fixed_z, 
                   Kx = length(do_x), Kz = length(fixed_z))
  fit <- sampling(model_fixed, data = d, 
                   iter = iter, warmup = warmup,
                   refresh = 0, chains = 1, cores = 1,
                   control = list(adapt_delta = 0.95, max_treedepth = 11),
                   save_warmup = FALSE, pars = "mean_ydox")
  ey_fixed <- get_posterior_mean(fit, "mean_ydox")
  
  fit <- sampling(model_pz, data = d, 
                   iter = iter, warmup = warmup,
                   refresh = 0, chains = 1, cores = 1,
                   control = list(adapt_delta = 0.95, max_treedepth = 11),
                   save_warmup = FALSE, pars = "mean_ydox")
  ey_pz <- get_posterior_mean(fit, "mean_ydox")
  
  fit <- sampling(model_pzx, data = d, 
                   iter = iter, warmup = warmup,
                   refresh = 0, chains = 1, cores = 1,
                   control = list(adapt_delta = 0.95, max_treedepth = 11),
                   save_warmup = FALSE, pars = "mean_ydox")
  ey_pzx <- get_posterior_mean(fit, "mean_ydox")
  
  fit <- sampling(model_constr, data = d, 
                   iter = iter, warmup = warmup,
                   refresh = 0, chains = 1, cores = 1,
                   control = list(adapt_delta = 0.95, max_treedepth = 11),
                   save_warmup = FALSE, pars = "mean_ydox")
  ey_constr <- get_posterior_mean(fit, "mean_ydox")
  
  d <- tibble::lst(n, y, x, u, do_x, Kx = length(do_x))
  fit <- sampling(model_true, data = d, 
                   iter = iter, warmup = warmup,
                   refresh = 0, chains = 1, cores = 1,
                   control = list(adapt_delta = 0.95, max_treedepth = 11),
                   save_warmup = FALSE, pars = "mean_ydox")
  ey_true <- get_posterior_mean(fit, "mean_ydox")
  
  fit <-  try(lavaan(sem_model, data.frame(y, x, z, w), 
                      int.ov.free = TRUE, int.lv.free = TRUE, 
                      auto.var = TRUE, 
                      auto.cov.lv.x = TRUE, 
                      auto.cov.y = TRUE,
                      std.ov = FALSE, 
                      meanstructure = TRUE,
                      estimator = "ML", optim.method = "BFGS",
                      start = "simple"), silent = TRUE)
  if(!inherits(fit, "try-error")) {
    ey_sem <- coef(fit)["y~1"] + coef(fit)["y~x"] * do_x
  } else ey_sem <- rep(NA, length(do_x))
  
  fit_z <- lm(z ~ w)
  weights <- dnorm(z, mean(z), sd(z)) / dnorm(z, predict(fit_z), sigma(fit_z))
  fit_y <- lm(y ~ x + z, weights = weights)
  
  fit_zx <- lm(z ~ x)
  pred_z <- predict(fit_zx, newdata = data.frame(x = do_x))
  ey_cwo <- predict(fit_y, newdata = data.frame(x = do_x, z = pred_z))
  
  methods <- c("Constrained",
               "SEM", 
               "CWO",
               "Z = 0", 
               "P(Z)", 
               "P(Z | X = x)", 
               "Fully observed")
  
  data.frame(
    value = c(ey_constr, ey_sem, ey_cwo, ey_fixed, ey_pz, ey_pzx, ey_true), 
    x = do_x,
    method = rep(factor(methods, levels = methods, ordered = TRUE), 
                 each = length(do_x)), 
    n = factor(n, levels = c(1,3,5) * 100))
}

set.seed(123)

# note the parallelisation

nsim <- 1000
plan(multicore, workers = 90L)
estimates <- vector("list", nsim)

for(i in 1:nsim) {
  estimates[[i]] <- future({estimate_once(100)})
}
estdf100 <-  purrr::map_dfr(lapply(estimates, value), dplyr::as_data_frame, .id = "iter")

for(i in 1:nsim) {
  estimates[[i]] <- future({estimate_once(300)})
}
estdf300 <-  purrr::map_dfr(lapply(estimates, value), dplyr::as_data_frame, .id = "iter")

for(i in 1:nsim) {
  estimates[[i]] <- future({estimate_once(500)})
}
estdf500 <-  purrr::map_dfr(lapply(estimates, value), dplyr::as_data_frame, .id = "iter")


sumr <- rbind(estdf100, estdf300, estdf500) %>% 
  group_by(method, x, n) %>% 
  summarise(mean = mean(value, na.rm = TRUE),
            se = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))),
            lwr = mean - 2*se,
            upr = mean + 2*se)
saveRDS(sumr, file = "summary_gaussian.rds")
