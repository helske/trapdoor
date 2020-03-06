library(dplyr)
library(rstan)
rstan_options(auto_write = TRUE)

model_ezx <- stan_model("models/gaussian_montecarlo_ezx.stan")
model_ez <- stan_model("models/gaussian_montecarlo_ez.stan")

n <- 100
set.seed(1)
v <- 1 + rnorm(n)
u <- 1 + rnorm(n)
w <- 1 + u + v + rnorm(n)
z <- 1 + w + rnorm(n)
x <- 1 + z + v + 0.1 * rnorm(n)
y <- 1 + x + u + 0.1 * rnorm(n)

do_x <- 0
d <- tibble::lst(n, y, w, z, x, do_x, N = 500)

fit_ezx <- sampling(model_ezx, data = d,
                iter = 26000, warmup = 1000,
                chains = 4, cores = 4,
                refresh = 100, save_warmup = FALSE)


fit_ez <- sampling(model_ez, data = d,
                    iter = 26000, warmup = 1000,
                    chains = 4, cores = 4,
                    refresh = 100, save_warmup = FALSE)

save(fit_ezx, fit_ez, file = "results/gaussian_mc.rds")


