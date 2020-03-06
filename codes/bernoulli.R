library(dplyr)
library(foreach)

run_loop <- function(nsim, n) {
  
  methods <- c("Z = 0",
               "Z = 1", 
               "P(Z)",
               "P(Z | X = x)",
               "Fully observed") %>% 
    factor(., ., ordered = TRUE)
  
  
  simulate_once <- function(n, k) {
    
    v <- rbinom(n, 1, prob = 0.5)
    u <- rbinom(n, 1, prob = 0.5)
    w <- rbinom(n, 1, prob = 0.4 * u + 0.4 * v)
    z <- rbinom(n, 1, prob = 0.4 + 0.4 * w)
    x <- rbinom(n, 1, prob = 0.4 * z + 0.4 * v)
    y <- rbinom(n, 1, prob = 0.4 * x + 0.4 * u)
    
    d <- data.frame(value = NA, 
                    y = 0:1,
                    do_x = rep(0:1, each = 2),
                    method = rep(methods, each = 4),
                    n = n, replication = k)
    
    pw1 <- mean(w)
    pw0 <- 1 - pw1
    
    for(i in 1:4) {
      Z <- 0
      pxw0 <- mean(x[z == Z & w == 0] == d$do_x[i]) * pw0
      pxw1 <- mean(x[z == Z & w == 1] == d$do_x[i]) * pw1
      
      pyw0 <- mean(y[x == d$do_x[i] & z == Z & w == 0] == d$y[i])
      pyw1 <- mean(y[x == d$do_x[i] & z == Z & w == 1] == d$y[i])
      
      d$value[i] <- (pyw0 * pxw0 + pyw1 * pxw1) / (pxw0 + pxw1)
    }
    for(i in 5:8) {
      Z <- 1
      
      pxw0 <- mean(x[z == Z & w == 0] == d$do_x[i]) * pw0
      pxw1 <- mean(x[z == Z & w == 1] == d$do_x[i]) * pw1
      
      pyw0 <- mean(y[x == d$do_x[i] & z == Z & w == 0] == d$y[i])
      pyw1 <- mean(y[x == d$do_x[i] & z == Z & w == 1] == d$y[i])
      
      d$value[i] <- (pyw0 * pxw0 + pyw1 * pxw1) / (pxw0 + pxw1)
    }
    
    for(i in 9:12) {
      mz <- mean(z)
      d$value[i] <- weighted.mean(c(d$value[i - 8], d$value[i - 4]), c(1 - mz, mz))
    }
    for(i in 13:16) {
      mz <- mean(z[x == d$do_x[i]])
      d$value[i] <- weighted.mean(c(d$value[i - 12], d$value[i - 8]), c(1 - mz, mz))
    }
    for(i in 17:20) {
      d$value[i] <- 
        mean(y[x == d$do_x[i] & u == 0] == d$y[i]) * mean(u == 0) + 
        mean(y[x == d$do_x[i] & u == 1] == d$y[i]) * mean(u == 1)
    }
    
    d
  }
  estimates <- 
    foreach(i = 1:nsim, .packages = "foreach", .combine = "rbind", .inorder = FALSE) %do% {
      simulate_once(n, i)
    }
  estimates
}


## varying n
set.seed(123)
nsim <- 1e5
b_n100 <- run_loop(nsim, 100)
b_n300 <- run_loop(nsim, 300)
b_n500 <- run_loop(nsim, 500)

sumr <- rbind(b_n100, b_n300, b_n500) %>% 
  group_by(replication, n) %>% 
  filter(!anyNA(value)) %>%
  ungroup() %>%
  group_by(method, y, do_x, n) %>% 
  summarise(mean = mean(value),
            se = sd(value) / sqrt(n()),
            lwr = mean - 2*se,
            upr = mean + 2*se,
            valid = n())
saveRDS(sumr, file = "../results/summary_bernoulli.rds")
