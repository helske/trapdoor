library(dplyr)
library(ggplot2)
library(ggthemes)
library(patchwork)
theme_set(theme_bw(base_size = 11))
size <- 0.5

# Bernoulli experiment
sumr <- readRDS("results/summary_bernoulli.rds")

sumr %>% 
  mutate(n = factor(n), 
         var = paste0("y = ", y, ", x = ", do_x)) %>%
  ggplot(aes(y=mean, x = n, colour = method)) + 
  geom_hline(aes(yintercept = mean), 
             data = data.frame(
               var = paste0("y = ", c(0,0,1,1), ", x = ", 0:1), 
               mean = c(0.8, 0.4, 0.2, 0.6)), linetype = "dashed", size = 0.5*size) +
  geom_linerange(aes(ymin = lwr, ymax = upr), 
                 position = position_dodge(0.75), size = size) + 
  geom_errorbar(aes(ymin = mean, ymax = mean), width = 0.5, size = size,
                position = position_dodge(0.75)) + 
  scale_colour_few() + 
  scale_y_continuous("P(Y = y | do(X = x))", minor_breaks = NULL) +
  scale_x_discrete("Sample size") +
  facet_wrap(~var, scales = "free_y", ncol = 2) +
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.5 * size))
ggsave("paper/figures/bernoulli.pdf", width = 16, height = 12, unit = "cm")  

# Gaussian experiment

sumr <- readRDS("results/summary_gaussian.rds")

sumr$method <- factor(sumr$method, levels = levels(sumr$method)[c(4:6, 1, 2, 3, 7)])
sumr$method <- forcats::fct_recode(sumr$method, "E(Z)" = "P(Z)", "E(Z | X = x)" = "P(Z | X = x)")
sumr %>% 
  ggplot(aes(y = mean, x = n, colour = method)) + 
  geom_hline(aes(yintercept = mean), 
             data = data.frame(
               x = c(0, 3, 6, 9), 
               mean = c(2, 5, 8, 11)), linetype = "dashed", size = 0.5 * size) +
  geom_linerange(aes(ymin = lwr, ymax = upr), 
                 position = position_dodge(0.75), size = size) + 
  geom_errorbar(aes(ymin = mean, ymax = mean), width = 0.5, size = size,
                position = position_dodge(0.75)) + 
  scale_colour_manual(values = few_pal()(8)[c(1, 3, 4, 6,7,8, 5)]) + 
  #scale_shape_manual(values = c(5, rev(few_shape_pal()(5)))) +
  scale_y_continuous("E(Y | do(X = x))", minor_breaks = NULL) +
  scale_x_discrete("Sample size") +
  facet_wrap(~ x, scales = "free_y", ncol = 2,
             labeller = purrr::partial(label_both, sep = " = ")) +
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.5 * size))

ggsave("paper/figures/gaussian.pdf", width = 16, height = 20, unit = "cm")  

# Gaussian Monte Carlo experiment

load("results/gaussian_mc.rds")


out_pzx <- rstan::extract(fit_ezx, pars = c("mean_analytical", "mean_montecarlo", "yrep"))
out_pz <- rstan::extract(fit_ez, pars = c("mean_analytical", "mean_montecarlo", "yrep"))

methods <- c("z ~ E(Z | X = 0) (Monte Carlo)",
             "z ~ E(Z)  (Monte Carlo)",
             "z ~ E(Z | X = 0)",
             "z ~ E(Z)")

means <- data.frame(value = c(out_pzx$mean_montecarlo, out_pz$mean_montecarlo,
                              out_pzx$mean_analytical, out_pz$mean_analytical),
                    method = rep(factor(methods, levels = methods, ordered = TRUE),
                                 each = length(out_pzx$mean_analytical)))

samples <- data.frame(value = c(out_pzx$yrep, out_pz$yrep),
                      method = rep(factor(methods[1:2], levels = methods[1:2], ordered = TRUE),
                                   each = length(out_pzx$mean_analytical)))

p1 <- ggplot(samples, aes(value)) +
  geom_density(aes(fill = method, linetype = method, colour = method), 
               alpha = 0.5, bw = 0.2) +
  scale_x_continuous("Y | do(X = 0)") +
  scale_y_continuous("Posterior density", expand = expansion(mult=c(0,0.03))) +
  geom_vline(xintercept = 2, linetype = "dashed", size = 0.5 * size) +
  scale_fill_few() + 
  scale_linetype_manual(values = c("solid", "solid", "dashed", "dotted")) +
  scale_colour_manual(values = c("black", NA, "black", "black")) +
  coord_cartesian(xlim = c(-1.5, 6)) +
  theme(legend.position = "none",#c(0.2,0.88), 
        legend.title = element_blank(), 
        panel.grid.major.x = element_line(size = 0.5 * size),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank())

p2 <- ggplot(means, aes(value)) +
  geom_density(aes(fill = method, linetype = method, colour = method), 
               alpha = 0.5, bw = 0.2) +
  scale_x_continuous("E(Y | do(X = 0))") +
  scale_y_continuous("Posterior density", expand = expansion(mult=c(0,0.03))) +
  geom_vline(xintercept = 2, linetype = "dashed", size = 0.5 * size) +
  scale_fill_few() + 
  scale_linetype_manual(values = c("solid", "solid", "dashed", "dotted")) +
  scale_colour_manual(values = c("black", NA, "black", "black")) +
  coord_cartesian(xlim = c(-1.5, 6)) +
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        panel.grid.major.x = element_line(size = 0.5 * size),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank())

p1 + p2 + plot_layout(ncol = 1)

ggsave("paper/figures/gaussian_mc.pdf", width = 16, height = 16, unit = "cm")  


# Life course data

library(haven)
library(dplyr)
library(rstan)

daF2076e <- read_spss("codes/daF2076e.por") %>%
  select(
    SES = T14,
    grade_avg = K82,
    ITPA = L17,
    gender = T3,
    education_level_2002 = BV2_1,
    income = BV4_1,
  ) %>% mutate(gender = factor(gender, labels = c("Male", "Female")),
               education = recode_factor(factor(education_level_2002),
                                         "2" = "Secondary",
                                         "3" = "Secondary",
                                         "4" = "Secondary",
                                         "5" = "Lowest/lower tertiary",
                                         "6" = "Lowest/lower tertiary",
                                         "7" = "Higher tertiary",
                                         "8" = "Higher tertiary", .ordered = TRUE),
               grade = (grade_avg-4) / 6,
               SES = recode_factor(as.character(SES), 
                                   "3" = "Low",
                                   "2" = "Middle",
                                   "1" = "High",
                                   .ordered = TRUE))

d <- daF2076e %>% select(income, education, grade, SES, ITPA, gender)
d <- d[complete.cases(d),]

do_x <- factor(levels(d$education),
               levels = levels(d$education),
               ordered = TRUE)

fit_ez <- readRDS("results/lifecourse_fit_ez.rds")
fit_ezx <- readRDS("results/lifecourse_fit_ezx.rds")

z <- c("z ~ E(Z | X = x, S = s, G = g)", 
       "z ~ E(Z | S = s, G = g)")
mean_ez <- extract(fit_ez, "eydox")[[1]]
mean_ezx <- extract(fit_ezx, "eydox")[[1]]

eydox <- data.frame(value = c(mean_ezx, mean_ez), 
                    do_x = rep(do_x, each = nrow(mean_ezx)),
                    z = rep(factor(z, levels = z), 
                            each = prod(dim(mean_ezx))))

sample_ez <- extract(fit_ez, "ydox_sample")[[1]]
sample_ezx <- extract(fit_ezx, "ydox_sample")[[1]]
ydox <- data.frame(value = c(sample_ezx, sample_ez), 
                   do_x = rep(do_x, each = nrow(sample_ezx)),
                   z = rep(factor(z, levels = z), 
                           each = prod(dim(sample_ezx))))

xrange <- c(17000, 36000)
eydox %>% 
  ggplot(aes(x = value, fill = do_x, linetype = do_x)) + 
  geom_density(alpha = 0.5, bw = 300) + 
  scale_fill_few() + 
  coord_cartesian(xlim = xrange, ylim = c(0, 5.5e-4)) +
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        panel.grid.major.x = element_line(size = 0.5 * size),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.title = element_text(size = 16)) +
  scale_y_continuous(name = "Posterior density", 
                     expand = expansion(mult=c(0,0.03))) +
  scale_x_continuous(name = "Income") + 
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  facet_wrap(~ z, ncol = 1)

ggsave("paper/figures/fsd_posterior_mean.pdf", width = 16, height = 16, unit = "cm")  

xrange <- c(0, 90000)
ydox %>% 
  ggplot(aes(x = value, fill = do_x, linetype = do_x)) + 
  geom_density(alpha = 0.5, bw = 3000) + 
  scale_fill_few() + 
  coord_cartesian(xlim = xrange, ylim = c(0, 4.2e-5)) +
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        panel.grid.major.x = element_line(size = 0.5 * size),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.title = element_text(size = 16)) +
  scale_y_continuous(name = "Posterior density", 
                     expand = expansion(mult=c(0,0.03))) +
  scale_x_continuous(name = "Income", expand = expansion(mult=c(0,0))) + 
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  facet_wrap(~ z, ncol = 1)

ggsave("paper/figures/fsd_posterior.pdf", width = 16, height = 16, unit = "cm")  
