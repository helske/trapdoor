library(dplyr)
library(ggplot2)
library(ggthemes)
library(patchwork)
theme_set(theme_bw(base_size = 11))
size <- 0.5

# Bernoulli experiment
sumr <- readRDS("results/summary_bernoulli.rds")
levels(sumr$method) <- 
  c("z = 0", "z = 1", "P(Z)", "P(Z | X = x)", "Fully observed")
sumr %>% filter(y == 0) %>%
  mutate(n = factor(n), 
    var = paste0("y = ", y, ", x = ", do_x)) %>%
  ggplot(aes(y=mean, x = n, colour = method)) + 
  geom_hline(aes(yintercept = mean), 
    data = data.frame(
      var = paste0("y = 0, x = ", 0:1), 
      mean = c(0.8, 0.4)), linetype = "dashed", size = 0.5*size) +
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
ggsave("paper/figures/bernoulli.pdf", width = 17, height = 8, unit = "cm")  

# Gaussian experiment

sumr <- readRDS("results/summary_gaussian.rds")

sumr$method <- factor(sumr$method, levels = levels(sumr$method)[c(4:6, 1, 2, 3, 7)])
sumr$method <- forcats::fct_recode(sumr$method, 
  "z = 0" = "Z = 0",  "E(Z)" = "P(Z)", "E(Z | X = x)" = "P(Z | X = x)")
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

ggsave("paper/figures/gaussian.pdf", width = 17, height = 12, unit = "cm")  

# Gaussian Monte Carlo experiment

load("results/gaussian_mc.rds")


out_ezx <- rstan::extract(fit_ezx, pars = c("mean_analytical", "mean_montecarlo", "yrep"))
out_ez <- rstan::extract(fit_ez, pars = c("mean_analytical", "mean_montecarlo", "yrep"))

methods <- c("z = E(Z | X = 0) (Monte Carlo)",
  "z = E(Z)  (Monte Carlo)",
  "z = E(Z | X = 0)",
  "z = E(Z)")

means <- data.frame(value = c(out_ezx$mean_montecarlo, out_ez$mean_montecarlo,
  out_ezx$mean_analytical, out_ez$mean_analytical),
  method = rep(factor(methods, levels = methods, ordered = TRUE),
    each = length(out_ezx$mean_analytical)))
means %>% group_by(method) %>% 
  summarise(mean = mean(value), sd = sd(value), 
    lwr = quantile(value, 0.025), upr = quantile(value, 0.975))

samples <- data.frame(value = c(out_ezx$yrep, out_ez$yrep),
  method = rep(factor(methods[1:2], levels = methods[1:2], ordered = TRUE),
    each = length(out_ezx$mean_analytical)))

p1 <- ggplot(samples, aes(value)) +
  geom_density(aes(fill = method, linetype = method, colour = method), 
    alpha = 0.5, bw = 0.2) +
  scale_x_continuous("P(Y | do(X = 0))") +
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

ggsave("paper/figures/gaussian_mc.pdf", width = 17, height = 12, unit = "cm")  


# Life course data

library(haven)
library(dplyr)
library(rstan)

daF2076e <- read_spss("daF2076e.por") %>%
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

fit_x <- readRDS("results/lifecourse_fit_withx.rds")
fit_z <- readRDS("results/lifecourse_fit_justz.rds")

print(fit_x, pars = c("mean_ydox", "median_ydox", "mcse_mean"))
print(fit_z, pars = c("mean_ydox", "median_ydox", "mcse_mean"))

z <- c("Z ~ P(Z | X = x, S = s, G = g)", 
  "Z ~ P(Z)")

sample_z <- extract(fit_z, "sample_ydox")[[1]]
sample_x <- extract(fit_x, "sample_ydox")[[1]]
ydox <- data.frame(value = c(sample_z, sample_x), 
  do_x = rep(do_x, each = nrow(sample_z)),
  z = rep(factor(z, levels = z), 
    each = prod(dim(sample_z))))

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

ggsave("paper/figures/fsd_posterior.pdf", width = 17, height = 12, unit = "cm")  

# Non-gaussian experiment

library(simhelpers)
df <- readRDS("results/nongaussian_experiment.rds")

# analytically obtained true values based on DGP
truth <- c(19687.3150968954, 24046.1409601164, 32458.8951632273)
names(truth) <- levels(df$x)
df$truth <- rep(truth, length = nrow(df))
df %>% filter(variable=="mean") %>%
  group_by(x, method) %>%
  do(calc_absolute(., estimates = value, true_param = truth)) %>% 
  select(x, method, rmse, bias, rmse_mcse, bias_mcse)
