library(haven)
library(dplyr)
library(rstan)
rstan_options(auto_write = TRUE)
# Data file 'daF2076e.por' is available for research at Aila,
# data service portal of the Finnish Social Science Data Archive
# See https://services.fsd.uta.fi/help?lang=en
 
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

# use only complete cases
d <- d[complete.cases(d),]
round(sapply(split(d$income, d$education), mean))
round(sapply(split(d$income, d$education), median))

d_stan <- list(income = d$income, 
               education = d$education, 
               grade = d$grade, 
               ses = d$SES, 
               gender = d$gender, 
               itpa = d$ITPA,
               M = 250,
               N = 250,
               n = nrow(d))

d_stan$education <- as.integer(d_stan$education) - 1
d_stan$gender <- as.integer(d_stan$gender) - 1
d_stan$ses <- as.integer(d_stan$ses)

# z = E(Z | X, S, G)
model_x <- stan_model("models/lifecourse_withx.stan")

fit_x <- sampling(model_x, data = d_stan, chains = 4, cores = 4,
                iter = 26000, warmup = 1000,
                refresh = 100, init = 0,
                save_warmup = FALSE)

# z = E(Z | S, G)
model_nox <- stan_model("models/lifecourse_nox.stan")
fit_nox <- sampling(model_nox, data = d_stan, chains = 4, cores = 4,
                    iter = 26000, warmup = 1000,
                    refresh = 100, init = 0,
                    save_warmup = FALSE)

saveRDS(fit_x, file = "results/lifecourse_fit_withx.rds")
saveRDS(fit_nox, file = "results/lifecourse_fit_nox.rds")