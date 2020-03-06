library(haven)
library(dplyr)
library(rstan)

# Data file 'daF2076e.por' is available for research at Aila,
# data service portal of the Finnish Social Science Data Archive
# See https://services.fsd.uta.fi/help?lang=en
 
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

# use only complete cases
d <- d[complete.cases(d),]

d_stan <- list(income = d$income, 
               education = d$education, 
               grade = d$grade, 
               ses = d$SES, 
               gender = d$gender, 
               itpa = d$ITPA,
               N = 250,
               M = 250,
               n = nrow(d))

d_stan$education <- as.integer(d_stan$education) - 1
d_stan$gender <- as.integer(d_stan$gender) - 1
d_stan$ses <- as.integer(d_stan$ses) - 1

# z = E(Z | X, S, G)
model_ezx <- stan_model("models/lifecourse_model_ezx.stan")
fit_ezx <- sampling(model_ezx, data = d_stan, chains = 4, cores = 4,
                iter = 26000, warmup = 1000,
                refresh = 100, init = 0,
                save_warmup = FALSE)

# z = E(Z | S, G)
model_ez <- stan_model("models/lifecourse_model_ez.stan")
fit_ez <- sampling(model_ez, data = d_stan, chains = 4, cores = 4,
                    iter = 26000, warmup = 1000,
                    refresh = 100, init = 0,
                    save_warmup = FALSE)

saveRDS(fit_ez, file = "results/lifecourse_fit_ez.rds")
saveRDS(fit_ezx, file = "results/lifecourse_fit_ezx.rds")