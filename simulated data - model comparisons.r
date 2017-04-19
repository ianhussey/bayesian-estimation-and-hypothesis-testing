# title: "Demonstration of bayesian estimation model comparisons using simulated data"
# author: "Ian Hussey (ian.hussey@ugent.be)"


# dependencies ------------------------------------------------------------


library(tidyverse)
library(ez) 
library(schoRsch)
library(brms)
library(parallel)


# simulate data -----------------------------------------------------------


# data simulated to have two main effects and an interaction effect.
# there are probably more efficient ways to do this, but it works.

# set seed for reproducability
set.seed(0)

# condition 1 data
condition         <- c(rep(1, 100))
T1                <- rnorm(100, mean = 0, sd = 1)
T2                <- rnorm(100, mean = 1, sd = 1)
generated_data_1  <- data.frame(condition, T1, T2)

# condition 2 data
condition         <- c(rep(2, 100))
T1                <- rnorm(100, mean = 0.4, sd = 1)
T2                <- rnorm(100, mean = 0.9, sd = 1)
generated_data_2  <- data.frame(condition, T1, T2)

# combined data
generated_data <- 
  bind_rows(generated_data_1, generated_data_2) %>%
  tibble::rownames_to_column(var = "participant") %>%
  gather(timepoint, dv, c(T1, T2)) %>%  # reshape to long format
  mutate(timepoint = as.factor(timepoint),
         condition = as.factor(condition),
         participant = as.factor(participant))


# bayesian estimation -----------------------------------------------------


# set working directory
setwd("/Users/Ian/git/brms-testing/")

# set prior
prior <- c(set_prior("cauchy(0, 0.707)", class = "b"))  # set all beta priors, i.e., not error prior, intercept prior (nb this requires non-lmer syntax!), or random effect prior


# null model --------------------------------------------------------------


fit_null <- brm(formula = dv ~ 0 + intercept,
                data = generated_data,
                family = gaussian(link = "identity"),
                prior = prior,
                iter = 2000,  
                chains = 4,  
                sample_prior = TRUE,  # to calculate BF
                cores = detectCores())  

# save/load model to save time
save(fit_null, file = "fit_null.rdata")
load(file = "fit_null.rdata")

summary(fit_null)


# condition model --------------------------------------------------------------


fit_condition <- brm(formula = dv ~ condition,
                     data = generated_data,
                     family = gaussian(link = "identity"),
                     prior = prior,
                     iter = 2000, 
                     chains = 4,
                     sample_prior = TRUE,  # to calculate BF
                     cores = detectCores())

# save/load model to save time
save(fit_condition, file = "fit_condition.rdata")
load(file = "fit_condition.rdata")

summary(fit_condition)


# timepoint model --------------------------------------------------------------


fit_timepoint <- brm(formula = dv ~ timepoint,
                     data = generated_data,
                     family = gaussian(link = "identity"),
                     prior = prior,
                     iter = 2000, 
                     chains = 4,  
                     sample_prior = TRUE,  # to calculate BF
                     cores = detectCores())

# save/load model to save time
save(timepoint, file = "timepoint.rdata")
load(file = "timepoint.rdata")

summary(fit_timepoint)


# interaction model --------------------------------------------------------------


fit_interaction <- brm(formula = dv ~ condition * timepoint,
                       data = generated_data,
                       family = gaussian(link = "identity"),
                       prior = prior,
                       iter = 2000, 
                       chains = 4,
                       sample_prior = TRUE,  # to calculate BF
                       cores = detectCores()) 

# save/load model to save time
save(interaction, file = "interaction.rdata")
load(file = "interaction.rdata")

summary(fit_interaction)


# model comparison --------------------------------------------------------


# Leave-One-Out (LOO) cross-validation (for model selection)

# WAIC (and indeed AIC) are alternatives to (Pareto-Smoothed Importance Sampling) LOO.
# Vehtari, Gelman & Gabry (2016) show that WAIC is asymptotically equal to (PSIS) LOO, 
# but the latter is more robust with weak priors or influential observations.

model_comparison <- LOO(fit_null, fit_timepoint, fit_condition, fit_interaction)

save(model_comparison, file = "model_comparison.rdata")
load(file = "model_comparison.rdata")

model_comparison

