# bayesian estimation and hypothesis testing of derivation experiment data
# author: ian hussey (ian.hussey@ugent.be)
# license: GPLv3+

# dependencies ------------------------------------------------------------


require(brms)
require(tidyverse)
require(parallel)
library(beepr)


# get data ----------------------------------------------------------------


# study 3 data - more reasonble bayes factor
setwd("/Users/Ian/Dropbox/Work/Projects/Alphabet soup/Hussey & Hughes - Derivation study/OSF - Transitive relations and implicit attitudes/Experiment 2/")

# acquire data
data_df <- 
  read.csv("Processed data/processed data - long format IAT data.csv") %>%
  mutate(rate = 1000/rt,
         block = ifelse(block == 1, "IAT block 1", "IAT block 2"),
         block = as.factor(block),
         condition = as.factor(ifelse(condition == 1, "Condition 1", "Condition 2")),
         participant = as.factor(participant)) %>%
  schoRsch::outlier(dv = "rt", 
                    res.name = "rt_outlier",
                    todo = "na", 
                    upper.z = 2.5, 
                    lower.z = -2.5) %>%
  schoRsch::outlier(dv = "rate", 
                    res.name = "rate_outlier",
                    todo = "na", 
                    upper.z = 2.5, 
                    lower.z = -2.5)

# set wd
setwd("/Users/Ian/git/Bayesian estimation and hypothesis testing/")


# distribution ------------------------------------------------------------


# needs thought here. ask antonio for his dist fit script

# remove NAs
valid_data <- data_df %>%
  filter(!is.na(rate))

# generate normal data with the same mean and sd as the data 
simulated_norm_data <- 
  rnorm(nrow(valid_data), 
        mean = mean(valid_data$rate), 
        sd = sd(valid_data$rate)) %>%
  as.data.frame() 
# change col name
colnames(simulated_norm_data) <- c("simulated_rate") 

# plot the generated and real data for fit comparison
ggplot() +
  geom_density(data = simulated_norm_data,
               aes(x = simulated_rate), 
               alpha = 0.1, colour = "red") +
  geom_density(data = valid_data,
               aes(x = rate), 
               alpha = 0.1, colour = "blue") +
  theme_minimal()


# prior -------------------------------------------------------------------


# set priors
possible_priors <- get_prior(formula = rate ~ block * condition + (1 | participant),
                             data = data_df, 
                             family = gaussian(link = "identity"))
possible_priors

# set prior to same as defaults used in JASP, i.e., cauchy with r=-.707 on "fixed" effects and r=1 on "random"
prior <- c(set_prior("cauchy(0, 0.707)", class = "b"),  # "fixed" effect
           set_prior("cauchy(0, 1.0)", class = "sd"))    # "random" effect


# model comparisons -------------------------------------------------------


# # null
# fit_null <- brm(formula = rate ~ (1 | participant),
#                 data = data_df,
#                 family = gaussian(link = "identity"),
#                 prior = prior,
#                 iter = 2000,  # default 
#                 warmup = 500, 
#                 chains = 4,  # default
#                 #cov_ranef = "participant",  # for setting random effects?
#                 sample_prior = TRUE,
#                 cores = detectCores())
# 
# # block
# fit_block <- brm(formula = rate ~ block + (1 | participant),
#                  data = data_df,
#                  family = gaussian(link = "identity"),
#                  prior = prior,
#                  iter = 2000,  # default 
#                  warmup = 500, 
#                  chains = 4,  # default
#                  #cov_ranef = "participant",  # for setting random effects?
#                  sample_prior = TRUE,
#                  cores = detectCores())
# 
# # both
# fit_both <- brm(formula = rate ~ condition + block + (1 | participant),
#                      data = data_df,
#                      family = gaussian(link = "identity"),
#                      prior = prior,
#                      iter = 2000,  # default 
#                      warmup = 500, 
#                      chains = 4,  # default
#                      #cov_ranef = "participant",  # for setting random effects?
#                      sample_prior = TRUE,
#                      cores = detectCores())
# 
# # condition
# fit_condition <- brm(formula = rate ~ condition + (1 | participant),
#                      data = data_df,
#                      family = gaussian(link = "identity"),
#                      prior = prior,
#                      iter = 2000,  # default 
#                      warmup = 500, 
#                      chains = 4,  # default
#                      #cov_ranef = "participant",  # for setting random effects?
#                      sample_prior = TRUE,
#                      cores = detectCores())

# interaction - complete pooling 
fit_full_complete <- brm(formula = rate ~ block * condition + (1 | participant),
                         data = data_df,
                         family = gaussian(link = "identity"),
                         prior = prior,
                         iter = 2000,  # default 
                         warmup = 500, 
                         chains = 4,  # default
                         #cov_ranef = "participant",  # for setting random effects?
                         sample_prior = TRUE,
                         cores = detectCores())

save(fit_full_complete, file = "deriv2_fit_full_complete.rdata")


# interaction - partial pooling
fit_full_partial <- brm(formula = rate ~ block * condition + (1 + block | participant),
                        data = data_df,
                        family = gaussian(link = "identity"),
                        prior = prior,
                        iter = 2000,  # default 
                        warmup = 500, 
                        chains = 4,  # default
                        #cov_ranef = "participant",  # for setting random effects?
                        sample_prior = TRUE,
                        cores = detectCores())

save(fit_full_partial, file = "deriv2_fit_full_partial.rdata")


# interaction - partial pooling sep slopes
fit_full_partial_sepslopes <- brm(formula = rate ~ block * condition + (1 + block || participant),
                                  data = data_df,
                                  family = gaussian(link = "identity"),
                                  prior = prior,
                                  iter = 2000,  # default 
                                  warmup = 500, 
                                  chains = 4,  # default
                                  #cov_ranef = "participant",  # for setting random effects?
                                  sample_prior = TRUE,
                                  cores = detectCores())

save(fit_full_partial_sepslopes, file = "deriv2_fit_full_partial_sepslopes.rdata")


# model comparison
partial_full_model_comparison <- LOO(fit_full_complete, fit_full_partial, fit_full_partial_sepslopes)
save(partial_full_model_comparison, file = "deriv2_partial_full_model_comparison.rdata")
partial_full_model_comparison

beep(sound = 8)


# results -----------------------------------------------------------------


# plot diagnostics
plot(fit_full_partial_sepslopes, ask = FALSE)

# summary of posterior parameters
summary(fit_full_partial_sepslopes)  # HDI for interaction does not contain zero, but BF test below

plot(marginal_effects(fit_full_partial), ask = FALSE)

# hypothesis test via BF
H1 <- fit_full_partial_sepslopes %>% 
  hypothesis(hypothesis = "blockIATblock2:conditionCondition2 = 0", alpha=.05)  # Evid.Ratio is BF01 i.e., for this hypothesis
H1
# convert BF01 to BF10:
1/H1$hypothesis$Evid.Ratio

# plot prior vs posterior
p1 <- plot(H1, plot = FALSE, theme = theme_get())[[1]]  
p1 + scale_x_continuous(limits = c(-.5, .5)) +
  theme_minimal() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = H1$hypothesis$`l-95% CI`) +
  geom_vline(xintercept = H1$hypothesis$`u-95% CI`)


