# demonstration of bayesian estimation via brms using simulated data
# includes both estimation of posterior parameters and also BF hypothesis tests

# remaining quesitons:
# - do HDI intervals and BF always agree?
# - pp_check() ??
# - plot(marginal_effects()) ??
# - shinystan app has lots of tests I need to read up on
# - needs model building via LOO - antonio has a working implementation of this
# - might need prior robustness analysis. think more about prior generally.
# - predictions?
# - mean, median, or mode posterior?
# - how to use posterior of this as prior for next study?
# - how to implement sequential testing like this?
# - what are the units returned? 
# - report mean, median or modal estimated value?


# dependencies ------------------------------------------------------------


library(ez) 
library(schoRsch)
require(brms)
require(tidyverse)
require(parallel)

modal_value <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


# generate data -----------------------------------------------------------


# data simulated to have two main effects and an interaction effect.
# there are probably more efficient ways to do this, but it works.

# set seed for reproducability
set.seed(0)

# condition 1 data
condition <- c(rep(1, 100))
T1 <- rnorm(100, mean = 0, sd = 1)
T2 <- rnorm(100, mean = 1, sd = 1)
generated_data_1 <- data.frame(condition, T1, T2)

# condition 2 data
condition <- c(rep(2, 100))
T1 <- rnorm(100, mean = 0.4, sd = 1)
T2 <- rnorm(100, mean = 0.9, sd = 1)
generated_data_2 <- data.frame(condition, T1, T2)

# combined data
generated_data <- 
  bind_rows(generated_data_1, generated_data_2) %>%
  tibble::rownames_to_column(var = "participant") %>%
  gather(timepoint, dv, c(T1, T2)) %>%  # reshape to long format
  mutate(timepoint = as.factor(timepoint),
         condition = as.factor(condition),
         participant = as.factor(participant))


# plot --------------------------------------------------------------------


# summarise data for plotting
summary_data <- generated_data%>%
  group_by(timepoint, condition) %>%
  summarize(mean_dv = mean(dv))

ggplot(summary_data, aes(x = timepoint, y = mean_dv, group = condition, colour = condition)) +
  geom_line() +
  theme_minimal()


# check differences via freq anova ----------------------------------------


# check main and interaction effects via frequentist analysis for comparison
fit_1 <- ezANOVA(data = generated_data,
                 dv = dv,
                 within = timepoint,
                 between = condition,
                 wid = participant,
                 type = 3,
                 detailed = TRUE)

anova_out(fit_1, 
          print = TRUE, 
          sph.cor = "GG", 
          mau.p = 0.05,
          etasq = "partial", 
          dfsep = ", ")


# bayesian estimation -----------------------------------------------------


# set working directory
setwd("/Users/Ian/git/brms-testing/")

# set prior
prior <- c(set_prior("cauchy(0, 0.707)", class = "b"))  # set all beta priors, i.e., not error prior, intercept prior (nb this requires non-lmer syntax!), or random effect prior

# fit model
fit_generated <- brm(formula = dv ~ timepoint * condition,
                     data = generated_data,
                     family = gaussian(link = "identity"),
                     #family = cumulative(link = "logit"),  # ordinal alternative
                     prior = prior,
                     iter = 2000,  # default 
                     chains = 4,  # default
                     sample_prior = TRUE,  # can be used to calculate BF?
                     cores = detectCores())  # use MCMC

# save/load model to save time
save(fit_generated, file = "fit_generated")
load(file = "fit_generated")

# plot diagnostics
plot(fit_generated, ask = FALSE)

# summary of posterior parameters
summary(fit_generated)  # HDI for interaction does not contain zero, but BF test below

# sample from posterior distribution
samples <- posterior_samples(fit_generated, "b")

# there is disagreement as to whehter one should report mean, median, or modal estiamted value. compare:
mean(samples$`b_timepointT2:condition2`)
median(samples$`b_timepointT2:condition2`)
modal_value(samples$`b_timepointT2:condition2`)

# make predictions 
#predict(fit_generated)

# hypothesis test via BF
H1 <- fit_generated %>% 
  hypothesis(hypothesis = "timepointT2:condition2 = 0", alpha=.05)  # Evid.Ratio is BF01 i.e., for this hypothesis
H1
# convert BF01 to BF10:
1/H1$hypothesis$Evid.Ratio

# plot
#plot(H1)  # very wide range, probably due to use of cauchy prior? 

# plot(hypothesis) returns a ggplot object that can be modified:
p1 <- plot(H1, plot = FALSE, theme = theme_get())[[1]]  
p1 + scale_x_continuous(limits = c(-2, 2)) +
  theme_minimal() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = H1$hypothesis$`l-95% CI`) +
  geom_vline(xintercept = H1$hypothesis$`u-95% CI`)
  # attempt to shade the HDI, but not working:
  #geom_area(aes(x = ifelse(x > H1$hypothesis$`l-95% CI` & x < H1$hypothesis$`u-95% CI`, x, 0)), fill = "red") 

# shinystan for diagnostics, plotting, etc.
launch_shiny(fit_generated)

