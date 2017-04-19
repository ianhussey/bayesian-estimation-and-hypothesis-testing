# title: "Demonstration of bayesian estimation and hypothesis testing using simulated data, via both Bayes Factors and Highest Density Intervals"
# author: "Ian Hussey (ian.hussey@ugent.be)"

# remaining questions:
# - do HDI intervals and BF always agree?
# - pp_check() needs some thought
# - plot(marginal_effects()) needs some thought
# - predictions needs some thought
# - shinystan shiny app has lots of checks that aren't used here
# - report mean or median posterior?
# - how to use posterior of one model as prior for a follow up study?
# - how to implement a sequential testing workflow here?


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
  rownames_to_column(var = "participant") %>%
  gather(timepoint, dv, c(T1, T2)) %>%  # reshape to long format
  mutate(timepoint = as.factor(timepoint),
         condition = as.factor(condition),
         participant = as.factor(participant))


# plot --------------------------------------------------------------------


# rehape and summarize data for plotting
summary_data <- generated_data%>%
  group_by(timepoint, condition) %>%
  summarize(mean_dv = mean(dv))

ggplot(summary_data, aes(x = timepoint, y = mean_dv, group = condition, colour = condition)) +
  geom_line() +
  theme_minimal()


# check differences via freq anova ----------------------------------------


# check main and interaction effects via frequentist analysis for the sake of familiarity
fit_1 <- 
  ezANOVA(data = generated_data,
          dv = dv,
          within = timepoint,
          between = condition,
          wid = participant,
          type = 3,
          detailed = TRUE) %>%
  anova_out(print = FALSE)

fit_1[3]


# bayesian estimation -----------------------------------------------------


# set working directory
setwd("/Users/Ian/git/brms-testing/")

# set prior
prior <- c(set_prior("cauchy(0, 0.707)", class = "b"))  # set all beta priors, i.e., not error prior, intercept prior (nb this requires non-lmer syntax!), or random effect prior

# fit model
fit_generated <- brm(formula = dv ~ timepoint * condition,
                     data = generated_data,
                     family = gaussian(link = "identity"),
                     prior = prior,
                     iter = 2000,
                     chains = 4,
                     sample_prior = TRUE,  # to calculate BF
                     cores = detectCores())

# save/load model to save time
save(fit_generated, file = "fit_generated")
load(file = "fit_generated")

# plot diagnostics
plot(fit_generated, ask = FALSE)

# hypothesis testing via HDI - summary of posterior parameters
# assess whether Highest Density Interval (HDI) contain zero? Also refered to as the (Bayesian) Credibility Interval (CI).
summary(fit_generated)

# sample from posterior distribution
samples <- posterior_samples(fit_generated, "b")

# there is disagreement as to whether one should report mean or median estimated value. compare:
mean(samples$`b_timepointT2:condition2`)
median(samples$`b_timepointT2:condition2`)

# make predictions 
#predict(fit_generated)

# hypothesis test via BF
H1 <- fit_generated %>% 
  hypothesis(hypothesis = "timepointT2:condition2 = 0", alpha=.05)  # Evid.Ratio is BF01 i.e., for this hypothesis
H1
# convert BF01 to BF10 (because in this case 01 is null and 10 is alternate):
1/H1$hypothesis$Evid.Ratio

# plot
#plot(H1)  # very wide range, probably due to use of cauchy prior. adjust the range manually:

# plot(hypothesis) returns a ggplot object that can be modified:
p1 <- plot(H1, plot = FALSE, theme = theme_get())[[1]]  
p1 + scale_x_continuous(limits = c(-2, 2)) +
  theme_minimal() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = H1$hypothesis$`l-95% CI`) +
  geom_vline(xintercept = H1$hypothesis$`u-95% CI`)
  # attempt to shade the HDI, but not working:
  #geom_area(aes(x = ifelse(x > H1$hypothesis$`l-95% CI` & x < H1$hypothesis$`u-95% CI`, x, 0)), fill = "red") 

# shinystan for diagnostics, plotting, etc. many more options here.
launch_shiny(fit_generated)

