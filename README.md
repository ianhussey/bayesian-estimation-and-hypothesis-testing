# Bayesian estimation and hypothesis testing

## Author

Ian Hussey (ian.hussey@ugent.be)

## License

GPLv3+

## Description and purpose

Workflow for Bayesian estimation and hypothesis testing using [brms](https://cran.r-project.org/web/packages/brms/index.html). Also includes an implementation of model comparison via [Leave One Out](http://www.stat.columbia.edu/~gelman/research/unpublished/loo_stan.pdf) method. No new code, just a working workflow.

brms uses lm/lme4 syntax for model specification, making it easy to switch from frequentist (generalized) linear modelling to Bayesian.

Prior is intended to mirror that employed in JASP, i.e., a Cauchy distribution with scaling factor of .707 on population effects ("fixed effects" in frequentist language) and a scaling factor of 1.0 on group effects ("random effects" in frequentist language).

