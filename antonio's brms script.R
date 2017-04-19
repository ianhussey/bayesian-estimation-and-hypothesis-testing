# RT ANALYSIS
  
# Linear mixed-effects models with random intercepts and fixed slopes
# (in my data random slopes never converge).   
# The random factor is participants.   
# Compute all models of interest. In my case:   
# 1. null model
# 2. main effect of SOA
# 3. main effect of cue
# 4. both main effects (SOA + cue)
# 5. full model (SOA + cue + SOA:cue)   

num.chains <- 4 # number of MCMC chains
num.iter <- 2000 # number of MCMC iterations
num.burnin <- 500 # number of burn-in samples
num.thin <- 1 # thinning rate

# how to set priors?
# first, check on which parameters you can set the priors
# example on model with only main effect of SOA:
priors.model.SOA <- get_prior(RT~SOA+(1|Participant),data=data,family="lognormal")
# you can set priors on all intercepts of fixed effects, standard deviations of random effects, and sigma (noise) 
# use different priors on fixed effects, to see whether the choice of priors dramatically affects the posterior distribution (sensitivity analysis)

# extremely wide normal distribution
priors.Normal <- c(set_prior("normal(0,100)",class="b"), # priors for fixed effects
                   set_prior("cauchy(0,2)",class="sd")) # priors for sd of random effects

# lasso prior, a shrinkage prior that places double exponential priors on fixed effects
priors.Lasso <- c(set_prior("lasso(1,scale=1)",class="b"), # scaled because SOA and MainConds are not on the same scale
                  set_prior("cauchy(0,2)",class="sd")) # priors for sd of random effects

# Cauchy distribution, wide scaling factor
priors.wideCauchy <- c(set_prior("cauchy(0,1)",class="b"), # priors for fixed effects
                       set_prior("cauchy(0,1)",class="sd")) # priors for sd of random effects

# Cauchy distribution, medium scaling factor
priors.midCauchy <- c(set_prior("cauchy(0,.707)",class="b"), # priors for fixed effects
                      set_prior("cauchy(0,.707)",class="sd")) # priors for sd of random effects

# Cauchy distribution, narrow scaling factor
priors.narrCauchy <- c(set_prior("cauchy(0,.5)",class="b"), # priors for fixed effects
                       set_prior("cauchy(0,.5)",class="sd")) # priors for sd of random effects

priors <- list(priors.Normal,priors.Lasso,priors.wideCauchy,priors.midCauchy,priors.narrCauchy) # include prior settings in a list (for loop)
priors.names <- c("wideNormal","lasso","wideCauchy","midCauchy","narrCauchy") # set prior names (useful when saving models to file)

## example: model with main effect of SOA
for(iprior in 1:length(priors.names)) {
  model.SOA <- brm(RT~SOA+(1|Participant),
                    data=data, # data
                    family="lognormal", # posterior distribution
                    prior=priors[[iprior]], # priors
                    inits="random", # initial parameter values in the MCMC chains
                    chains=num.chains, # number of MCMC chains
                    iter=num.iter, # number of MCMC iterations
                    warmup=num.burnin, # number of burn-in samples
                    thin=num.thin, # thinning rate
                    algorithm="sampling", # type of sampling algorithm (MCMC)
                    cores=cores.beast, # number of processor cores to use when running MCMC chains in parallel
                    seed=seed.smorfia) # use RNG seed specified above
  saveRDS(model.SOA,file=paste0("model_SOA_",priors.names[iprior],".rds")) # save results
}

## model diagnostics
# plot MCMC chains
model.SOA.chains <- plot(model.SOA,ask=FALSE) # model with main effect of SOA

# plot posterior predictive checks of full model
# random effects of participant
pp_check(model.SOA,nsamples=NULL,type="stat_grouped",group="Participant")
# fixed effects of SOA
pp_check(model.SOA,nsamples=NULL,type="stat_grouped",group="SOA")
plot(marginal_effects(model.SOA),ask=FALSE)
dev.off()

# leave-one-out (LOO) cross-validation (for model selection)
compare.loo <- LOO(model.null,model.SOA,model.Cue,model.mains,model.full)
saveRDS(compare.loo,file=paste0("model_comparison_",priors.names[iprior],".rds")) # save results

# load model comparisons
compare.loo.wideCauchy <- readRDS("model_comparison_wideCauchy.rds")
compare.loo.midCauchy <- readRDS("model_comparison_midCauchy.rds")
compare.loo.narrCauchy <- readRDS("model_comparison_narrCauchy.rds")
# the best model seems to be the full model, regardless of the chosen prior

# which is the best prior?
compare.loo.wideCauchy$model.full$looic-compare.loo.midCauchy$model.full$looic
compare.loo.wideCauchy$model.full$looic-compare.loo.narrCauchy$model.full$looic
compare.loo.midCauchy$model.full$looic-compare.loo.narrCauchy$model.full$looic
# they are basically identical... midCauchy seems *slightly* better than the other two

# when your full model is the best one:
model.full.win <- readRDS("model_full_midCauchy.rds") # load null model
data$RT.est <- predict(model.full.win)[,1] # add estimated RTs to original dataframe
data.long <- melt(data,id.vars=.(Participant,Age,Sex,Handedness,Condition,Cue,SOA,MergedConds,MainConds,Resp),variable.name="datatype",value.name="RT") # convert to long format
levels(data.long$datatype) <- c("observed","estimated") # change names of factor levels
summary.data.est <- summarySEwithin(data.long,measurevar="RT",withinvars=c("SOA","MainConds","datatype")) # summarize data
summary.data.est$SOA <- factor(summary.data.est$SOA,levels=c("-217","-150","-83","-17","0","+17","+83","+150","+217")) # re-order factor levels as before (for plotting)
summary.data.est$Condition <- paste0(summary.data.est$MainConds,"_",summary.data.est$datatype) # create variable with merged names of conditions and datapoints (for plotting)

# plot observed and estimated data
# all conditions
ggplot(data=summary.data.est,aes(x=SOA,y=RT,group=Condition,color=Condition)) + # basic graph
  geom_line(aes(linetype=Condition), # separate lines per condition
            size=1.3) + # line: thickness
  scale_linetype_manual(values=rep(c("dotted","solid"),3)) + # different line types per condition
  scale_color_manual(values=rep(c("red","black","blue"),each=2)) + # different line colors per conditions
  geom_errorbar(aes(ymin=RT-ci,ymax=RT+ci), # error bars (mean +/- 95% confidence intervals)
                size=1, # error bars: thickness
                width=.2) + # error bars: width
  geom_point(aes(shape=Condition), # points at each SOA 
             fill="white", # point shape: color fill
             size=3) + # point shape: size
  scale_shape_manual(values=rep(c(21,22),3)) + # estimated -> circle; observed -> square
  scale_x_discrete("SOA", # x-axis: title
                   limits=levels(summary.data.est$SOA)) + # x-axis: SOAs (positive are on the right) 
  ylab("ms") + # y-axis: title
  scale_y_continuous(breaks=seq(300,800,100)) + # x-axis: tick marks
  geom_hline(yintercept=seq(300,800,100), # reference lines
             linetype="dotted", # line: type
             colour="#999999", # line: color
             size=.8, # line: thickness
             alpha=.5) + # line: transparency
  ggtitle("RTs") + # graph title
  theme_pander(base_size=20,pc="white",lp=c(.5,.2)) # custom theme

# ggsave(filename="RTs_obs_est_allConds.jpg",path=paste0(wd,"graphs/"),width=8,height=8,units="in",dpi=600)

#################################################################################################
# TO DO:
#       hypothesis testing (not done with BFs but 95% HDI)
#                         tempHyp <- hypothesis(model.full,"SOAM150:MainCondsHoriz1st - Intercept > 0")
#                         plot(tempHyp)
#       prediction in new dataset based on model
#       predict(toj_pilot,newdata=toj_annelies)
#################################################################################################


