## # Macroevolutionary analysis with microevolutionary priors
## This tutorial will guide you through the use of soon to be released bayou v2.0. The goal of bayou
## is to fit Ornstein-Uhlenbeck models to phylogenies without a priori hypotheses about regime placement. 
## bayou v2.0 adds flexibility for implementation of customized models. 
## First, let's install and load the R packages we will use for this exercise:
##+echo=FALSE
setwd("~/repos/EQG2016/bayou/R/")
##+eval=FALSE
devtools::install_github("uyedaj/bayou", ref="dev")
install.packages("treeplyr")

##+echo=TRUE
library(bayou)
library(treeplyr)

## # Preparing data for analysis with the R package treeplyr.
## We will use data from the Species 360 database of normal animal reference ranges for a variety of 
## traits (formerly known as the International Species Information System = ISIS). The data file has
## been slightly modified here, but the original data file is available at: 
## http://www2.isis.org/support/MEDARKS/Pages/Reference%20Ranges.aspx
sp360 <- readRDS("../data/species360.rds")
tbl_df(sp360)

## The tree file we are using is a mashup of several trees spliced together that includes 29,301 species of
## vertebrate. First we're going to make sure we resolve polytomies and set zero length branches to some 
## small value. 
tree <- read.tree("../data/tetrapods.tre")
tree <- multi2di(tree, random=FALSE)
tree$edge.length[tree$edge.length==0] <- .Machine$double.eps

## Matching the tree to the data is easy and fast in treeplyr, simply use the command make.treedata(). The 
## resulting list includes a pruned tree 'phy' and a pruned data table in the same order as the tree 'dat':
td <- make.treedata(tree, sp360)

## Let's choose a trait to study. In this case, let's use Weight from 1.8-2.2 years old (this isn't the best
## dataset for body size, but we'll use it for now).
tdW <- filter(td, !is.na(Weight_1.82.2yearsage) & Weight_1.82.2yearsage > 0) %>% 
  mutate(., lnMass = log(Weight_1.82.2yearsage), lnMass.SD = Weight_1.82.2yearsage.SD/Weight_1.82.2yearsage,N = Weight_1.82.2yearsage.N) %>%
  dplyr::select(., lnMass, lnMass.SD, N)

## Visualize the trait on the tree (always do this!!)
phenogram(tdW$phy, tdW[['lnMass']], spread.labels=FALSE, fsize=0.4)

## Run a bayou analysis by first setting up a prior distribution for parameters. Let's split up the data frame now
## that everything is order and make our prior distribution.
tree <- tdW$phy
dat <- tdW[[1]]
ME <- tdW[['lnMass.SD']]/sqrt(tdW[['N']])

priorOU <- make.prior(tree, 
                      dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", 
                                 dk="cdpois", dtheta="dnorm"),
                        param=list(dk=list(lambda=15, kmax=200), dsb=list(bmax=1, prob=1), 
                                   dtheta=list(mean=mean(dat), sd=1.5*sd(dat)))
                     )

## We can simulate some values from the prior distribution to test to make sure our prior function works.
startpars <- priorSim(priorOU, tree, plot=TRUE)$pars[[1]]
priorOU(startpars)

## Now that we have a prior function, let's build an MCMC object that will help us run our analyses.
set.seed(1)
mcmcOU <- bayou.makeMCMC(tree, dat, SE=ME, prior=priorOU, new.dir="../output/runs/modelOU/", outname="modelOU_r001", plot.freq=NULL)

## To run the analysis, use the 'run' function in the mcmcOU object and specify the number of generations. 
## If this amount of generations is insufficient, the MCMC chain will restart from its previous point (as long as
## the output files have not been removed.)
mcmcOU$run(10000)

## Load the chain back into R and clean up the output files, saving the chain as an RDS file. 
chainOU <- mcmcOU$load(saveRDS=TRUE)

## Set the burnin proportion for plotting with coda to check the MCMC chains
chainOU <- set.burnin(chainOU, 0.3)
summary(chainOU)
plot(chainOU)

## ## Visualization
## Visualizing a tree this large is difficult on a screen as small as ours, so let's output the file to a pdf with
## larger dimensions, so that we can read the tip labels. 
##+eval=FALSE
pdf("../output/runs/modelOU/ModelOUHeatMap.pdf", height=20, width=8)
##+eval=TRUE
#Posterior probabilities proportional to circle diameter:
plotSimmap.mcmc(chainOU, burnin=0.3, cex=0.5)
#Branch color proportional to posterior mode value for theta
plotBranchHeatMap(tree, chainOU, variable="theta", burnin=0.3, cex=0.5, lwd=2)
#Phenogram with branches with posterior probabilities > 0.3 colored and a density plot of the posterior of all thetas:
phenogram.density(tree, dat, burnin=0.3, chainOU, pp.cutoff=0.3)
##+eval=FALSE
dev.off()

##+eval=FALSE
pdf("../output/runs/modelOU/ShiftSummaries.pdf", height=8, width=8)
##+eval=TRUE
shiftsumOU <- shiftSummaries(chainOU, mcmcOU, pp.cutoff=0.5)
plotShiftSummaries(shiftsumOU)
##+eval=FALSE
dev.off()

## Or you can plot samples from the posterior distribution:
par(mfrow=c(3,3))
for(i in seq(floor(0.3*length(chainOU$gen)), length(chainOU$gen), length.out=9)){
  samplepars <- pull.pars(i, chainOU, model=mcmcOU$model.pars)
  tr <- pars2simmap(samplepars, tree)
  plotRegimes(tr$tree, pal=rainbow, lwd=2, type="phylogram", direction="upwards", show.tip.label=FALSE)
}
par(mfrow=c(1,1))

## # Alternative parameterizations
## Two alternative parameterizations of the OU model are built into bayou. First, is a parameterization where 
## priors can be placed directly on phylogenetic half-life and stationary variance, rather than alpha and sigma^2.
## For example, let's say we want to have a mildly informative prior on the phylogenetic half-life--say a log-normal
## distribution:
par.halflife <- list(meanlog=3, sdlog=2.5)
#Draw a bunch of samples from this distribution:
samp <- rlnorm(10000, par.halflife$meanlog, par.halflife$sdlog)
hist(log(samp,10), breaks=100)
abline(v=log(c(1,max(branching.times(tree))),10), col="red", lwd=2, lty=2)
## Notice that there is about equal density of prior probability on the half-life being greater than tree height (rightmost
## red line) as there is below 1 million years (leftmost red line). The exact quantiles of this distribution are:
qlnorm(c(0.025, 0.25, 0.5, 0.75, 0.975), meanlog=par.halflife$meanlog, sdlog=par.halflife$sdlog)

## Second, we'll set the prior on the stationary variance of the OU process from the Blunderbuss model, as this seems to be
## the "niche width" we expect to see on million year timescales. We will center it at the estimate from the multiple-burst
## model from Uyeda et al. 2011:
par.Vy <- list(meanlog=log(0.0958), sdlog=0.2)
hist(rlnorm(10000, par.Vy$meanlog, par.Vy$sdlog))


priorBB <- make.prior(tree, 
                     dists=list(dhalflife="dlnorm", dVy="dlnorm", 
                                dk="cdpois", dsb="dsb", dtheta="dnorm"),
                     param=list(dhalflife=par.halflife,
                                dVy=par.Vy,
                                dk=list(lambda=15, kmax=200), dsb=list(bmax=1, prob=1), 
                                dtheta=list(mean=mean(dat), sd=1.5*sd(dat))),
                     model="OUrepar"
)


## We then run everything as we did before.

startparsBB <- priorSim(priorBB, tree, plot=TRUE)$pars[[1]]
priorBB(startparsBB)

## Now that we have a prior function, let's build an MCMC object that will help us run our analyses. 
set.seed(1)
mcmcBB <- bayou.makeMCMC(tree, dat, SE=ME, model="OUrepar", prior=priorBB, new.dir="../output/runs/modelBB/", outname="modelBB_r001", plot.freq=NULL)
mcmcBB$run(10000)

chainBB <- mcmcBB$load(saveRDS=TRUE)
chainBB <- set.burnin(chainBB, 0.3)
summary(chainBB)
plot(chainBB)

##+eval=FALSE
pdf("../output/runs/modelBB/ModelBBHeatMap.pdf", height=20, width=8)
##+eval=TRUE
#Posterior probabilities proportional to circle diameter:
plotSimmap.mcmc(chainBB, burnin=0.3, cex=0.5)
#Branch color proportional to posterior mode value for theta
plotBranchHeatMap(tree, chainBB, variable="theta", burnin=0.3, cex=0.5, lwd=2)
#Phenogram with branches with posterior probabilities > 0.3 colored and a density plot of the posterior of all thetas:
phenogram.density(tree, dat, burnin=0.3, chainBB, pp.cutoff=0.3)
##+eval=FALSE
dev.off()

## ## Quantitative Genetics Model
## We can also fit a model that follows the Quantitative Genetics parameterization. You can check the specified prior
## distributions on your own, but these are informative priors that specify moderate to high heritability, realistic
## phenotypic variances, large uncertainty regarding the strength of selection and reasonable effective population 
## sizes for entire species. 
par.h2 <- list(shape1=10, shape2=10)
par.P <- list(meanlog=log(0.12), sdlog=0.2)
par.w2 <- list(meanlog=log(100), sdlog=2.5)
par.Ne <- list(meanlog=log(500000), sdlog=2.5)

## We should rescale the branch lengths to correspond roughly to generation time. However, here we will assume that for
## most of the history of birds and mammals, the generation time has been around 2 year/gen.
QGtree <- tree
QGtree$edge.length <- QGtree$edge.length/2

priorQG <- make.prior(QGtree, 
                      dists=list(dh2="dbeta", dP="dlnorm",
                                 dw2="dlnorm", dNe="dlnorm",
                                 dk="cdpois", dtheta="dnorm"),
                      param=list(dh2=par.h2,
                                 dP=par.P,
                                 dw2=par.w2,
                                 dNe=par.Ne,
                                 dk=list(lambda=30, kmax=200), dsb=list(bmax=1, prob=1), 
                                 dtheta=list(mean=mean(dat), sd=1.5*sd(dat))),
                      model="QG"
)

## Note that this model has difficulty fitting if the starting point is a poor fit. So rather than drawing from the prior distribution,
## we will start with shifts chosen by previous analyses:
endparBB <- pull.pars(length(chainBB$gen), chainBB, mcmcBB$model.pars)
startpars <- priorSim(priorQG, QGtree)$pars[[1]]
startpars <- c(startpars[1:4], endparBB[3:8])
startpars$loc <- startpars$loc/2

set.seed(100)
mcmcQG <- bayou.makeMCMC(QGtree, dat, SE=ME, model="QG", startpar=startpars, prior=priorQG, new.dir="../output/runs/modelQG/", outname="modelQG_r001", plot.freq=NULL)
mcmcQG$run(10000)

chainQG <- mcmcQG$load(saveRDS=TRUE)
chainQG <- set.burnin(chainQG, 0.3)
summary(chainQG)
plot(chainQG)

##+eval=FALSE
pdf("../output/runs/modelQG/ModelQGHeatMap.pdf", height=20, width=8)
##+eval=TRUE
#Posterior probabilities proportional to circle diameter:
plotSimmap.mcmc(chainQG, burnin=0.3, cex=0.5)
#Branch color proportional to posterior mode value for theta
plotBranchHeatMap(QGtree, chainQG, variable="theta", burnin=0.3, cex=0.5, lwd=2)
#Phenogram with branches with posterior probabilities > 0.3 colored and a density plot of the posterior of all thetas:
phenogram.density(QGtree, dat, burnin=0.3, chainQG, pp.cutoff=0.3)
##+eval=FALSE
dev.off()


## # Model Comparison
## Alternative parameterizations, shift locations, and priors can be compared using Bayes Factors. This requires estimation 
## of the marginal likelihood, which can be difficult. bayou uses stepping-stone sampling to estimate the marginal likelihoods. 
## To estimate marginal likelihoods, using the '$steppingstone' function in the mcmc object. For this exercise, we will do a much
## shorter run than is recommended. If you have multiple cores available on your machine, you can make use of these to run the
## stepping stone analysis in parallel and conduct the analysis much faster. 

library(foreach)
library(doParallel)
registerDoParallel(cores=2)

# Choose the steps from the reference distribution to the posterior distribution, this should be 20-100 steps for a full analysis, 
# but we will use only 5.Note that for each step, an MCMC chain must be run, meaning that this step could take substantially longer
# than the original analysis.

Bk <- qbeta(seq(0,1, length.out=5), 0.3,1)
ssOU <- mcmcOU$steppingstone(10000, chainOU, Bk, burnin=0.3, plot=FALSE)
ssBB <- mcmcBB$steppingstone(10000, chainBB, Bk, burnin=0.3, plot=FALSE)
ssQG <- mcmcQG$steppingstone(10000, chainQG, Bk, burnin=0.3, plot=FALSE)

mlnL <- c("OU"=ssOU$lnr, "BB"=ssBB$lnr, "QG"=ssQG$lnr)
mlnL

## # Allometric models
## What if there is a known (or unknown) relationship between the trait of interest and another predictor variable? For example,
## we may be interested in a relationship between trait known to vary with body size, but consider the possibility that the relationship
## with body size itself varies over macroevolutionary time. Here, instead of having a single optimum that changes upon a regime shift,
## it is possible to have both the slope and intercept of the relationship change at once. bayou v2.0 allows you to include these additional
## predictors and test for shifts in the scaling between a trait and its predictors. 

## Well studied allometries such as the relationships between morphological variables (e.g. brain-body size) or metabolic rate and body size
## can be studied using this method. But let's do a relationship from the species 360 database that is less explored, the evolution of body 
## temperature. 
tdTb <- mutate(td, TbK = BodyTemperature + 273.15, lnGlucose = log(GLUCOSE),
                lnMass = log(Weight_1.82.2yearsage),
                  TbK.SD = BodyTemperature.SD,
                    N = BodyTemperature.N) %>% 
                  dplyr::select(., TbK, TbK.SD, N, lnMass, lnGlucose) %>% filter(., !is.na(TbK), !is.na(lnGlucose), 
                                                                                 is.finite(lnMass), TbK<320, !is.na(lnMass))

tdTb <- mutate(tdTb, ME=ifelse(N==1 | TbK.SD==0, mean(TbK.SD), TbK.SD/sqrt(N)))
tdTb
#Note, I'm filtering out data with mass missing purely to make the dataset smaller and more tractable, you could
#run this analysis on the full dataset without any problem.

plot(tdTb[['lnGlucose']], tdTb[['TbK']])
## Let's create a separate table for the predictors, and scale them to a common scale. 
tree <- tdTb$phy
dat <- tdTb[['TbK']]
ME <- tdTb[['ME']]
pred <- mutate(tdTb$dat, lnGlu = scale(lnGlucose), lnM = scale(lnMass)) %>%  dplyr::select(., lnGlu, lnM)
lnGlu <- pred[1]
lnM <- pred[2]

## We are going to test 3 models in this analysis: Global intercepts & slopes (11), Separate intercepts & global slope (N1), and separate
## intercepts & slopes (NN).
prior.11 <- make.prior(tree, plot.prior = FALSE, 
                           dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnGlu="dnorm",
                                      dsb="fixed", dk="fixed", dtheta="dnorm"), 
                           param=list(dbeta_lnGlu=list(mean=0, sd=1),
                                     dtheta=list(mean=mean(dat), sd=1.5*sd(dat))),
                           fixed=list(k=0, sb=numeric(0), loc=numeric(0))
)

prior.N1 <- make.prior(tree, plot.prior = FALSE, 
                       dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnGlu="dnorm",
                                  dsb="dsb", dk="cdpois", dtheta="dnorm"), 
                       param=list(dbeta_lnGlu=list(mean=0, sd=1),
                                  dk=list(lambda=15, kmax=200),
                                  dsb=list(bmax=1, prob=1),
                                  dtheta=list(mean=mean(dat), sd=1.5*sd(dat)))
)


prior.NN <- make.prior(tree, plot.prior = FALSE, 
                       dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnGlu="dnorm",
                                  dsb="dsb", dk="cdpois", dtheta="dnorm"), 
                       param=list(dbeta_lnGlu=list(mean=0, sd=1),
                                  dk=list(lambda=15, kmax=200), 
                                  dsb=list(bmax=1, prob=1),
                                  dtheta=list(mean=mean(dat), sd=1.5*sd(dat)))
)

cache <- bayou:::.prepare.ou.univariate(tree, dat, SE=ME, pred=pred)
## Manually set tuning parameters:
D11 = list(alpha=1, sig2=1, beta_lnGlu=0.5, k=1, theta=0.5, slide=1)
DN1 = list(alpha=1, sig2=1, beta_lnGlu=0.5, k=1, theta=2, slide=1)
DNN = list(alpha=1, sig2=1, beta_lnGlu=1, k=c(1,1), theta=2, slide=1)

set.seed(1)
model.11 <- makeBayouModel(dat ~ lnGlu, rjpars = c(), cache=cache, prior=prior.11, D=D11)
model.N1 <- makeBayouModel(dat ~ lnGlu, rjpars = c("theta"), cache=cache, prior=prior.N1, D=DN1)
model.NN <- makeBayouModel(dat ~ lnGlu, rjpars = c("theta", "lnGlu"), cache=cache, prior=prior.NN, D=DNN)

## Make MCMC objects:
mcmc.11 <- bayou.makeMCMC(tree, dat, pred=lnGlu, SE=ME, model=model.11$model, prior=prior.11, startpar=model.11$startpar, new.dir="../output/runs/TbK_11/", outname="model11_r001", plot.freq=NULL)
mcmc.N1 <- bayou.makeMCMC(tree, dat, pred=lnGlu, SE=ME, model=model.N1$model, prior=prior.N1, startpar=model.N1$startpar, new.dir="../output/runs/TbK_N1/", outname="modelN1_r001", plot.freq=NULL)
mcmc.NN <- bayou.makeMCMC(tree, dat, pred=lnGlu, SE=ME, model=model.NN$model, prior=prior.NN, startpar=model.NN$startpar, new.dir="../output/runs/TbK_NN/", outname="modelNN_r001", plot.freq=NULL)

mcmc.11$run(10000)
mcmc.N1$run(10000)
mcmc.NN$run(10000)

chain.11 <- mcmc.11$load(saveRDS=TRUE)
chain.N1 <- mcmc.N1$load(saveRDS=TRUE)
chain.NN <- mcmc.NN$load(saveRDS=TRUE)

chain.11 <- set.burnin(chain.11, 0.3)
chain.N1 <- set.burnin(chain.N1, 0.3)
chain.NN <- set.burnin(chain.NN, 0.3)
plot(chain.N1)
summary(chain.N1)

##+eval=FALSE
pdf("../output/runs/TbK_N1/ShiftSummaries_N1.pdf", height=8, width=8)
##+eval=TRUE
shiftsumsN1 <- shiftSummaries(chain.N1, mcmc.N1, pp.cutoff=0.5)
plotShiftSummaries(shiftsumsN1, lwd=2)
##+eval=FALSE
dev.off()

##+eval=FALSE
pdf("../output/runs/TbK_NN/ShiftSummaries_NN.pdf", height=8, width=8)
##+eval=TRUE
shiftsumsNN <- shiftSummaries(chain.NN, mcmc.NN, pp.cutoff=0.5)
plotShiftSummaries(shiftsumsNN, lwd=2)
##+eval=FALSE
dev.off()

## We can do model comparison as well with allometric models. Often, these reversible jump models 
## can be hard to get to converge. It may be wise to fix the shifts with high posterior probabilities and do model
## selection on these models, rather than doing the fully reversible-jump analysis.

# Choose the steps from the reference distribution to the posterior distribution, this should be 20-100 steps for a full analysis, 
# but we will use only 5.Note that for each step, an MCMC chain must be run, meaning that this step could take substantially longer
# than the original analysis.
Bk <- qbeta(seq(0,1, length.out=5), 0.3,1)
registerDoParallel(cores=2)
set.seed(1)
ss.11 <- mcmc.11$steppingstone(10000, chain.11, Bk, burnin=0.3)
ss.N1 <- mcmc.N1$steppingstone(10000, chain.N1, Bk, burnin=0.3)
ss.NN <- mcmc.NN$steppingstone(10000, chain.NN, Bk, burnin=0.3)

ss.TbK <- list(ss.11, ss.N1, ss.NN)
sapply(ss.TbK, plot)
setNames(sapply(ss.TbK, function(x) x$lnr), c("11", "N1", "NN"))
