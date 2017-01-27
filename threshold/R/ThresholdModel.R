## # Threshold models

## # PHYLIP & RPhylip
## Follow the directions to install phylip and Rphylip. 
#+echo=FALSE
setwd("~/repos/threshold/R/")
#+eval=FALSE
devtools::install_github("liamrevell/Rphylip/Rphylip")
library(Rphylip)
library(treeplyr)
library(geiger)
library(phytools)

## Let's start by simulating data under multivariate Brownian Motion with a known covariance structure.
# Simulate the tree
ntaxa <- 100
set.seed(1)
tree <- sim.bdtree(b=1, stop="taxa", n=ntaxa)
tree$edge.length <- tree$edge.length/max(branching.times(tree))

# Set up the covariance matrix for changes in the underlying liabilities
R <- matrix(c(  1,   0.7, -0.9,
                0.7,     1, -0.4,
                -0.9,   -0.4,   1), ncol=3, nrow=3, byrow=TRUE)
colnames(R) <- rownames(R) <- c("A", "B", "C")
# Simulate the data using 'sim.char()'. 
dat <- sim.char(tree, par=R, nsim=1, model="BM")[,,1]
colnames(dat) <- c("A", "B", "C")
pairs(dat)
## Just for kicks, realize that we don't really need a fancy function to simulate multivariate Brownian motion on a
## tree, we can just use matrix algebra. 
library(MASS)
root <- c(rep(0, ntaxa), rep(0, ntaxa), rep(0, ntaxa))
bigR <- kronecker(R, vcv.phylo(tree))
dat2 <- mvrnorm(1, root, bigR)
dat2 <- matrix(dat2, ncol=3)
pairs(dat2)

## Match the tree and to the liabilities and create discretized variables. 
td <- make.treedata(tree, dat)
X <- as.matrix(td[c("A", "B", "C")])
rownames(X) <- td$phy$tip.label

## A pretty way to plot the liabilities:
fancyTree(td$phy, type="scattergram", X=X)

td <- mutate(td, a = as.numeric(A > median(A)), b = as.numeric(B > median(B)), c = as.numeric(C > median(C)))
tree <- td$phy
traits <- as.data.frame(td$dat)
rownames(traits) <- td$phy$tip.label

# Plot out the discrete traits
par(mfrow=c(1,1))
plot(tree, show.tip.label=FALSE)
tiplabels(pch=22, bg=traits[['a']])
tiplabels(pch=22, bg=traits[['b']], adj=c(0.53,0.5))
tiplabels(pch=22, bg=traits[['c']], adj=c(0.56,0.5))

## # Rphylip
## Now let's try to estimate parameters for the threshold model using Rphylip & threshML. Note that if you have
## not put threshML in your system path, or in your working directory, you will need to specify where it is.
mod <- Rthreshml(tree, traits[,c('a', 'b', 'c')], types=c("discrete", "discrete", "discrete"), nchain=2, ngen=10000, 
                 lrtest=TRUE, path = "~/repos/EQG2016/threshold/")

mod$Covariance_matrix
R

## We can also include continuous traits that may be correlated to the liabilities (or the liabilities themselves!).
mod2 <- Rthreshml(tree, traits[,c('a', 'A', 'C')], types=c("discrete", "continuous", "continuous"), nchain=2, ngen=10000, 
                  lrtest=TRUE, path = "~/repos/EQG2016/threshold/")

## Note: there is in error in the labeling of the covariance matrix, see if you can't figure out what it should be.
mod2$Covariance_matrix
R

## # phytools & threshBayes
## Now we'll use phytools threshBayes function. I will rely heavily on Liam Revell's tutorial for this, so 
## all credit goes to him. His tutorial can be found here: http://www.phytools.org/eqg/Exercise_6.2/
library(coda)
ngen <- 20000
X <- as.matrix(traits[,c('A', 'C')])
mod.AC <- threshBayes(tree, X, types=c("cont", "cont"), ngen=ngen)
mcmc.AC <- mcmc(mod.AC$par)
plot(mcmc.AC)

X <- as.matrix(traits[,c('a', 'C')])
mod.aC <- threshBayes(tree, X, types=c("disc", "cont"), ngen=ngen)
mcmc.aC <- mcmc(mod.aC$par)
plot(mcmc.aC)

X <- as.matrix(traits[,c('a', 'c')])
mod.ac <- threshBayes(tree, X, types=c("disc", "disc"), ngen=ngen)
mcmc.ac <- mcmc(mod.ac$par)
plot(mcmc.ac)

par(mfrow=c(1,2))
postburnin <- floor(0.5*nrow(mcmc.AC)):nrow(mcmc.AC)
plot(density(mcmc.AC[postburnin,'r']), main="Posterior for correlation", xlab= "Correlation", xlim=c(-1,1))
lines(density(mcmc.aC[postburnin,'r']), col="blue")
lines(density(mcmc.ac[postburnin,'r']), col="green")

plot(0,0, type="n", main="Trace", xlab= "Generation", ylab="Correlation", ylim=c(-1,1), xlim=c(0, length(postburnin)))
lines(mcmc.AC[postburnin,'r'])
lines(mcmc.aC[postburnin,'r'], col="blue")
lines(mcmc.ac[postburnin,'r'], col="green")

## # MCMCglmmRAM
library(MCMCglmmRAM)
traits$animal <- factor(tree$tip.label)
prior1 <- list(R=list(V=diag(2)*1e-15, fix=1), G=list(G1=list(V=diag(2), nu=0.002)))
mglmm.AC <- MCMCglmm(cbind(C, A)~trait-1, random=~us(trait):animal, rcov=~us(trait):units,
                     pedigree=tree, reduced=TRUE, data=traits, prior=prior1, 
                     pr=TRUE, pl=TRUE, family=c("gaussian", "gaussian"), thin=1)
summary(mglmm.AC)

prior2 <- list(R=list(V=diag(2)*1e-15, fix=1), G=list(G1=list(V=diag(2), nu=0.002, fix=2)))
mglmm.aC <- MCMCglmm(cbind(C, a)~trait-1, random=~us(trait):animal, rcov=~us(trait):units,
                  pedigree=tree, reduced=TRUE, data=traits, prior=prior2, 
                    pr=TRUE, pl=TRUE, family=c("gaussian", "threshold"), thin=1)
summary(mglmm.aC)

prior3 <- list(R=list(V=diag(2)*1e-15, fix=1), G=list(G1=list(V=diag(2), nu=0.002)))
mglmm.ac <- MCMCglmm(cbind(c, a)~trait-1, random=~corg(trait):animal, rcov=~corg(trait):units,
                     pedigree=tree, reduced=TRUE, data=traits, prior=prior3, 
                     pr=TRUE, pl=TRUE, family=c("threshold", "threshold"), thin=1)
summary(mglmm.ac)


prior3 <- list(R=list(V=diag(3)*1e-15, fix=1), G=list(G1=list(V=diag(3), nu=0.002)))
mglmm.abc <- MCMCglmm(cbind(a,b,c)~trait-1, random=~corg(trait):animal, rcov=~corg(trait):units,
                     pedigree=tree, reduced=TRUE, data=traits, prior=prior3, 
                     pr=TRUE, pl=TRUE, family=c("threshold", "threshold", "threshold"), thin=1)
summary(mglmm.abc)
