## # Threshold models

## # PHYLIP & RPhylip
## Follow the directions to install phylip and Rphylip. 
## You will need the following packages:
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

## Set up the covariance matrix for changes in the underlying liabilities. Choose the values you want to 
## use below. Remember the matrix must be positive definite (how do you check?) 
covAB <- ??
covAC <- ??
covBC <- ??
R <- matrix(c(  1,   covAB, covAC,
              covAB,     1, covBC,
              covAC,  covBC,  1), ncol=3, nrow=3, byrow=TRUE)
colnames(R) <- rownames(R) <- c("A", "B", "C")
R 
eigen(R)

## Simulate the data using 'sim.char()'. 
dat <- sim.char(tree, par=R, nsim=1, model="BM")[,,1]
colnames(dat) <- c("A", "B", "C")
pairs(dat)

## Challenge: Just for kicks, realize that we don't really need a fancy function to 
## simulate multivariate Brownian motion on a tree, we can just use matrix algebra. 
## Here's how to start:
library(MASS)
## Here are our parameters:
root <- c(rep(0, ntaxa), rep(0, ntaxa), rep(0, ntaxa))
bigR <- kronecker(R, vcv.phylo(tree)) 
## The Kronecker product expands our 3 x 3 rate matrix and our 100 x 100 phylogenetic matrix
## into a 300 x 300 block matrix. 

## We are now ready to draw our random tip data. What function will we use? 
## (Hint: how are the data distributed on the phylogeny?)

dat2 <- ??
dat2 <- matrix(dat2, ncol=3)
pairs(dat2) ## This should look like 'pairs(dat)'

## Match the tree and to the liabilities and create discretized variables. 
td <- make.treedata(tree, dat)

## A pretty way to plot the liabilities:
X <- as.matrix(td[c("A", "B", "C")])
rownames(X) <- td$phy$tip.label
fancyTree(td$phy, type="scattergram", X=X)

## We now convert the liabilities into binary characters:
td <- mutate(td, a = as.numeric(A > median(A)), b = as.numeric(B > median(B)), c = as.numeric(C > median(C)))
tree <- td$phy
traits <- as.data.frame(td$dat)
rownames(traits) <- td$phy$tip.label

## Plot out the discrete traits
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
library(coda) # This package lets us plot MCMC results nicely, as well as lots of diagnostic features
ngen <- 20000
## First, let's run two continuous traits:
X <- as.matrix(traits[,c('A', 'C')])
mod.AC <- threshBayes(tree, X, types=c("cont", "cont"), ngen=ngen)
mcmc.AC <- mcmc(mod.AC$par)
plot(mcmc.AC)

## Now, change run again with one discrete and one continuous:
X <- as.matrix(traits[,c('a', 'C')])
mod.aC <- threshBayes(tree, X, types=c("disc", "cont"), ngen=ngen)
mcmc.aC <- mcmc(mod.aC$par)
plot(mcmc.aC)

## Finally, run with two discrete characters:
X <- as.matrix(traits[,c('a', 'c')])
mod.ac <- threshBayes(tree, X, types=c("disc", "disc"), ngen=ngen)
mcmc.ac <- mcmc(mod.ac$par)
plot(mcmc.ac)

## In all cases, the correlations should be the same. What happens? What is the result of 
## discretization on our estimate of the correlation?
par(mfrow=c(1,2))
postburnin <- floor(0.5*nrow(mcmc.AC)):nrow(mcmc.AC)
plot(density(mcmc.AC[postburnin,'r']), main="Posterior for correlation", xlab= "Correlation", xlim=c(-1,1))
lines(density(mcmc.aC[postburnin,'r']), col="blue")
lines(density(mcmc.ac[postburnin,'r']), col="green")

plot(0,0, type="n", main="Trace", xlab= "Generation", ylab="Correlation", ylim=c(-1,1), xlim=c(0, length(postburnin)))
lines(mcmc.AC[postburnin,'r'])
lines(mcmc.aC[postburnin,'r'], col="blue")
lines(mcmc.ac[postburnin,'r'], col="green")

## Now, use the function "effectiveSize" from the coda package on the postburnin sample of r. Has your mcmc chain converged?

## # MCMCglmmRAM
## Jarrod Hadfield has developed a more efficient algorithm for estimating the threshold model using 
## GLMM's. This is essentially the same as our animal model from quantitative genetics, but with the
## heritability set to 1, but instead of pedigree, we use our tree. 
install.packages("MCMCglmmRAM_2.22.tar.gz", repos=NULL, type="source")
library(MCMCglmmRAM)
## The data frame must have a column named "animal" and it MUST BE A FACTOR (not character).
traits$animal <- factor(tree$tip.label)
prior1 <- list(R=list(V=diag(2)*1e-15, fix=1), G=list(G1=list(V=diag(2), nu=0.002)))
mglmm.AC <- MCMCglmm(cbind(A, C)~trait-1, random=~us(trait):animal, rcov=~us(trait):units,
                     pedigree=tree, reduced=TRUE, data=traits, prior=prior1, 
                     pr=TRUE, pl=TRUE, family=c("gaussian", "gaussian"), thin=1)
summary(mglmm.AC)

prior2 <- list(R=list(V=diag(2)*1e-15, fix=1), G=list(G1=list(V=diag(2), nu=0.002, fix=2)))
mglmm.aC <- MCMCglmm(cbind(a, C)~trait-1, random=~us(trait):animal, rcov=~us(trait):units,
                     pedigree=tree, reduced=TRUE, data=traits, prior=prior2, 
                     pr=TRUE, pl=TRUE, family=c("threshold", "gaussian"), thin=1)
summary(mglmm.aC)

prior3 <- list(R=list(V=diag(2)*1e-15, fix=1), G=list(G1=list(V=diag(2), nu=0.002)))
mglmm.ac <- MCMCglmm(cbind(a, c)~trait-1, random=~corg(trait):animal, rcov=~corg(trait):units,
                     pedigree=tree, reduced=TRUE, data=traits, prior=prior3, 
                     pr=TRUE, pl=TRUE, family=c("threshold", "threshold"), thin=1)
summary(mglmm.ac)

## Now compare the effective sizes for r from the mglmm objects. Hint: look in the '$VCV' 
## object. 

## Now, run all 3 discrete traits at once. 
prior4 <- ??
mglmm.abc <- ??
summary(mglmm.abc)

## Notice that we said pr=TRUE and pl=TRUE, this saves the posterior distribution of random effects
## and latent variables, respectively, so that we can recover the estimated liabilities at the nodes.
## First, let's get the liabilities at the tips:
postburnin <- floor(0.3*nrow(mglmm.ac$Liab)):nrow(mglmm.ac$Liab)
a.liab <- apply(mglmm.ac$Liab[postburnin,1:100], 2, mean)
c.liab <- apply(mglmm.ac$Liab[postburnin,101:200], 2, mean)

plot(a.liab, traits$A, pch=21,  bg=traits$a)
plot(c.liab, traits$C, pch=21,  bg=traits$c)

## Now get the liabilities at the nodes.
a.nodeLiab <- apply(mglmm.ac$Sol[postburnin, 1:tree$Nnode], 2, mean)

## And we can plot them on the tree:
plot(tree)
a.pie <- (a.nodeLiab-min(a.nodeLiab))/diff(range(a.nodeLiab))
a.pie <- cbind(a.pie, 1-a.pie)
nodelabels(pie=a.pie, piecol = c("white", "black"), cex=0.5)
tiplabels(pch=21, bg = c("white", "black")[traits$a+1])

