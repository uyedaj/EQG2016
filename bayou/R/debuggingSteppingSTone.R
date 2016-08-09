devtools::load_all("~/repos/bayou/bayou_1.0/")
setwd("~/repos/EQG2016/bayou/R/")
library(treeplyr)

## # Preparing data for analysis with the R package treeplyr.
## We will use data from the Species 360 database of normal animal reference ranges for a variety of 
## traits (formerly known as the International Species Information System = ISIS). The data file has
## been slightly modified here, but the original data file is available at: 
## http://www2.isis.org/support/MEDARKS/Pages/Reference%20Ranges.aspx
sp360 <- readRDS("../data/species360.rds")
tbl_df(sp360)

## The tree file we are using is a mashup of several trees spliced together that includes 29,301 species of
## vertebrate. 
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
# Run a bayou analysis by first setting up a prior distribution for parameters. First we're going to make sure
## we resolve polytomies and set zero length branches to some small value. I've also found it is easier to set 
## priors if we scale the data. Note that to understand the parameter values, you must back-translate them using
## the scaling coefficients.
tree <- tdW$phy
dat <- tdW[[1]]
ME <- tdW[['lnMass.SD']]/sqrt(tdW[['N']])


prior <- readRDS("../output/runs/modelBB/priorBB.rds")
mcmc <- readRDS("../output/runs/modelBB/mcmcBB.rds")
chain <- readRDS("../output/runs/modelBB/chainBB.rds")

SE=0; ngen=10000; samp=100; chunk=100; control=NULL; tuning=NULL; new.dir=TRUE; plot.freq=NULL; outname="bayou"; ticker.freq=1000; tuning.int=NULL; moves=NULL; control.weights=NULL; lik.fn=NULL; plot.fn=NULL
model <- "OUrepar"; plot.fn <- NULL; pred=NULL; startpar <- NULL
Bk <- qbeta(seq(0,1, length.out=5), 0.3,1)

C_weightmatrix(cache, )
