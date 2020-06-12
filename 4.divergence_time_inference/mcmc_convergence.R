setwd("~/Documents/OneDrive - Universidade do Porto/BIOLOGIA/Doutoramento/Chapter 1 Phylogenomics/2-Pipeline/MCMCtree")

# code from https://github.com/mariodosreis/divtime/blob/master/R/analysis.R#L37

install.packages("coda")
library(coda)


mcmc1<-read.table("mcmc_9k_123_CDS_tree5_run1.txt",header=T)
mcmc2<-read.table("mcmc_9k_123_CDS_tree5_run2.txt",header=T)

# to check for convergence of the MCMC runs, we calculate the posterior
# means of times for each run, and plot them against each other
t.mean1 <- apply(mcmc1[,2:10], 2, mean) * 100
t.mean2 <- apply(mcmc2[,2:10], 2, mean) * 100
# good convergence is indicated when the points fall on the y = x line.

par(mfrow=c(2,2))
plot(t.mean1, t.mean2, main="a) Posterior times, r 1 vs. r 2"); abline(0, 1)
# notice that ancient times (t_n11 and t_n12) have small ESS
# trace plots are useful to visualise the MCMC and split problems

# trace plots are useful to visualise the MCMC and split problems
plot(mcmc1$t_n17, ty='l', main="b) trace of t_n17")
plot(mcmc1$t_n18, ty='l', main="b) trace of t_n18")
plot(mcmc1$t_n19, ty='l', main="b) trace of t_n19")
plot(mcmc1$t_n20, ty='l', main="b) trace of t_n20")
plot(mcmc1$t_n21, ty='l', main="b) trace of t_n21")
plot(mcmc1$t_n22, ty='l', main="b) trace of t_n22")
plot(mcmc1$t_n23, ty='l', main="b) trace of t_n23")
plot(mcmc1$t_n24, ty='l', main="b) trace of t_n24")
plot(mcmc1$t_n25, ty='l', main="b) trace of t_n25")
plot(mcmc1$t_n26, ty='l', main="b) trace of t_n26")
plot(mcmc1$t_n27, ty='l', main="b) trace of t_n27")

# we can calculate the effective sample sizes (ESS) of the parameters
# (you need to have the coda package installed for this to work)
mean.mcmc <- apply(mcmc1[,-1], 2, mean)
ess.mcmc <- apply(mcmc1[,-1], 2, coda::effectiveSize)
var.mcmc <- apply(mcmc1[,-1], 2, var)
se.mcmc <- sqrt(var.mcmc / ess.mcmc)
cbind(mean.mcmc, ess.mcmc, var.mcmc, se.mcmc)

mean.mcmc <- apply(mcmc2[,-1], 2, mean)
ess.mcmc <- apply(mcmc2[,-1], 2, coda::effectiveSize)
var.mcmc <- apply(mcmc2[,-1], 2, var)
se.mcmc <- sqrt(var.mcmc / ess.mcmc)
cbind(mean.mcmc, ess.mcmc, var.mcmc, se.mcmc)
