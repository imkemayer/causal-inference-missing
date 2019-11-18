# ACIC 2019
# Simulating high dimensional data based on epilepsy dataset
# Susan Gruber, sgruber@putnamds.com
# July 5, 2019
# Epilepsy dataset from UCI repository
# http://archive.ics.uci.edu/ml/datasets/Epileptic+Seizure+Recognition#
# Andrzejak RG, Lehnertz K, Rieke C, Mormann F, David P, Elger CE (2001)
#  Indications of nonlinear deterministic and finite dimensional structures in time series of brain 
#  electrical activity: Dependence on recording region and brain state, Phys. Rev. E, 64, 061907
# Step 1: modify the covariates to make them better suited for the Data challenge
# Step 2:  Develop four DGPs based on the modified covariate set, W.
#  		    and generate 100 datasets for each, with sample size = 500
# 				(Consider the 30,000 original observations to be the source population)
#     Mod 1. Simple, straightforward, parametric main terms models  
# 	   Mod 2. Treatment and outcome are complex functions of measured covariates
#     Mod 3.  Poor overlap
#     Mod 4: Treatment effect heterogeneity, complex models.
#------------------------------------------------------------------------------
#  R 3.6.0 changed the "sample" function
# This re-sets it to the old behavior, same as when the ACIC 2019 datasets
# were generated
if (getRversion() < '3.6.0') {
	RNGkind(sample.kind="Rounding")
}
##------------------------------------------------------------------
# Step 1. Read in the data and create modified covariate set
#------------------------------------------------------------------
 d <- read.csv("data.csv", header = TRUE)  
# drop ID in column 1, and outcome (seizure) in column 180
Y <- as.integer(d[,180] == 1)  # 20% are ones = 2300 of them.
d <- d[,-c(1,180)]
colnames(d) <- paste0("V", 1:ncol(d))
# sample from the time series to remove closely correlated columns
keepCoef<- seq(1, 178, 4)
n.coef <- length(keepCoef) # 45 of them
n <- nrow(d)
#------------------------------------------------------------------
# Step 2:  Develop four DGPs 
#------------------------------------------------------------------
#------------------------------------------------------------------
# Mod 1: treatment and outcome are main term functions of the 45 coefficients
#----------------------------------------------
set.seed(6)
beta.A.mod1 <- runif(n.coef, -.04, .094) / apply(d[,keepCoef], 2, sd)
logitA.mod1 <- -.4 + as.matrix(d[,keepCoef]) %*% beta.A.mod1
beta.Y.mod1 <- 3 * beta.A.mod1  # * sample(c(-1, 1), length(beta.A.mod1), replace = TRUE)
logit.drs.mod1 <- -.7 + as.matrix(d[,keepCoef]) %*% beta.Y.mod1
A.mod1 <- rbinom(n, 1, plogis(logitA.mod1))
Y.mod1 <- rbinom(n, 1, plogis(-3*A.mod1 + logit.drs.mod1))
# True ATE in the population
psi0.mod1 <- mean(  plogis(-3 + logit.drs.mod1) - plogis(logit.drs.mod1))
# Generate datasets
d.mod1 <- data.frame(Y= NA, A= NA, d)
set.seed(1)
niter <- 100
n.b <- 1500
for (i in mod1_files){
	b <- sample(1:n, n.b, replace = TRUE)
	d.mod1$A[b] <-rbinom(n.b, 1, plogis(logitA.mod1[b]))
	d.mod1$Y[b]<- rbinom(n.b, 1, plogis(-3*d.mod1$A[b]+ logit.drs.mod1[b]))
	write.csv(d.mod1[b,], file =  paste0("epilepsyMod1", i , ".csv"), row.names = FALSE)
}
#------------------------------------------------------------------
# Mod 2: model misspecification
#------------------------------------------------------------------
W.ratio <- log(abs(d[,keepCoef[-n.coef]+1]) + 1) / ( abs(d[,keepCoef[-n.coef] + 3]) + 1)  # between -6.5 and 6.5
# 11500 x 44
 random.cols <- c(1, 51, 81, 111, 151)
W.interact <- d[, random.cols] * d[,random.cols + 10]
temp.mod2 <- as.matrix(cbind(d[,c(random.cols, random.cols+10)],W.ratio, W.interact))
colnames(temp.mod2) <- paste0("V", 1:ncol(temp.mod2))
set.seed(21)
beta.A.mod2 <- runif(ncol(temp.mod2), -.02, .09) / apply(temp.mod2, 2, sd)
logitA.mod2 <- -.5 + temp.mod2 %*% beta.A.mod2
beta.Y.mod2 <-2.4* beta.A.mod2 
logit.drs.mod2 <- -1.8 + as.matrix(temp.mod2[,-1]) %*% beta.Y.mod2[-1] + as.matrix(d[,c(150, 160)]) %*% c(-.04, .06)
# true ATE in the population
psi0.mod2 <- mean(  plogis(3 + logit.drs.mod2) - plogis(logit.drs.mod2))  # .2165881
# generate datasets
d.mod2 <- data.frame(Y= NA, A= NA, d)
set.seed(2)
niter <- 100
n.b <- 1500
for (i in mod2_files){
	b <- sample(1:n, n.b, replace = TRUE)
	d.mod2$A[b] <-rbinom(n.b, 1, plogis(logitA.mod2)[b])
	d.mod2$Y[b]<- rbinom(n.b, 1, plogis(3*d.mod2$A[b] + logit.drs.mod2[b]))
	write.csv(d.mod2[b,], file = paste0("epilepsyMod2", i , ".csv"), row.names = FALSE)
}
#------------------------------------------------------------------
# Mod 3: treatment effect heterogeneity
#------------------------------------------------------------------
set.seed(3)
beta.A.mod3 <- runif(n.coef, -.03, .1) / apply(d[,keepCoef], 2, sd)
logitA.mod3 <- -.5 + as.matrix(d[,keepCoef]) %*% beta.A.mod3
beta.Y.mod3 <- 2 * beta.A.mod3 
logit.drs.mod3 <- -1 + as.matrix(d[,keepCoef]) %*% beta.Y.mod3
set.seed(23)
A.mod3 <- rbinom(n, 1, plogis(logitA.mod3))
Y.mod3 <- rbinom(n, 1, plogis( A.mod3 * (.5 + .6 *  (d[,keepCoef[30]] < mean(d[,keepCoef[30]])) + .8 * (d[,keepCoef[4]] < 0)) + logit.drs.mod3))
# true ATE in the population
psi0.mod3 <- mean(plogis( (.5 + .6 * (d[,keepCoef[30]] < mean(d[,keepCoef[30]])) + .8 * (d[,keepCoef[4]] < 0)) + logit.drs.mod3) - plogis(logit.drs.mod3))  
# Generate datasets
d.mod3 <- data.frame(Y= NA, A= NA, d)
set.seed(3)
niter <- 100
n.b <- 2000
for (i in mod3_files){
	b <- sample(1:n, n.b, replace = TRUE)
	d.mod3$A[b] <-rbinom(n.b, 1, plogis(logitA.mod3)[b])
	d.mod3$Y[b]<-  rbinom(n.b, 1, plogis( d.mod3$A * (.5 + .6 *  (d[,keepCoef[30]] < mean(d[,keepCoef[30]])) + .8 * (d[,keepCoef[4]] < 0)) + logit.drs.mod3)[b])
	write.csv(d.mod3[b,], file =  paste0("epilepsyMod3", i , ".csv"), row.names = FALSE)
}
#------------------------------------------------------------------
# mod 4: do treatment heterogeneity along with  IVs.
#------------------------------------------------------------------

set.seed(40)
beta.A.mod4 <- runif(ncol(temp.mod2), -.1, .12) / apply(temp.mod2, 2, sd)
logitA.mod4 <- -.1 + temp.mod2 %*% beta.A.mod4 
beta.Y.mod4 <- 2 * beta.A.mod4
beta.Y.mod4[1:5] <- 0  # create IVs
logit.drs.mod4 <- -1.8 + as.matrix(temp.mod2) %*% beta.Y.mod4+ as.matrix(d[,c(150, 160)]) %*% c(-.005, -0.02)
set.seed(21)
# true ATE in the population
psi0.mod4 <- mean(plogis( 2 + .01 * d[,160] + logit.drs.mod4) - plogis(logit.drs.mod4))    # 0.2916274
# Generate datasets
d.mod4 <- data.frame(Y= NA, A= NA, d)
set.seed(4)
niter <- 100
n.b <- 2000
for (i in mod4_files){
	b <- sample(1:n, n.b, replace = TRUE)
	d.mod4$A[b] <-rbinom(n.b, 1, plogis(logitA.mod4)[b])
	d.mod4$Y[b]<-  rbinom(n.b, 1, plogis(d.mod4$A * (2 + .01 * d[,160]) + logit.drs.mod4)[b])
	write.csv(d.mod4[b,], file =  paste0("epilepsyMod4", i , ".csv"), row.names = FALSE)
}


