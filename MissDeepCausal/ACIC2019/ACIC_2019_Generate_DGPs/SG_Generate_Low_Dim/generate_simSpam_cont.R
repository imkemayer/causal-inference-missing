# ACIC 2019
# Simulate low dimensional data based on spam email dataset
# (continuous outcome)
# Susan Gruber, sgruber@putnamds.com
# July 5, 2019
# Generate 100 datasets from each DGP
# 2. Sources:
   # (a) Creators: Mark Hopkins, Erik Reeber, George Forman, Jaap Suermondt
        # Hewlett-Packard Labs, 1501 Page Mill Rd., Palo Alto, CA 94304
   # (b) Donor: George Forman (gforman at nospam hpl.hp.com)  650-857-7835
   # (c) Generated: June-July 1999
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
# 
##------------------------------------------------------------------
# Step 1. Read in the data and create modified covariate set
#------------------------------------------------------------------
	d <- read.csv("spambase/spambase.csv", header =  FALSE)
	# columns 1- 54 are continuous word frequencies (percentages between 0 and 100) 
	# columns 55-57 count capital letters in various ways.
	# discard them, and  replace with  a binary indicator of more than the mean number of capital letters
	# column 58 is the indicator of whether message is spam or not
	# Create datasets of size 500 x 25 and let people analyze them.
	# Keep 20 from 1- 54, and also 55
	d$A  <- as.integer(d[,57] > mean(d[,57]))  # 24% are above the mean - think of this as treatment
	keep <- c(58,59,  5:25, 55)
	d <- d[,keep]
	d$x22 <- log(d$x22) 	
	colnames(d) <- c("Y", "A", paste0("V", 1:(ncol(d)-2)))	
#------------------------------------------------------------------
# Step 2:  Develop four DGPs 
#------------------------------------------------------------------
	n <- nrow(d)
# Generate treatment so it is less problematic
# ACIC 2019
# Simulate low dimensional data based on spam email dataset
# (binary outcome)
# Susan Gruber, sgruber@putnamds.com
# July 5, 2019
# Generate 100 datasets from each DGP
# 2. Sources:
   # (a) Creators: Mark Hopkins, Erik Reeber, George Forman, Jaap Suermondt
        # Hewlett-Packard Labs, 1501 Page Mill Rd., Palo Alto, CA 94304
   # (b) Donor: George Forman (gforman at nospam hpl.hp.com)  650-857-7835
   # (c) Generated: June-July 1999
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
# 
##------------------------------------------------------------------
# Step 1. Read in the data and create modified covariate set
#------------------------------------------------------------------
	d <- read.csv("spambase/spambase.csv", header =  FALSE)
	# columns 1- 54 are continuous word frequencies (percentages between 0 and 100) 
	# columns 55-57 count capital letters in various ways.
	# discard them, and  replace with  a binary indicator of more than the mean number of capital letters
	# column 58 is the indicator of whether message is spam or not
	# Create datasets of size 500 x 25 and let people analyze them.
	# Keep 20 from 1- 54, and also 55
	d$A  <- as.integer(d[,57] > mean(d[,57]))  # 24% are above the mean - think of this as treatment
	keep <- c(58,59,  5:25, 55)
	d <- d[,keep]
	d$x22 <- log(d$x22) 	
	colnames(d) <- c("Y", "A", paste0("V", 1:(ncol(d)-2)))	
#------------------------------------------------------------------
# Step 2:  Develop four DGPs 
#------------------------------------------------------------------
	n <- nrow(d)
#------------------------------------------------------------------
# Modification 1 : main terms parametric models work
#------------------------------------------------------------------
set.seed(100)
beta.mod1 <- c(-2.8, .33*colMeans(d[,-c(1:2, 21, 22, 24)])/ (apply(d[,-c(1:2, 21, 22, 24)], 2, sd) + rnorm(ncol(d)-5, .05, .03)))
g1W.mod1 <- plogis(cbind(1, as.matrix(d[,-c(1:2, 21, 22, 24)])) %*% beta.mod1)
# 
A.mod1 <- rbinom(n, 1, g1W.mod1)
beta.Q.mod1 <- runif(ncol(d), -3, 3)
beta.Q.mod1[2] <- 4  # coef on A
EY0.mod1 <- cbind(1, as.matrix(d[,-(1:2)])) %*% beta.Q.mod1[-2]
Y.mod1 <- EY0.mod1 + beta.Q.mod1[2] * A.mod1 + rnorm(n, 0, sd = 1.2)
# true population ATE
psi0.mod1 <- beta.Q.mod1[2]
d.mod1 <- d
# Generate datasets
set.seed(15)
niter <- 100
n.b <- 500
for (i in mod1_files){
	b <- sample(1:n, n.b, replace = TRUE)
	d.mod1$A[b] <- rbinom(n.b, 1, g1W.mod1[b])
	d.mod1$Y[b] <- EY0.mod1[b] + beta.Q.mod1[2] * d.mod1$A[b] + rnorm(n.b, 0, sd = 1.2)
	write.csv(d.mod1[b,], file = paste0("spam_contMod1", i , ".csv"), row.names = FALSE)
}

#------------------------------------------------------------------
# Modification 2.  Model misspecification
# Outcome  is  a complex function of measured covariates. 
# Pscore same as Mod1
#------------------------------------------------------------------
set.seed(5)
g1W.mod2 <- plogis(qlogis(g1W.mod1) * 1.15)
A.mod2 <-  rbinom(n, 1, g1W.mod2)
beta.Q.mod2 <- runif(ncol(d), -3, 12)
drs.mod2 <-  cbind(1, as.matrix(d[,-(1:4)])) %*% beta.Q.mod2[-(2:4)]  - 6* (d[,12] + 3) / (d[,22] + 1) + 3* d[,23] > 0
Y.mod2 <- 10 + 2.5 *A.mod2 + drs.mod2 + rnorm(n)
# true population ATE
psi0.mod2 <- 2.5
# Generate datasets
d.mod2 <- d
set.seed(16)
niter <- 100
n.b <- 500
for (i in mod2_files){
	b <- sample(1:n, n.b, replace = TRUE)
	d.mod2$A[b] <- rbinom(n.b, 1, g1W.mod2[b])
	d.mod2$Y[b] <- 10 + 2.5 *d.mod2$A[b] + drs.mod2[b] + rnorm(n.b)
	write.csv(d.mod2[b,], file = paste0("spam_contMod2", i , ".csv"), row.names = FALSE)
}

#------------------------------------------------------------------
# Modification 3. Treatment heterogeneity and
# non-linear pscore 
#------------------------------------------------------------------
set.seed(3)
W.mod3 <- d[,c(3,10, 14, 17, 22, 23)]
W.mod3.int <- cbind(int1 = (W.mod3[,1] + 2) * (W.mod3[,6] + 1), int2 = (W.mod3[,1] - 2) * (W.mod3[2] + 1) * (W.mod3[,5] + 3))
W.mod3.ratio <-  (W.mod3[,6] + 3) / (W.mod3[,1] + 1)
beta.mod3 <- runif(7, -.1, .32) / c(apply(W.mod3.int, 2, sd), sd(W.mod3.ratio), apply(W.mod3[,c(1, 3:5)], 2, sd))
logit.A.mod3 <- -.7 + as.matrix(cbind(W.mod3.int, W.mod3.ratio, W.mod3[,c(1, 3, 4,5)])) %*% beta.mod3

A.mod3 <- rbinom(n, 1, plogis(logit.A.mod3))
beta.Q.mod3 <- rnorm(ncol(d), 1, 2)
drs.mod3 <- cbind(1, as.matrix(d[,-(1:2)])) %*% beta.Q.mod3[-2]
Y.mod3 <- -3 +  3* A.mod3 * (1 - d[,21]  + 2.2* d[,10]  ) + drs.mod3 + rnorm(n)
# true population ATE
psi0.mod3 <- mean( 3* (1 - d[,21]  + 2.2* d[,10]  ) )
# Generate datasets
d.mod3 <- d
set.seed(17)
niter <- 100
n.b <- 500
for (i in mod3_files){
	b <- sample(1:n, n.b, replace = TRUE)
	d.mod3$A[b] <- rbinom(n.b, 1, plogis(logit.A.mod3[b]))
	d.mod3$Y[b] <- -3 +  3* d.mod3$A[b] * (1 - d[b,21]  + 2.2* d[b,10]  ) + drs.mod3[b] + rnorm(n.b)
	write.csv(d.mod3[b,], file = paste0("spam_contMod3", i , ".csv"), row.names = FALSE)
}
#------------------------------------------------------------------
# Modification 4. Model misspecification
# treatment and outcome  are a complex function of measured covariates
#------------------------------------------------------------------
set.seed(10)
W4.1 <- unlist((2*d[,13]+5 + (d[,22])) / (d[14]+2))
W4.2 <- cut(d[,17], 8, labels = FALSE)
W4.2 <- 1.1 * W4.2  - 3 * (W4.2 == 3) - .8 * (W4.2 == 2)
W4.3 <-  sin(d[,10] + d[,16]+1)
# Generate PS as function of d[,10:14] and W4.1 W4.2
logit.A.mod4 <- 1.65   - .3 * d[,12]  - .25 * d[,14] - .5 * W4.1 + .22*W4.2 + .58*W4.3
A.mod4  <- rbinom(n, 1, plogis(logit.A.mod4))
# outcome is function of d[,8:12] and W4.1, W4.2
drs.Y.mod4 <- 12 + 3*d[,8] + 1.5*d[,7] + 4.1*d[,10] - 1.5*d[,11] + 10 *d[,12]  + 3*W4.3 + 4*W4.2
Y.mod4 <-  A.mod4 * (5 + 2*log(W4.2 * W4.1 )) + drs.Y.mod4 + rnorm(n)
# true population ATE
psi0.mod4 <- mean(5 + 2*log(W4.2 * W4.1 ))
# Generate datasets
d.mod4 <- d
set.seed(18)
niter <- 100
n.b <- 500
for (i in mod4_files){
	b <- sample(1:n, n.b, replace = TRUE)
	d.mod4$A[b] <- rbinom(n.b, 1, plogis(logit.A.mod4[b]))
	d.mod4$Y[b] <- d.mod4$A[b] * (5 + 2*log(W4.2 * W4.1 )[b]) + drs.Y.mod4[b] + rnorm(n.b)
	write.csv(d.mod4[b,], file = paste0("spam_contMod4", i , ".csv"), row.names = FALSE)
}