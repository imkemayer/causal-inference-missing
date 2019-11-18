# ACIC 2019
# Simulating low dimensional data based on credit card fraud dataset
# Susan Gruber, sgruber@putnamds.com
# July 5, 2019
# Yeh, I. C., & Lien, C. H. (2009). The comparisons of data mining techniques for the predictive accuracy of probability of default of credit card clients. Expert Systems with Applications, 36(2), 2473-2480. 
# https://archive.ics.uci.edu/ml/datasets/default+of+credit+card+clients
# Step 1: modify the covariates to make them better suited for the Data challenge
# Step 2: model the outcome ange pscore in the original data
# Step 3:  Develop four DGPs based on the modified covariate set, W.
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
# Step 1. Create modified covariate set
# And estimate effect of marital status on defaulting on next payment  
# in the real dataset to use as a starting point.
#------------------------------------------------------------------
d <- read.csv("default of credit card clients.csv", skip = 1)

 # Drop ID column
 d <- d[,-1]
 # change sex and marital status to binary 0/1
 d$sex[d$sex == 2] <- 0  # 1 male, 0 female
 d$marriage[d$marriage != 1] <- 0
# column 24 is the outcome - default on next payment or not
# column 4 is treatment - marital status
d <- d[,c(24, 4, 1:3, 5:23)]  # 30000 x 24
 colnames(d)[1:2] <- c("Y","A")
 
 #------------------------------------------------------------------
# Step 2: model the outcome ange pscore in the original data
# Parametric outcome regression 
m.or <- glm(Y ~ ., data = d, family = "binomial")

# Model the Pscore
g <- glm(A ~ ., data = d[,-1], family = "binomial")

#------------------------------------------------------------------
# Step 3:  Develop four DGPs based on the modified covariate set, W.
# Simulate the outcome and treatment assignment 
#------------------------------------------------------------------
# Modification 1 : very close to the logistic regression models fitted to the actual data
set.seed(10)
n <- nrow(d)
d.mod1 <- d
colnames(d.mod1)[-(1:2)] <- paste0("V", 1:(ncol(d)-2))

# Generate treatment that is a little less extreme
beta <- coef(g)
g1W.mod1 <- plogis(cbind(1, as.matrix(d.mod1[,-(1:2)])) %*% (beta )) * .95 
d.mod1$A <- rbinom(n, 1, g1W.mod1)
# Generate the outcome
beta.Q.mod1 <- coef(m.or)
beta.Q.mod1[2] <- 0
logit.drs.mod1 <-  cbind(1, as.matrix(d[,-1])) %*% beta.Q.mod1
d.mod1$Y <- rbinom(n, 1, plogis(d.mod1$A + logit.drs.mod1 ))
# True ATE in the population
psi0.mod1 <- mean(plogis(1 + logit.drs.mod1)  - plogis(logit.drs.mod1))
# Generate the datasets
set.seed(1)
n.b <- 500
for (i in mod1_files){
	b <- sample(1:n, n.b, replace = TRUE)
	d.mod1$A[b] <-rbinom(n.b, 1, g1W.mod1[b])
	d.mod1$Y[b]<-  rbinom(n.b, 1, plogis(d.mod1$A + logit.drs.mod1)[b])
	write.csv(d.mod1[b,], file = paste0("creditCardMod1", i , ".csv"), row.names = FALSE)
}

####
# Modification 2. Model misspecification
# treatment and outcome  are a complex function of measured covariates
# (easiest way is to create unobserved covariates that are then used in simple models to simulate the data.)
age.cycle <- sin(d$age - mean(d$age))
 risk <- log((d$education + 1) / d$limit)
 paid <- rowMeans((d[,13:18] - d[,19:24]) > (.3 * d[,13:18]))
 
 ps.mod2 <- plogis(.5 + .2*risk + .8 * age.cycle  + .02 * d$sex * d$age + paid)
 
logit.drs.mod2 <-  1 + .2 *risk + 2 * age.cycle  + .5 * d$sex - .6 * paid
d.mod2 <-d 
colnames(d.mod2)[-(1:2)] <- paste0("V", 1:(ncol(d)-2))

d.mod2$A <- rbinom(n, 1, ps.mod2)
d.mod2$Y <- rbinom(n, 1,plogis(d.mod2$A  + logit.drs.mod2))
# True ATE in the population
psi0.mod2 <- mean(plogis(logit.drs.mod2 + 1) - plogis(logit.drs.mod2))
# Generate the datasets
set.seed(2)
niter <- 100
n.b <- 500
for (i in mod2_files){
	b <- sample(1:n, n.b, replace = TRUE)
	d.mod2$A[b] <-rbinom(n.b, 1, ps.mod2[b])
	d.mod2$Y[b]<-  rbinom(n.b, 1, plogis(d.mod2$A + logit.drs.mod2)[b])
	write.csv(d.mod2[b,], file = paste0("creditCardMod2", i , ".csv"), row.names = FALSE)
}

# Mod3. Poor overlap.  Go back to mod1. Keep same outcome, but change pscore model so that 
# it is more extreme because of IVs
set.seed(3)
d.mod3 <- d
colnames(d.mod3[-(1:2)]) <- paste0("V", 1:(ncol(d)-2))
beta.mod3 <- c(-1.8, runif(ncol(d)-2, -.1, .2) / apply(d[,-(1:2)], 2, sd))
beta.mod3[12:16] <- beta.mod3[12:16]* .7

g1W.mod3 <-  as.vector(plogis(cbind(1, as.matrix(d.mod3[,-(1:2)])) %*% (beta.mod3 ))) 
 set.seed(10) 
beta.Q.mod3 <- beta.Q.mod1
beta.Q.mod3[c(2:3, 18:24)] <- 0
logit.drs.mod3 <-  cbind(1, as.matrix(d[,-(1)])) %*% beta.Q.mod3
# True ATE in the population
psi0.mod3 <- mean( plogis(1 + logit.drs.mod3 ) - plogis(logit.drs.mod3))
# Generate the datasets
set.seed(3)
colnames(d.mod3)[-(1:2)] <- paste0("V", 1:(ncol(d)-2))
niter <- 100
n.b <- 500
for (i in mod3_files){
	b <- sample(1:n, n.b, replace = TRUE)
	d.mod3$A[b] <-rbinom(n.b, 1, g1W.mod3[b])
	d.mod3$Y[b]<-  rbinom(n.b, 1, plogis(d.mod3$A + logit.drs.mod3)[b])
	write.csv(d.mod3[b,], file = paste0("creditCardMod3", i , ".csv"), row.names = FALSE)
}

# Modification 4: Treatment effect heterogeneity + potential model misspecification
d.mod4 <- d
logit.g1W.mod4 <- 1 + .2 * d$sex + .6 * age.cycle - .25 * log(abs(d$bill_amt1) + 1) -  1 * (d$pay_2 < 0)
g1W.mod4 <- plogis(logit.g1W.mod4) # .05 to 0.95
logit.drs.mod4 <-  -1 - .2 *risk + .2 * age.cycle  + 1 * d.mod4$sex - .25 * log(abs(d.mod4$bill_amt1) + 1) - 1 * 10^-5 * d.mod4$bill_amt5
set.seed(10)
f.young <- ecdf(d.mod4$age)
young <- f.young(d.mod4$age)
# True ATE in the population
psi0.mod4 <- mean(plogis(1 + .2 * young  + .2*age.cycle + logit.drs.mod4 ) - plogis(logit.drs.mod4))
# Generate the datasets
set.seed(4)
colnames(d.mod4)[-(1:2)] <- paste0("V", 1:(ncol(d)-2))
niter <- 100
n.b <- 500
for (i in mod4_files){
	b <- sample(1:n, n.b, replace = TRUE)
	d.mod4$A[b] <-rbinom(n.b, 1, g1W.mod4[b])
	d.mod4$Y[b]<-  rbinom(n.b, 1, plogis(d.mod4$A * (1 + .2 * young  + .2*age.cycle) + logit.drs.mod4)[b])
	write.csv(d.mod4[b,], file = paste0("creditCardMod4", i , ".csv"), row.names = FALSE)
}


 