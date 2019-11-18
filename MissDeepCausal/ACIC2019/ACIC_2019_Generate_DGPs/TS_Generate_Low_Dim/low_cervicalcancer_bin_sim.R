low_cervicalcancer_bin_sim<-function(seed17,seed18,seed19,seed20){

#
# read in csv data file from https://archive.ics.uci.edu/ml/datasets/Cervical+cancer+%28Risk+Factors%29
dat<-na.omit(read.csv("risk_factors_cervical_cancer.csv",header=T))
# remove HPV as it is a very strong predictor for cervical cancer
dat<-dat[,-13]
#dim(dat)
#names(dat)
#rename outcome column 
outcome<-"Cancer"
names(dat)[which(names(dat)==outcome)]<-"Y"
#rename exposure column 
exposure<-"STDs"
names(dat)[which(names(dat)==exposure)]<-"A"
dat<-dat[,c(which(names(dat)=="Y"),which(names(dat)=="A"),which(names(dat)!="A"&names(dat)!="Y"))]
m.or <- glm(Y ~ ., data = dat, family = "binomial")
#fitting outcome model and pedicting counterfactual outcomes  
Q <- cbind(QAW = predict(m.or, type = "response"),
           Q0W = predict(m.or, newdata = data.frame(A = 0, dat[,-2]), type = "response"),
           Q1W = predict(m.or, newdata = data.frame(A = 1, dat[,-2]), type = "response"))
# Treatment effect in the original data
#diff(colMeans(Q[,-1]))  # 0.0287616
# Model the Pscore
g <- glm(A ~ ., data = dat[,-1], family = "binomial")
g1W <- predict(g, type = "response")
#summary(g1W)  
# Simulate the outcome and treatment assignment for each experimental dataset
# Modification 1 : very close to the logistic regression models fitted to the 
# actual data, however, higher prevalence of the exposure
n <- nrow(dat)
d.mod1 <- dat
# Generate treatment so that it is a little less problematic
beta <- coef(g)
g1W.mod1 <- as.numeric(plogis(cbind(1, as.matrix(d.mod1[,-(1:2)])) %*% (beta)) * 2.5)+0.05
#summary(g1W.mod1)
#0.09304 0.22010 0.25210 0.29330 0.30350 0.92100
set.seed(seed17)
d.mod1$A <- rbinom(n, 1, g1W.mod1)
# Use glm model as truth for this simulation
beta.Q.mod1 <- coef(m.or)
beta.Q.mod1[1] <- -3.2
beta.Q.mod1[2] <- 0
logit.drs.mod1 <-  cbind(1, as.matrix(dat[,-1])) %*% beta.Q.mod1
d.mod1$Y <- rbinom(n, 1, plogis(d.mod1$A*0.5 + logit.drs.mod1 ))
psi0.mod1 <- mean(plogis(0.5 + logit.drs.mod1)  - plogis(logit.drs.mod1))
#psi0.mod1
# 0.1066965
#table(d.mod1$A,d.mod1$Y)

#rename remaining variables  (covariates / confounders)
names(d.mod1)[which(names(d.mod1)!="A"&names(d.mod1)!="Y")]<-paste("V",1:(length(names(d.mod1))-2),sep="")

#simcc.mod1 <- calcMetrics(d.mod1, Y1 = plogis(0.5 + logit.drs.mod1), Y0 = plogis(logit.drs.mod1),	
#                          ps.full = g1W.mod1, simName = "cervical_cancer_mod1")
rownames(d.mod1)<-1:nrow(d.mod1)

DGP<-"lowdim_cervicalcancer_bin_mod1"
ATE<-psi0.mod1
filename<-paste("low",seed17,".csv",sep="")
gendat<-data.frame(DGP,filename,ATE)
write.csv(d.mod1,paste("C:/Users/Tibs/Dropbox/sim/simdat/low",seed17,".csv",sep=""),row.names = FALSE)
write.table(gendat,paste("C:/Users/Tibs/Dropbox/sim/gendatsummary",".csv",sep=""),append=T,row.names = FALSE,col.names = FALSE,sep=",")

#for(i in 1:100)
#{
#set.seed(i)
#ix<-sample(1:nrow(d.mod1),500,F) 
#write.csv(d.mod1[ix,],paste("C:/Users/Tibs/Dropbox/sim/simdata/cervicalcancer/lowdim_cervicalcancer_bin_mod1_",i,".csv",sep=""),row.names = FALSE)
#}

###################################################################################

# new variables
age.cycle <- sin(0.15*(dat$Age - mean(dat$Age))) # age cycle
#plot(smooth.spline(dat$Age,age.cycle),type="l")
risk.smoke <- sqrt(log(dat$Smokes*10+dat$Packyears+1))  # risk based on education and limit
#plot(dat$Packyears,risk.smoke)
risk.HC<-sqrt(log(dat$HC*2+dat$HCyears+1)+1)-1
#plot(dat$HCyears,risk.HC)


#PS model based on new transformed variables and interaction terms
g1W.mod2 <- plogis(.8 +log(1.2)*dat$Sexualpartners +.5*age.cycle  + 
                    risk.smoke + risk.HC-0.5*sqrt(dat$Pregnancies)-log(dat$First.sexual.intercourse))*0.9+0.05
#summary(g1W.mod2)
#hist(g1W.mod2)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.08177 0.13200 0.17010 0.21400 0.24460 0.93760 

#outcome model based on new transformed variables and existing covariates
logit.drs.mod2 <-  1 +risk.smoke + 0.5* +risk.HC+2*age.cycle-log(dat$First.sexual.intercourse)+
  sqrt(dat$Sexualpartners)-sqrt(dat$Pregnancies) 
d.mod2 <-dat 
set.seed(seed18)
d.mod2$A <- rbinom(n, 1, g1W.mod2)
d.mod2$Y <- rbinom(n, 1,plogis(d.mod2$A  + logit.drs.mod2))
# and the truth is:
psi0.mod2 <- mean(plogis(logit.drs.mod2 + 1) - plogis(logit.drs.mod2))
#psi0.mod2
# 0.1470471

#rename remaining variables  (covariates / confounders)
names(d.mod2)[which(names(d.mod2)!="A"&names(d.mod2)!="Y")]<-paste("V",1:(length(names(d.mod2))-2),sep="")

#simcc.mod2 <- calcMetrics(d.mod2, Y1 = plogis(1 + logit.drs.mod2), Y0 = plogis(logit.drs.mod2),	
#                          ps.full = g1W.mod2, simName = "cervical_cancer_mod2")

rownames(d.mod2)<-1:nrow(d.mod2)

DGP<-"lowdim_cervicalcancer_bin_mod2"
ATE<-psi0.mod2
filename<-paste("low",seed18,".csv",sep="")
gendat<-data.frame(DGP,filename,ATE)
write.csv(d.mod2,paste("C:/Users/Tibs/Dropbox/sim/simdat/low",seed18,".csv",sep=""),row.names = FALSE)
write.table(gendat,paste("C:/Users/Tibs/Dropbox/sim/gendatsummary",".csv",sep=""),append=T,row.names = FALSE,col.names = FALSE,sep=",")

#for(i in 101:200)
#{
#  set.seed(i)
#  ix<-sample(1:nrow(d.mod2),500,F) 
#  write.csv(d.mod2[ix,],paste("C:/Users/Tibs/Dropbox/sim/simdata/cervicalcancer/lowdim_cervicalcancer_bin_mod2_",i,".csv",sep=""),row.names = FALSE)
#}

###################################################################################

# Mod3. Poor overlap.  Go back to mod1. Keep same outcome, but change pscore model so that 
# it is more extreme because of IVs
d.mod3 <- dat
beta.mod3 <- coef(g)
beta.mod3[3]<--0.1
beta.mod3[4]<-0.25
beta.mod3[5]<-0.05
beta.mod3[10]<-2
g1W.mod3 <- as.vector(plogis(cbind(1, as.matrix(d.mod3[,-(1:2)])) %*% (beta.mod3 )))
#summary(g1W.mod3)  
# 0.08717 0.62000 0.71300 0.71380 0.81670 0.99080

set.seed(seed19)
d.mod3$A <- rbinom(n, 1, g1W.mod3)
#library(ROCR)
#res.rocr <- prediction(g1W.mod3, labels = d.mod3$A)
#performance(res.rocr, measure = "auc")@y.values[[1]]  #0.71

#regnoIV<-glm(g1W.mod3~.,data=d.mod3[,-c(1,2,3,4,5,10)], family ="binomial")
#g1W.mod3woIV<-predict(regnoIV, type ="response")
#summary(g1W.mod3woIV)

# setting some covariate effects (on Y) to 0 to introcuce instruments
beta.Q.mod3 <- beta.Q.mod1
beta.Q.mod3[1]<--3 # increase prevalence of Y
beta.Q.mod3[c(4,5,6,11)] <- 0
logit.drs.mod3 <-  cbind(1, as.matrix(dat[,-1])) %*% beta.Q.mod3
d.mod3$Y <- rbinom(n, 1, plogis(d.mod3$A*0.8 + logit.drs.mod3 ))
psi0.mod3 <- mean( plogis(0.8 + logit.drs.mod3 ) - plogis(logit.drs.mod3))
#psi0.mod3
# 0.1207257

#rename remaining variables  (covariates / confounders)
names(d.mod3)[which(names(d.mod3)!="A"&names(d.mod3)!="Y")]<-paste("V",1:(length(names(d.mod3))-2),sep="")

#simcc.mod3 <- calcMetrics(d.mod3, Y1 = plogis(0.8 + logit.drs.mod3), Y0 = plogis(logit.drs.mod3),	
#                          ps.full = g1W.mod3,ps.noIV=g1W.mod3woIV, simName = "cervical_cancer_mod3")

rownames(d.mod3)<-1:nrow(d.mod3)
DGP<-"lowdim_cervicalcancer_bin_mod3"
ATE<-psi0.mod3
filename<-paste("low",seed19,".csv",sep="")
gendat<-data.frame(DGP,filename,ATE)
write.csv(d.mod3,paste("C:/Users/Tibs/Dropbox/sim/simdat/low",seed19,".csv",sep=""),row.names = FALSE)
write.table(gendat,paste("C:/Users/Tibs/Dropbox/sim/gendatsummary",".csv",sep=""),append=T,row.names = FALSE,col.names = FALSE,sep=",")


#for(i in 201:300)
#{
#  set.seed(i)
#  ix<-sample(1:nrow(d.mod3),500,F) 
#  write.csv(d.mod3[ix,],paste("C:/Users/Tibs/Dropbox/sim/simdata/cervicalcancer/lowdim_cervicalcancer_bin_mod3_",i,".csv",sep=""),row.names = FALSE)
#}



#######################################################################
# Modification 4: Treatment effect heterogeneity + potential model misspecification
d.mod4 <- dat

# propensity score with manually specified coefficients
logit.g1W.mod4 <- 2.5 + abs(age.cycle) - .25*dat$Smokes - log(dat$Pregnancies+1)-dat$First.sexual.intercourse/(dat$Sexualpartners*8)

g1W.mod4 <- plogis(logit.g1W.mod4) # .05 to 0.95
#summary(g1W.mod4)
#0.1248  0.5954  0.7160  0.6829  0.8077  0.9386


# Q model with manually specified coefficients
logit.drs.mod4 <- -4 + age.cycle + dat$Smokes - log(dat$Pregnancies+1)+log(dat$First.sexual.intercourse+dat$Sexualpartners)
#summary(logit.drs.mod4)

set.seed(seed20)
d.mod4$A <- rbinom(n, 1, g1W.mod4)
# transformation of age to [0,1] using ecdf
agetrans <- log(rank(d.mod4$Age))/6
# linear predictor for Y with age - treatment interaction 
d.mod4$Y <- rbinom(n, 1, plogis(d.mod4$A  + .1*d.mod4$A * agetrans + .1 * d.mod4$A * age.cycle +   logit.drs.mod4 ))
psi0.mod4 <- mean(plogis(1 + .1 * agetrans  + .1*age.cycle + logit.drs.mod4 ) - plogis(logit.drs.mod4))
#psi0.mod4
# [1] 0.1567597

#c(psi0.mod4, mean(d.mod4$Y[d.mod4$A == 1]) - mean(d.mod4$Y[d.mod4$A == 0]))

# C statistic
#res.rocr.mod4 <- prediction(g1W.mod4, labels = d.mod4$A)
#performance(res.rocr.mod4, measure = "auc")@y.values[[1]]  # 0.69

names(d.mod4)[which(names(d.mod4)!="A"&names(d.mod4)!="Y")]<-paste("V",1:(length(names(d.mod4))-2),sep="")

#simcc.mod4 <- calcMetrics(d.mod4, Y1 = plogis(1 +.1*agetrans+.1*age.cycle +logit.drs.mod4), Y0 = plogis(logit.drs.mod4),	
#                          ps.full = g1W.mod4, simName = "cervical_cancer_mod4")

rownames(d.mod4)<-1:nrow(d.mod4)
DGP<-"lowdim_cervicalcancer_bin_mod4"
ATE<-psi0.mod4
filename<-paste("low",seed20,".csv",sep="")
gendat<-data.frame(DGP,filename,ATE)
write.csv(d.mod4,paste("C:/Users/Tibs/Dropbox/sim/simdat/low",seed20,".csv",sep=""),row.names = FALSE)
write.table(gendat,paste("C:/Users/Tibs/Dropbox/sim/gendatsummary",".csv",sep=""),append=T,row.names = FALSE,col.names = FALSE,sep=",")

#for(i in 301:400)
#{
#  set.seed(i)
#  ix<-sample(1:nrow(d.mod4),500,F) 
#  write.csv(d.mod4[ix,],paste("C:/Users/Tibs/Dropbox/sim/simdata/cervicalcancer/lowdim_cervicalcancer_bin_mod4_",i,".csv",sep=""),row.names = FALSE)
#}
}