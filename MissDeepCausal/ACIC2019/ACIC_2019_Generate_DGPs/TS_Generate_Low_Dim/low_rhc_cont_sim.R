low_rhc_cont_sim<-function(seed1,seed2,seed3,seed4)
 {

# read in csv data file from http://biostat.mc.vanderbilt.edu/wiki/pub/Main/DataSets/rhc.html
dat<-na.omit(read.csv("rhc.csv",header=T))
#dim(dat)
#names(dat)

#rename outcome column 
outcome<-"aps1"
names(dat)[which(names(dat)==outcome)]<-"Y"

#rename exposure column 
exposure<-"ca"
names(dat)[which(names(dat)==exposure)]<-"A"
dat$A<-(dat$A!="No")*1
dat$income<-as.numeric(dat$income)-1
dat$sex<-as.numeric(dat$sex)-1
dat$swang1<-as.numeric(dat$swang1)-1
dat$dnr1<-as.numeric(dat$dnr1)-1
dat$race<-as.numeric(dat$race=="white")*1
dat$death<-as.numeric(dat$death=="Yes")*1
dat$ninsclas<-(dat$ninsclas=="Private")*1
dat<-dat[,c(which(names(dat)=="Y"),which(names(dat)=="A"),which(names(dat)!="A"&names(dat)!="Y"))]

m.or <-lm(Y ~ ., data =dat)
   
#fitting outcome model and pedicting counterfactual outcomes  

Q <- cbind(QAW = predict(m.or),
           Q0W = predict(m.or, newdata = data.frame(A = 0, dat[,-2])),
           Q1W = predict(m.or, newdata = data.frame(A = 1, dat[,-2])))


# Treatment effect in the original data
#diff(colMeans(Q[,-1]))  # 2.560737 


# Model the Pscore
g <- glm(A ~ ., data = dat[,-1], family = "binomial")
g1W <- predict(g, type = "response")
#summary(g1W)  
#hist(g1W)

# Simulate the outcome and treatment assignment for each experimental dataset
# Modification 1 : very close to the logistic regression models fitted to the 
# actual data, however, higher prevalence of the exposure
n <- nrow(dat)
d.mod1 <- dat
# Generate treatment so that it is a little less problematic
beta <- coef(g)

g1W.mod1 <- as.numeric(plogis(cbind(1, as.matrix(d.mod1[,-(1:2)])) %*% beta) * 1.1)+0.05
#summary(g1W.mod1)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.05476 0.18590 0.27790 0.31010 0.40640 0.92390 

set.seed(seed1)
d.mod1$A <- rbinom(n, 1, g1W.mod1)
# Use glm model as truth for this simulation
beta.Q.mod1 <- coef(m.or)
beta.Q.mod1[1] <- -3.2
beta.Q.mod1[2] <- 0
lp.drs.mod1 <-  cbind(1, as.matrix(dat[,-1])) %*% beta.Q.mod1
psi0.mod1 <- 4
d.mod1$Y <- rnorm(n, d.mod1$A*psi0.mod1 + lp.drs.mod1,sd(m.or$residuals))

#rename remaining variables  (covariates / confounders)
names(d.mod1)[which(names(d.mod1)!="A"&names(d.mod1)!="Y")]<-paste("V",1:(length(names(d.mod1))-2),sep="")

#simcc.mod5 <- calcMetrics(d.mod1, Y1 = (psi0.mod1 + lp.drs.mod1), Y0 = lp.drs.mod1,	
#                          ps.full = g1W.mod1, simName = "rhc2_mod1")

rownames(d.mod1)<-1:nrow(d.mod1)
#for(i in 1:100)
#{
#  set.seed(i)
#  ix<-sample(1:nrow(d.mod1),1200,F) 
#  write.csv(d.mod1[ix,],paste("C:/Users/Tibs/Dropbox/sim/simdata/rhc2/lowdim_rhc2_cont_mod1_",i,".csv",sep=""),row.names = FALSE)
#}
DGP<-"lowdim_rhc_cont_mod1"
ATE<-psi0.mod1
filename<-paste("low",seed1,".csv",sep="")
gendat<-data.frame(DGP,filename,ATE)
write.csv(d.mod1,paste("C:/Users/Tibs/Dropbox/sim/simdat/low",seed1,".csv",sep=""),row.names = FALSE)
write.table(gendat,paste("C:/Users/Tibs/Dropbox/sim/gendatsummary",".csv",sep=""),append=T,row.names = FALSE,col.names = FALSE,sep=",")

###################################################################################

# new variables
age.cycle <- sin(0.1*(dat$age - mean(dat$age))) # age cycle
#plot(smooth.spline(dat$age,age.cycle),type="l")

risk.1<- log((dat$swang1+dat$wblc1/(dat$hrt1+1))+1) # risk based on education and limit
#hist(risk.1)
risk.2<-(rowSums(dat[,24:28]))/((rowSums(dat[,24:28]))+100)
#hist(risk.2)
#plot(rowSums(dat[,24:28]),risk.2)

#PS model based on new transformed variables and interaction terms
g1W.mod2<- plogis(.8 -log(1.2)*sqrt(dat$pot1)+2.5*age.cycle  + 
                    risk.1 + 1.3*risk.2-0.5*(sqrt(dat$sod1)*log(dat$pafi1+1))/100)*0.9631
#summary(g1W.mod2)
#hist(g1W.mod2)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.08072 0.42250 0.83180 0.68460 0.91990 0.96220   


#outcome model based on new transformed variables and existing covariates
lp.drs.mod2 <- (log(sqrt(dat$pot1))+age.cycle+risk.1 + risk.2+0.25*(sqrt(dat$sod1))/log(dat$pafi1+1))*20 
set.seed(seed2)
d.mod2 <-dat 
d.mod2$A <- rbinom(n, 1, g1W.mod2)
psi0.mod2 <- 10
d.mod2$Y <- rnorm(n,d.mod2$A*psi0.mod2+lp.drs.mod2,sd(m.or$residuals))


#rename remaining variables  (covariates / confounders)
names(d.mod2)[which(names(d.mod2)!="A"&names(d.mod2)!="Y")]<-paste("V",1:(length(names(d.mod2))-2),sep="")

#simcc.mod6<- calcMetrics(d.mod2, Y1 = psi0.mod2 + lp.drs.mod2, Y0 = lp.drs.mod2,	
#                          ps.full = g1W.mod2, simName = "rhc2_mod2")

rownames(d.mod2)<-1:nrow(d.mod2)
DGP<-"lowdim_rhc_cont_mod2"
ATE<-psi0.mod2
filename<-paste("low",seed2,".csv",sep="")
gendat<-data.frame(DGP,filename,ATE)
write.csv(d.mod2,paste("C:/Users/Tibs/Dropbox/sim/simdat/low",seed2,".csv",sep=""),row.names = FALSE)
write.table(gendat,paste("C:/Users/Tibs/Dropbox/sim/gendatsummary",".csv",sep=""),append=T,row.names = FALSE,col.names = FALSE,sep=",")

#rownames(d.mod2)<-1:nrow(d.mod2)
#for(i in 101:200)
#{
#  set.seed(i)
#  ix<-sample(1:nrow(d.mod2),1200,F) 
#  write.csv(d.mod2[ix,],paste("C:/Users/Tibs/Dropbox/sim/simdata/rhc2/lowdim_rhc2_cont_mod2_",i,".csv",sep=""),row.names = FALSE)
#}
###################################################################################

# Mod3. Poor overlap.  Go back to mod1. Keep same outcome, but change pscore model so that 
# it is more extreme because of IVs
d.mod3 <- dat
beta.mod3 <- coef(g)
beta.mod3[2]<-0.05
beta.mod3[3]<-+0.1
beta.mod3[4]<--0.1
beta.mod3[5]<-+0.1
beta.mod3[6]<--0.1
beta.mod3[7]<--0.1
beta.mod3[8]<--0.01
g1W.mod3 <- as.vector(plogis(cbind(1, as.matrix(d.mod3[,-(1:2)])) %*% (beta.mod3 )))*0.80+0.15
#summary(g1W.mod3)  
#0.1500  0.2621  0.6905  0.5970  0.8922  0.9496
#hist(g1W.mod3)

#regnoIV<-glm(g1W.mod3~.,data=d.mod3[,-c(1:9)], family ="binomial")
#g1W.mod3woIV<-predict(regnoIV, type ="response")
#summary(g1W.mod3woIV)

set.seed(seed3)
d.mod3$A <- rbinom(n, 1, g1W.mod3)
#library(ROCR)
#res.rocr <- prediction(g1W.mod3, labels = d.mod3$A)
#performance(res.rocr, measure = "auc")@y.values[[1]]  #0.85

# setting some covariate effects (on Y) to 0 to introcuce instruments
beta.Q.mod3 <- beta.Q.mod1
beta.Q.mod3[1]<-450
beta.Q.mod3[c(2:9)] <- 0
lp.drs.mod3 <-cbind(1, as.matrix(dat[,-1])) %*% beta.Q.mod3
psi0.mod3 <- 2
d.mod3$Y <- rnorm(n,d.mod3$A*psi0.mod3 + lp.drs.mod3,sd(m.or$residuals))


#rename remaining variables  (covariates / confounders)
names(d.mod3)[which(names(d.mod3)!="A"&names(d.mod3)!="Y")]<-paste("V",1:(length(names(d.mod3))-2),sep="")

#simcc.mod7 <- calcMetrics(d.mod3, Y1 = psi0.mod3 + lp.drs.mod3, Y0 = lp.drs.mod3,	
#                          ps.full = g1W.mod3,ps.noIV=g1W.mod3woIV, simName = "rhc2_mod3")

rownames(d.mod3)<-1:nrow(d.mod3)
DGP<-"lowdim_rhc_cont_mod3"
ATE<-psi0.mod3
filename<-paste("low",seed3,".csv",sep="")
gendat<-data.frame(DGP,filename,ATE)
write.csv(d.mod3,paste("C:/Users/Tibs/Dropbox/sim/simdat/low",seed3,".csv",sep=""),row.names = FALSE)
write.table(gendat,paste("C:/Users/Tibs/Dropbox/sim/gendatsummary",".csv",sep=""),append=T,row.names = FALSE,col.names = FALSE,sep=",")

#rownames(d.mod3)<-1:nrow(d.mod3)
#for(i in 201:300)
#{
#  set.seed(i)
#  ix<-sample(1:nrow(d.mod3),1200,F) 
#  write.csv(d.mod3[ix,],paste("C:/Users/Tibs/Dropbox/sim/simdata/rhc2/lowdim_rhc2_cont_mod3_",i,".csv",sep=""),row.names = FALSE)
#}

#######################################################################
# Modification 4: Treatment effect heterogeneity + potential model misspecification
d.mod4 <- dat

# propensity score with manually specified coefficients
logit.g1W.mod4 <- 2.5 + 2*abs(age.cycle) - 1*log(dat$sex/(dat$edu+1)+1) - 0.8*log(dat$hrt1+1)-0.01*(dat$temp1/(dat$swang1+1)/10)

g1W.mod4 <- plogis(logit.g1W.mod4)-0.05 # .05 to 0.95
#summary(g1W.mod4)
#0.07303 0.31950 0.47430 0.45370 0.56900 0.93860
#hist(g1W.mod4)
# Q model with manually specified coefficients
lp.drs.mod4 <- ((age.cycle + dat$sex*dat$hrt1 - log(dat$temp1)+dat$temp1^2+dat$swang+dat$edu+dat$swang*dat$edu)-1500)/200
#summary(lp.drs.mod4)

set.seed(seed4)
d.mod4$A <- rbinom(n, 1, g1W.mod4)
# transformation of age to [0,1] using logistic function
agetrans <- plogis(dat$age/10)
# linear predictor for Y with age - treatment interaction 
d.mod4$Y <- rnorm(n, 2*d.mod4$A  + .5*d.mod4$A * agetrans + .2 * d.mod4$A * age.cycle +   lp.drs.mod4,sd(m.or$residuals))
psi0.mod4 <- mean(2 + .5 * agetrans  + .2*age.cycle + lp.drs.mod4 ) - mean(lp.drs.mod4)
#psi0.mod4
# [1] 2.52739

# C statistic
#res.rocr.mod4 <- prediction(g1W.mod4, labels = d.mod4$A)
#performance(res.rocr.mod4, measure = "auc")@y.values[[1]]  # 0.71

names(d.mod4)[which(names(d.mod4)!="A"&names(d.mod4)!="Y")]<-paste("V",1:(length(names(d.mod4))-2),sep="")

rownames(d.mod4)<-1:nrow(d.mod4)
DGP<-"lowdim_rhc_cont_mod4"
ATE<-psi0.mod4
filename<-paste("low",seed4,".csv",sep="")
gendat<-data.frame(DGP,filename,ATE)
write.csv(d.mod4,paste("C:/Users/Tibs/Dropbox/sim/simdat/low",seed4,".csv",sep=""),row.names = FALSE)
write.table(gendat,paste("C:/Users/Tibs/Dropbox/sim/gendatsummary",".csv",sep=""),append=T,row.names = FALSE,col.names = FALSE,sep=",")


#simcc.mod4 <- calcMetrics(d.mod4, Y1 = (2 +.5*agetrans +.2*age.cycle + lp.drs.mod4), Y0 =lp.drs.mod4,	
#                          ps.full = g1W.mod4, simName = "rhc2_mod4")


#rownames(d.mod4)<-1:nrow(d.mod4)
#for(i in 301:400)
#{
#  set.seed(i)
#  ix<-sample(1:nrow(d.mod4),1200,F) 
#  write.csv(d.mod4[ix,],paste("C:/Users/Tibs/Dropbox/sim/simdata/rhc2/lowdim_rhc2_cont_mod4_",i,".csv",sep=""),row.names = FALSE)
#}
}