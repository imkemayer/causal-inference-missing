low_studentperf2_cont_sim<-function(seed9,seed10,seed11,seed12){

# read in data from: https://archive.ics.uci.edu/ml/datasets/student+performance
dat<-read.table("student-por.csv",sep=";",header=TRUE)

#names(dat)
for(i in 1:ncol(dat))
{
if(i!=15&i!=30)
{
dat[,i]<-as.numeric(dat[,i])-1 
}
}

A<-dat[,4]
Y<-apply(dat[,31:33],1,mean)
dat<-cbind(Y,A,dat[,-c(4,31:33)])

m.or <- lm(Y ~ ., data =dat)
   
#fitting outcome model and pedicting counterfactual outcomes  

Q <- cbind(QAW = predict(m.or),
           Q0W = predict(m.or, newdata = data.frame(A = 0, dat[,-2])),
           Q1W = predict(m.or, newdata = data.frame(A = 1, dat[,-2])))

# Treatment effect in the original data
#diff(colMeans(Q[,-1]))  # 0.1734114


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

g1W.mod1 <- as.numeric(plogis(cbind(1, as.matrix(d.mod1[,-(1:2)])) %*% beta))*0.95
#summary(g1W.mod1)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.05716 0.50660 0.73960 0.66160 0.83650 0.91780

set.seed(seed9)
d.mod1$A <- rbinom(n, 1, g1W.mod1)
# Use glm model as truth for this simulation
beta.Q.mod1 <- coef(m.or)
beta.Q.mod1[2] <- 0
lp.drs.mod1 <-  cbind(1, as.matrix(dat[,-1])) %*% beta.Q.mod1
psi0.mod1<-2
d.mod1$Y <- rnorm(n,d.mod1$A*psi0.mod1 + lp.drs.mod1,sd(m.or$residuals))


#rename remaining variables  (covariates / confounders)
names(d.mod1)[which(names(d.mod1)!="A"&names(d.mod1)!="Y")]<-paste("V",1:(length(names(d.mod1))-2),sep="")

rownames(d.mod1)<-1:nrow(d.mod1)
DGP<-"lowdim_studentperf2_cont_mod1"
ATE<-psi0.mod1
filename<-paste("low",seed9,".csv",sep="")
gendat<-data.frame(DGP,filename,ATE)
write.csv(d.mod1,paste("C:/Users/Tibs/Dropbox/sim/simdat/low",seed9,".csv",sep=""),row.names = FALSE)
write.table(gendat,paste("C:/Users/Tibs/Dropbox/sim/gendatsummary",".csv",sep=""),append=T,row.names = FALSE,col.names = FALSE,sep=",")


#simcc.mod9 <- calcMetrics(d.mod1, Y1 = (psi0.mod1 + lp.drs.mod1), Y0 = lp.drs.mod1,
#                          ps.full = g1W.mod1, simName = "student_perf_mod1")


#for(i in 1:100)
#{
#  set.seed(i)
#  ix<-sample(1:nrow(d.mod1),500,F) 
#  write.csv(d.mod1[ix,],paste("C:/Users/Tibs/Dropbox/sim/simdata/studentperf/lowdim_studentperf_cont_mod1_",i,".csv",sep=""),row.names = FALSE)
#}


###################################################################################

# new variables
age.cycle <- (cos(dat$age+0.5*dat$age^2-2*dat$age^3 - mean(dat$age))) # age cycle
#plot(smooth.spline(dat$age,age.cycle),type="l")

risk.edu<-log((rowSums(dat[,8:11])-dat$famsize*3-dat$Pstatus*5)+10) # risk based on education 
#hist(risk.edu)

risk.ft<-(dat$freetime+1)/(dat$goout+1)
#hist(risk.ft)

#PS model based on new transformed variables and interaction terms
g1W.mod2 <- plogis(-2+risk.edu+0.25*age.cycle - risk.ft - 2*log(1/(dat$health+1)))*0.9+0.04
                    
#summary(g1W.mod2)
#hist(g1W.mod2)
# Min. 1st Qu.  Median    Mean 3rd Qu. Max. 
# 0.06962 0.63110 0.80110 0.71770 0.86610 0.91760 


#outcome model based on new transformed variables and existing covariates
#lp.drs.mod2 <- 3+sqrt(risk.edu)+age.cycle + risk.ft + log(dat$age/(dat$health+1))

lp.drs.mod2 <- (13+sqrt(risk.edu)+age.cycle + 10*sin(risk.ft) - log(dat$age/(dat$health+1)))/7+5
#summary(lp.drs.mod2)
d.mod2 <-dat 
set.seed(seed10)
d.mod2$A <- rbinom(n, 1, g1W.mod2)
psi0.mod2 <-2
d.mod2$Y <- rnorm(n, d.mod2$A*psi0.mod2+lp.drs.mod2,sd(m.or$residuals))
# and the truth is:

#rename remaining variables  (covariates / confounders)
names(d.mod2)[which(names(d.mod2)!="A"&names(d.mod2)!="Y")]<-paste("V",1:(length(names(d.mod2))-2),sep="")


#simcc.mod10 <- calcMetrics(d.mod2, Y1 = (psi0.mod2 + lp.drs.mod2), Y0 = lp.drs.mod2,
#                          ps.full = g1W.mod2, simName = "student_perf_mod2")


rownames(d.mod2)<-1:nrow(d.mod2)
DGP<-"lowdim_studentperf2_cont_mod2"
ATE<-psi0.mod2
filename<-paste("low",seed10,".csv",sep="")
gendat<-data.frame(DGP,filename,ATE)
write.csv(d.mod2,paste("C:/Users/Tibs/Dropbox/sim/simdat/low",seed10,".csv",sep=""),row.names = FALSE)
write.table(gendat,paste("C:/Users/Tibs/Dropbox/sim/gendatsummary",".csv",sep=""),append=T,row.names = FALSE,col.names = FALSE,sep=",")

#for(i in 101:200)
#{
#  set.seed(i)
#  ix<-sample(1:nrow(d.mod2),500,F) 
#  write.csv(d.mod2[ix,],paste("C:/Users/Tibs/Dropbox/sim/simdata/studentperf/lowdim_studentperf_cont_mod2_",i,".csv",sep=""),row.names = FALSE)
#}

###################################################################################

# Mod3. Poor overlap.  Go back to mod1. Keep same outcome, but change pscore model so that 
# it is more extreme because of IVs
d.mod3 <- dat
beta.mod3 <- coef(g)
beta.mod3[3]<-1.9
beta.mod3[4]<--0.3
beta.mod3[5]<--0.3
beta.mod3[6]<--0.3
beta.mod3[7]<-1.6
beta.mod3[8]<-0.6
g1W.mod3 <- as.vector(plogis(cbind(1, as.matrix(d.mod3[,-(1:2)])) %*% (beta.mod3 )))*0.9+0.05 
#summary(g1W.mod3)  # .05 to 0.95
# 0.0505  0.1174  0.4956  0.4996  0.8788  0.9482
#hist(g1W.mod3)

#regnoIV<-glm(g1W.mod3~.,data=d.mod3[,-c(1:9)], family ="binomial")
#g1W.mod3woIV<-predict(regnoIV, type ="response")
#summary(g1W.mod3woIV)
#0.03606 0.32690 0.52680 0.49960 0.68550 0.91130

set.seed(seed11)
d.mod3$A <- rbinom(n, 1, g1W.mod3)
#library(ROCR)
#res.rocr <- prediction(g1W.mod3, labels = d.mod3$A)
#performance(res.rocr, measure = "auc")@y.values[[1]]  #0.90


# setting some covariate effects (on Y) to 0 to introcuce instruments
beta.Q.mod3 <- beta.Q.mod1
beta.Q.mod3[1]<--0.5 # increase prevalence of Y
beta.Q.mod3[c(3:9)] <- 0
lp.drs.mod3 <-  cbind(1, as.matrix(dat[,-1])) %*% beta.Q.mod3
psi0.mod3<--2
d.mod3$Y <- rnorm(n, psi0.mod3*d.mod3$A + lp.drs.mod3,sd(m.or$residuals))


#rename remaining variables  (covariates / confounders)
names(d.mod3)[which(names(d.mod3)!="A"&names(d.mod3)!="Y")]<-paste("V",1:(length(names(d.mod3))-2),sep="")


#simcc.mod11 <- calcMetrics(d.mod3, Y1 = (psi0.mod3 + lp.drs.mod3), Y0 = lp.drs.mod3,
#                           ps.full = g1W.mod3,ps.noIV=g1W.mod3woIV,simName = "student_perf_mod3")

rownames(d.mod3)<-1:nrow(d.mod3)

DGP<-"lowdim_studentperf2_cont_mod3"
ATE<-psi0.mod3
filename<-paste("low",seed11,".csv",sep="")
gendat<-data.frame(DGP,filename,ATE)
write.csv(d.mod3,paste("C:/Users/Tibs/Dropbox/sim/simdat/low",seed11,".csv",sep=""),row.names = FALSE)
write.table(gendat,paste("C:/Users/Tibs/Dropbox/sim/gendatsummary",".csv",sep=""),append=T,row.names = FALSE,col.names = FALSE,sep=",")


#for(i in 201:300)
#{
#  set.seed(i)
#  ix<-sample(1:nrow(d.mod3),500,F) 
#  write.csv(d.mod3[ix,],paste("C:/Users/Tibs/Dropbox/sim/simdata/studentperf/lowdim_studentperf_cont_mod3_",i,".csv",sep=""),row.names = FALSE)
#}


#######################################################################
# Modification 4: Treatment effect heterogeneity + potential model misspecification
d.mod4 <- dat

# propensity score with manually specified coefficients
logit.g1W.mod4 <- 2*abs(age.cycle) + log(dat$sex/(dat$absences+1)+1) + dat$nursery*dat$failures - 3*risk.ft/risk.edu

g1W.mod4 <- (plogis(logit.g1W.mod4)+0.05)*0.90 # .05 to 0.92
#summary(g1W.mod4)
#0.04536 0.46760 0.59520 0.56170 0.68770 0.92220 
# Q model with manually specified coefficients
lp.drs.mod4 <- 4+2*abs(age.cycle) + exp(dat$sex/(dat$absences+1)+1) - 2*dat$nursery/(dat$failures+1) + risk.ft/risk.edu
#summary(lp.drs.mod4)

set.seed(seed12)
d.mod4$A <- rbinom(n, 1, g1W.mod4)
# transformation of age to [0,1] using logistic function
agetrans <- rank(dat$age)/100
# linear predictor for Y with age - treatment interaction 
d.mod4$Y <- rnorm(n,2*d.mod4$A +0.5*d.mod4$A * agetrans +lp.drs.mod4,sd(m.or$residuals)*0.8)
#summary(d.mod4$Y)
psi0.mod4 <- mean((2+0.5*agetrans+lp.drs.mod4) - lp.drs.mod4)
#psi0.mod4
# 3.625

# C statistic
#res.rocr.mod4 <- prediction(g1W.mod4, labels = d.mod4$A)
#performance(res.rocr.mod4, measure = "auc")@y.values[[1]]  # 0.71

names(d.mod4)[which(names(d.mod4)!="A"&names(d.mod4)!="Y")]<-paste("V",1:(length(names(d.mod4))-2),sep="")

#simcc.mod12 <- calcMetrics(d.mod4, Y1 = (2+0.5*agetrans+lp.drs.mod4), Y0 = lp.drs.mod4,
#                           ps.full = g1W.mod4,simName = "student_perf_mod4")

rownames(d.mod4)<-1:nrow(d.mod4)

DGP<-"lowdim_studentperf2_cont_mod4"
ATE<-psi0.mod4
filename<-paste("low",seed12,".csv",sep="")
gendat<-data.frame(DGP,filename,ATE)
write.csv(d.mod4,paste("C:/Users/Tibs/Dropbox/sim/simdat/low",seed12,".csv",sep=""),row.names = FALSE)
write.table(gendat,paste("C:/Users/Tibs/Dropbox/sim/gendatsummary",".csv",sep=""),append=T,row.names = FALSE,col.names = FALSE,sep=",")

#for(i in 301:400)
#{
#  set.seed(i)
#  ix<-sample(1:nrow(d.mod4),500,F) 
#  write.csv(d.mod4[ix,],paste("C:/Users/Tibs/Dropbox/sim/simdata/studentperf/lowdim_studentperf_cont_mod4_",i,".csv",sep=""),row.names = FALSE)
#}
}