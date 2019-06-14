# probably best run as a batch job on a unix/linux machine

# Simulation TO see how BART compares with other methods over repeated samples
# AND varying coefficients -- this setting does not ensure overlap

niters=1000

library(BayesTree)
set.seed(2659232)

source("functions.R")

#  THIS First part doesn't vary over sims 


load(file="sim.data")
obs <- imp1[!(imp1$treat==1 & imp1$momwhite==0),]

covs.cont.n=c("bw","b.head","preterm","birth.o","nnhealth","momage")
covs.cat.n=c("sex","twin","b.marr","mom.lths","mom.hs",	"mom.scoll","cig","first","booze","drugs","work.dur","prenatal","ark","ein","har","mia","pen","tex","was")
p=length(c(covs.cont.n,covs.cat.n))

### calculate pscores and weights for tot
Trt=obs$treat
form.qx=as.formula(obs[,c("treat",covs.cont.n,covs.cat.n)])
qx=glm(data=obs[,c("treat",covs.cont.n,covs.cat.n)],formula=form.qx,family=binomial)$fitted
wts=rep(1,nrow(obs))
# now *controls* get a weight of 1
wts[Trt==1]=(1-qx[Trt==1])/(qx[Trt==1])
# treated get weights equal to the probability of being untreated over
# probability of being treated (to weight them up to look like controls)

#### get data in right format for BART
xt=obs[,c(covs.cont.n,covs.cat.n,"treat")]
xt=as.matrix(xt)
xp1=xt[xt[,"treat"]==0,]
xp2=xp1
xp1[,ncol(xt)]=1
xp2[,ncol(xt)]=0
xp=rbind(xp1,xp2)

#### probably should change this to "nc" rather than nt to be clear
#### but its a hassle to change everything below so we'll stick with this!
nt=sum(1-obs$treat)

##################### now simulate outcome data
##### covariate data, X
covs.ols = c(covs.cont.n,covs.cat.n)
X = obs[,covs.ols]
#X = na.omit(X)
# now standardize the continuous variables
X[,covs.cont.n]=as.data.frame(t((t(X[,covs.cont.n])-unlist(lapply(X[,covs.cont.n],mean)))/sqrt(unlist(lapply(X[covs.cont.n],var)))))
# record numbers of units and covariates
N = nrow(X)
dimx = ncol(X)
Xmat = as.matrix(X)

### now create matrix of all interactions etc for third response surface
ytmp=rnorm(N)
mod.bal <- glm(formula=ytmp~(bw+b.head+preterm+birth.o+nnhealth+momage+sex+twin+b.marr+mom.lths+mom.hs+mom.scoll+cig+first+booze+drugs+work.dur+prenatal+ark+ein+har+mia+pen+tex+was)^2 + I(bw^2) + I(b.head^2) + I(preterm^2) + I(birth.o^2) + I(nnhealth^2) + I(momage^2),x=T,data=cbind.data.frame(Xmat))
coefs <- mod.bal$coef[-1]
XX <- mod.bal$x[,-1]
XX <- XX[,!is.na(coefs)]

nouts=3
os=c("YA","YB","YC")

### we'll save t.e. estimates, whether the interval covered, and
#length of interval
results.a=matrix(0,niters,18)
results.b=matrix(0,niters,18)
results.c=matrix(0,niters,18)
dimnames(results.a)=list(NULL,c("b.te","b.cov","b.cil","r.te","r.cov","r.cil","ps.te","ps.cov","ps.cil","ipw.te","ipw.cov","ipw.cil","tau.est","r2","b.wrong","r.wrong","ps.wrong","ipw.wrong"))
dimnames(results.b)=list(NULL,c("b.te","b.cov","b.cil","r.te","r.cov","r.cil","ps.te","ps.cov","ps.cil","ipw.te","ipw.cov","ipw.cil","tau.est","r2","b.wrong","r.wrong","ps.wrong","ipw.wrong"))
dimnames(results.c)=list(NULL,c("b.te","b.cov","b.cil","r.te","r.cov","r.cil","ps.te","ps.cov","ps.cil","ipw.te","ipw.cov","ipw.cil","tau.est","r2","b.wrong","r.wrong","ps.wrong","ipw.wrong"))

precision.a=matrix(0,niters,4)
dimnames(precision.a)=list(NULL,c("bart","reg","psm","ipw"))
precision.b=matrix(0,niters,4)
dimnames(precision.b)=list(NULL,c("bart","reg","psm","ipw"))
precision.c=matrix(0,niters,4)
dimnames(precision.c)=list(NULL,c("bart","reg","psm","ipw"))

rm(imp1)
save.image()

XXXmat=cbind(rep(1,N),XX)
rm(XX)

for(i in 1:niters){
if(i<=500){set.seed(565 + i*5)}
if(i>500){set.seed(7565 + i*5)}

### here's where things start to vary

############## RESPONSE SURFACES (3 versions)

# YA1 and YA0 are response surfaces corresponding to assignment to treatment
# and assignment to control, respectively;  these are both linear
sigy = 1
tau = 4
betaA = sample(c(0:4),dimx+1,replace=TRUE,prob=c(.5,.2,.15,.1,.05))
yahat = cbind(rep(1, N), Xmat) %*% betaA
YA0 = rnorm(N, yahat, sigy)
YA1 = rnorm(N, yahat+tau, sigy)
# YA is the vector of observed responses
YA = YA1; YA[Trt==0] = YA0[Trt==0]
tauAis=4
#
# YB1 and YB0 are response surfaces corresponding to assignment to treatment
# and assignment to control, respectively;  the former is non-linear
betaB = c(sample(c(.0,.1,.2,.3,.4),(dimx+1),replace=TRUE,prob=c(.6,.1,.1,.1,.1))) 
yb0hat = exp((cbind(rep(1, N), (Xmat+.5)) %*% betaB))
yb1hat = cbind(rep(1, N), (Xmat+.5)) %*% betaB 
offset = c(mean(yb1hat[Trt==0] - yb0hat[Trt==0])) - 4
yb1hat = cbind(rep(1, N), (Xmat+.5)) %*% betaB -offset
YB0 = rnorm(N, yb0hat, sigy)
YB1 = rnorm(N, yb1hat, sigy)
# try 1 set sigy to 2 
# YB is the vector of observed responses
YB = YB1; YB[Trt==0] = YB0[Trt==0]
tauBis = yb1hat[Trt==0] - yb0hat[Trt==0]
tauB = mean(tauBis)
#
# YC1 and YC0 are response surfaces corresponding to assignment to treatment
# and assignment to control, respectively;  these are both non-linear in X (lots
#  of interactions) and now non-parallel as well
sigy = 1
tau = 4
#
# main effects coefficients
betaC.m0 = sample(c(0,1,2),p+1,replace=T,prob=c(.6,.3,.1))
betaC.m1 = sample(c(0,1,2),p+1,replace=T,prob=c(.6,.3,.1))
# quadratic coefficients
#these we make pretty rare since they really represent 3-way interactions
betaC.q0 = sample(c(0,.5,1),ncol(XXXmat)-(p+1),replace=TRUE,prob=c(.8,.15,.05))
betaC.q1 = sample(c(0,.5,1),ncol(XXXmat)-(p+1),replace=TRUE,prob=c(.8,.15,.05))
#
betaC0 = c(betaC.m0,betaC.q0)
betaC1 = c(betaC.m1,betaC.q1)
yc0hat = (XXXmat) %*% betaC0
yc1hat = (XXXmat) %*% betaC1 
offset = c(mean(yc1hat[Trt==0] - yc0hat[Trt==0])) - 4
yc1hat = (XXXmat) %*% betaC1 - offset
YC0 = rnorm(N, yc0hat, sigy)
YC1 = rnorm(N, yc1hat, sigy)
# YC is the vector of observed responses
YC = YC1; YC[Trt==0] = YC0[Trt==0]
tauCis = yc1hat[Trt==0] - yc0hat[Trt==0]
tauC = mean(tauCis)
#
# generate true individual level
# generate sample treatment effects 
tauAs = mean(YA1[Trt==0] - YA0[Trt==0])
tauBs = mean(YB1[Trt==0] - YB0[Trt==0])
tauCs = mean(YC1[Trt==0] - YC0[Trt==0])

taus2 = c(tauAs,tauBs,tauCs)

rm(betaC0,betaC1,betaC.m0,betaC.m1,betaC.q0,betaC.q1,yc0hat,yc1hat,yb0hat,yb1hat)

results.a[i,13]=taus2[1]
results.b[i,13]=taus2[2]
results.c[i,13]=taus2[3]
  
######################
library(BayesTree)
bart2a <- bart(x.train=xt,   y.train=YA,  x.test=xp)
bart2b <- bart(x.train=xt,   y.train=YB,  x.test=xp)
bart2c <- bart(x.train=xt,   y.train=YC,  x.test=xp)
save.image()

########## results

# a
tmp=apply(bart2a$yhat.test[,1:nt]-bart2a$yhat.test[,(nt+1):(2*nt)],1,mean)
tmpa=mean(tmp)
results.a[i,1]=4-tmpa
#
ndraws=nrow(bart2a$yhat.test)
sd=sqrt(var(tmp))
ci=c(tmpa-1.96*sd,tmpa+1.96*sd)
results.a[i,2]=(ci[1]<4 & ci[2]>4)*1
results.a[i,3]=ci[2]-ci[1]
results.a[i,15]=(ci[1]<0 & ci[2]>0)*1

tmp2=apply(bart2a$yhat.test[,1:nt]-bart2a$yhat.test[,(nt+1):(2*nt)],2,mean)
precision.a[i,1] = sqrt(mean((tmp2-tauAis)^2))
rm(tmp,tmp2)

# b
tmp=apply(bart2b$yhat.test[,1:nt]-bart2b$yhat.test[,(nt+1):(2*nt)],1,mean)
tmpb=mean(tmp)
results.b[i,1]=4-tmpb
#
ndraws=nrow(bart2b$yhat.test)
sd=sqrt(var(tmp))
ci=c(tmpb-1.96*sd,tmpb+1.96*sd)
results.b[i,2]=(ci[1]<4 & ci[2]>4)*1
results.b[i,3]=ci[2]-ci[1]
results.b[i,15]=(ci[1]<0 & ci[2]>0)*1

tmp2=apply(bart2b$yhat.test[,1:nt]-bart2b$yhat.test[,(nt+1):(2*nt)],2,mean)
precision.b[i,1] = sqrt(mean((tmp2-tauBis)^2))
rm(tmp,tmp2)

# c
tmp=apply(bart2c$yhat.test[,1:nt]-bart2c$yhat.test[,(nt+1):(2*nt)],1,mean)
tmpc=mean(tmp)
results.c[i,1]=4-tmpc
#
ndraws=nrow(bart2c$yhat.test)
sd=sqrt(var(tmp))
ci=c(tmpc-1.96*sd,tmpc+1.96*sd)
results.c[i,2]=(ci[1]<4 & ci[2]>4)*1
results.c[i,3]=ci[2]-ci[1]
results.c[i,15]=(ci[1]<0 & ci[2]>0)*1

tmp2=apply(bart2c$yhat.test[,1:nt]-bart2c$yhat.test[,(nt+1):(2*nt)],2,mean)
precision.c[i,1] = sqrt(mean((tmp2-tauCis)^2))
rm(tmp,tmp2)

##### how did competitors do?
data2 = cbind.data.frame(X,Trt=Trt,YA=YA,YB=YB,YC)
tmpp=summary(lm(data2[,c("YA","Trt",covs.cat.n,covs.cont.n)]))
tmp=tmpp$coef[2,1:2]
results.a[i,4]=4-tmp[1]
ci=c(tmp[1]-1.96*tmp[2],tmp[1]+1.96*tmp[2])
results.a[i,5]=(ci[1]<4 & ci[2]>4)*1
results.a[i,6]=ci[2]-ci[1]
results.a[i,14]=tmpp$r.squared
results.a[i,16]=(ci[1]<0 & ci[2]>0)*1
precision.a[i,2] = sqrt(mean((tmp[1]-tauAis)^2))

tmpp=summary(lm(data2[,c("YB","Trt",covs.cat.n,covs.cont.n)]))
tmp=tmpp$coef[2,1:2]
results.b[i,4]=4-tmp[1]
ci=c(tmp[1]-1.96*tmp[2],tmp[1]+1.96*tmp[2])
results.b[i,5]=(ci[1]<4 & ci[2]>4)*1
results.b[i,6]=ci[2]-ci[1]
results.b[i,14]=tmpp$r.squared
results.b[i,16]=(ci[1]<0 & ci[2]>0)*1
precision.b[i,2] = sqrt(mean((tmp[1]-tauBis)^2))

tmpp=summary(lm(data2[,c("YC","Trt",covs.cat.n,covs.cont.n)]))
tmp=tmpp$coef[2,1:2]
results.c[i,4]=4-tmp[1]
ci=c(tmp[1]-1.96*tmp[2],tmp[1]+1.96*tmp[2])
results.c[i,5]=(ci[1]<4 & ci[2]>4)*1
results.c[i,6]=ci[2]-ci[1]
results.c[i,14]=tmpp$r.squared
results.c[i,16]=(ci[1]<0 & ci[2]>0)*1
precision.c[i,2] = sqrt(mean((tmp[1]-tauCis)^2))

#### now also compare to p-score matching (now with robust standard errors)
tmp=matrix(0,3,2)
tmp[1,]=pscores.fun(treat=(1-Trt),outs=YA,covs=as.matrix(X))
tmp[2,]=pscores.fun(treat=(1-Trt),outs=YB,covs=as.matrix(X))
tmp[3,]=pscores.fun(treat=(1-Trt),outs=YC,covs=as.matrix(X))

results.a[i,7]=4-(-tmp[1,1])
ci=c(-tmp[1,1]-1.96*tmp[1,2],-tmp[1,1]+1.96*tmp[1,2])
results.a[i,8]=(ci[1]<4 & ci[2]>4)*1
results.a[i,9]=ci[2]-ci[1]
results.a[i,17]=(ci[1]<0 & ci[2]>0)*1
precision.a[i,3] = sqrt(mean((-tmp[1,1]-tauAis)^2))

results.b[i,7]=4-(-tmp[2,1])
ci=c(-tmp[2,1]-1.96*tmp[2,2],-tmp[2,1]+1.96*tmp[2,2])
results.b[i,8]=(ci[1]<4 & ci[2]>4)*1
results.b[i,9]=ci[2]-ci[1]
results.b[i,17]=(ci[1]<0 & ci[2]>0)*1
precision.b[i,3] = sqrt(mean((-tmp[1,1]-tauBis)^2))

results.c[i,7]=4-(-tmp[3,1])
ci=c(-tmp[3,1]-1.96*tmp[3,2],-tmp[3,1]+1.96*tmp[3,2])
results.c[i,8]=(ci[1]<4 & ci[2]>4)*1
results.c[i,9]=ci[2]-ci[1]
results.c[i,17]=(ci[1]<0 & ci[2]>0)*1
precision.c[i,3] = sqrt(mean((-tmp[1,1]-tauCis)^2))

#### now also compare to inverse probability weighting
#### need to use wls.all2 to get robust standard errors

tmpp=wls.all2(X=cbind(rep(1,N),Trt,as.matrix(X)),w=wts,Y=YA,treat=Trt)
tmp=c(tmpp[3],sqrt(tmpp[2]))
results.a[i,10]=4-tmp[1]
ci=c(tmp[1]-1.96*tmp[2],tmp[1]+1.96*tmp[2])
results.a[i,11]=(ci[1]<4 & ci[2]>4)*1
results.a[i,12]=ci[2]-ci[1]
results.a[i,18]=(ci[1]<0 & ci[2]>0)*1
precision.a[i,4] = sqrt(mean((tmp[1]-tauAis)^2))

tmpp=wls.all2(X=cbind(rep(1,N),Trt,as.matrix(X)),w=wts,Y=YB,treat=Trt)
tmp=c(tmpp[3],sqrt(tmpp[2]))
results.b[i,10]=4-tmp[1]
ci=c(tmp[1]-1.96*tmp[2],tmp[1]+1.96*tmp[2])
results.b[i,11]=(ci[1]<4 & ci[2]>4)*1
results.b[i,12]=ci[2]-ci[1]
results.b[i,18]=(ci[1]<0 & ci[2]>0)*1
precision.b[i,4] = sqrt(mean((tmp[1]-tauBis)^2))

tmpp=wls.all2(X=cbind(rep(1,N),Trt,as.matrix(X)),w=wts,Y=YC,treat=Trt)
tmp=c(tmpp[3],sqrt(tmpp[2]))
results.c[i,10]=4-tmp[1]
ci=c(tmp[1]-1.96*tmp[2],tmp[1]+1.96*tmp[2])
results.c[i,11]=(ci[1]<4 & ci[2]>4)*1
results.c[i,12]=ci[2]-ci[1]
results.c[i,18]=(ci[1]<0 & ci[2]>0)*1
precision.c[i,4] = sqrt(mean((tmp[1]-tauCis)^2))
}
save.image()
