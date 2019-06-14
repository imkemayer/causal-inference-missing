
set.seed(3847293)

############# first load the data
load("example.data")
source("functions.R")

covs.cont=c("bw","momage","nnhealth","birth.o","parity","moreprem","cigs","alcohol","ppvt.imp")
covs.cat=c("bwg","female","mlt.birt","b.marry","livwho","language","whenpren","drugs","othstudy","mom.lths","mom.hs","mom.coll","mom.scoll","site1","site2","site3","site4","site5","site6","site7","site8","momblack","momhisp","momwhite","workdur.imp")
covs=c(covs.cont,covs.cat)
ncovs=length(covs)

treat=ihdp$treat

#### FIRST OLS
use=ihdp[,c("iqsb.36","treat",covs)]
summary(lm(use))$coef[1:3,]

#               Estimate   Std. Error  t value     Pr(>|t|)
#(Intercept) 73.89299419 11.184169606 6.606927 6.808134e-11
#treat        9.18586845  1.046878812 8.774529 8.868121e-18
#bw           0.00823738  0.001829119 4.503468 7.590583e-06

mod=lm(ihdp[,c("iqsb.36","treat","ncdctt",covs)])
summary(mod)$coef[1:3,]
#(Intercept) 76.15030228 11.0281306 6.90509618 9.647532e-12
#treat        0.07225855  2.0259079 0.03566724 9.715558e-01
#ncdctt       3.18134984  0.6086712 5.22671338 2.159356e-07
plot(x=ihdp$ncdctt[!is.na(ihdp$iqsb.36)],y=mod$resid)

ihdp$ncdct2=ihdp$ncdctt^2

mod=lm(ihdp[,c("iqsb.36","treat","ncdctt","ncdct2",covs)])
summary(mod)$coef[1:4,]

#              Estimate Std. Error   t value     Pr(>|t|)
#(Intercept) 77.1011272 11.0720524 6.9635804 6.523248e-12
#treat        1.3115311  2.3956719 0.5474586 5.842036e-01
#ncdctt       0.9924712  2.3388263 0.4243458 6.714182e-01
#ncdct2       0.4983410  0.5141301 0.9692896 3.326692e-01
plot(x=ihdp$ncdctt[!is.na(ihdp$iqsb.36)],y=mod$resid)

## but to mimic BART perhaps we should just see if we can pick
## up the quadratic in the treatment group

summary(lm(ihdp[,c("iqsb.36","treat","ncdctt",covs)],subset=treat==1))$coef[1:3,]
#                Estimate   Std. Error  t value     Pr(>|t|)
#(Intercept) 52.619964931 17.923920485 2.935740 3.573176e-03
#ncdctt       3.309514633  0.634634578 5.214835 3.344955e-07
#bw           0.009482722  0.003006627 3.153940 1.766619e-03

summary(lm(ihdp[,c("iqsb.36","treat","ncdctt","ncdct2",covs)],subset=treat==1))$coef[1:4,]
#                Estimate   Std. Error   t value    Pr(>|t|)
#(Intercept) 54.901871282 18.365142095 2.9894607 0.003016209
#ncdctt       1.899274620  2.502071050 0.7590810 0.448375360
#ncdct2       0.321614983  0.551914200 0.5827264 0.560497050
#bw           0.009520583  0.003010495 3.1624640 0.001717772

### THEN SET-UP AND RUN BART

## first for iq at age 3
use=ihdp[,c("iqsb.36","treat","ncdctt",covs)]

xt=as.matrix(na.omit(use)[,-1])
xp1=as.matrix(use[use$treat==1,-1])
xp2=xp1
xp2[,c(1:2)]=0
xp=rbind(xp1,xp2)

nt=sum(use$treat==1)

y=as.numeric(na.omit(use)[,1])

library(BayesTree)
bart.tot <- bart(x.train=xt,   y.train=y,  x.test=xp)
save.image()

#### results
# first just effect of treatment on treated
tmpa=mean(bart.tot$yhat.test.mean[1:nt])-mean(bart.tot$yhat.test.mean[(nt+1):(2*nt)])
tmpa
# [1] 8.78
ndraws=nrow(bart.tot$yhat.test)
tmp=apply(bart.tot$yhat.test[,1:nt]-bart.tot$yhat.test[,(nt+1):(2*nt)],1,mean)
sd=sqrt(var(tmp)) # 
ci=c(tmpa-1.96*sd,tmpa+1.96*sd)
#[1]  6.8 10.7

pp.draws1 <- bart.tot$yhat.test[,1:nt]
pp.draws0 <- bart.tot$yhat.test[,(nt+1):(2*nt)]

# now let's plot out response surface as a function of ncdctt

ci.fun <- function(a){
c(quantile(a,.025),quantile(a,.975))
}

cis1 <- apply(pp.draws1,2,ci.fun)
cis0 <- apply(pp.draws0,2,ci.fun)

tmp=bart.tot$yhat.test[,1:nt]-bart.tot$yhat.test[,(nt+1):(2*nt)]
tes=apply(tmp,2,mean)
te.cis=apply(tmp,2,ci.fun)

postscript("ihdp.bart.ps", height=4, width=8.5, horizontal=T)
par(mfrow=c(1,2))
plot(lowess(xp1[,2],bart.tot$yhat.test.mean[1:nt]),pch=20,xlab="Number of CDC days (100)",   ylab="IQ at age 3",col="red",ylim=c(75,105))
points(lowess(xp1[,2],bart.tot$yhat.test.mean[(nt+1):(2*nt)]),pch=20)

plot(lowess(xp1[,2],tes),xlab="Number of CDC days (100)",ylab="Conditional treatment effects",ylim=c(-5,20))
points(lowess(xp1[,2],te.cis[1,]),type="l",lty=2)
points(lowess(xp1[,2],te.cis[2,]),type="l",lty=2)
abline(h=0)
dev.off()


