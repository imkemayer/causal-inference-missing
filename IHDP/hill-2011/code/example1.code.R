### try fitting various matching methods to the data and compare to bart
### save treatment effect estimates and s.e.'s (approx for matching) for all
### save balance summaries for unmatched and all matched

set.seed(3847293)
library(MatchIt)

############# first load the data
load("example.data")
source("functions.R")

covs.cont=c("bw","momage","nnhealth","birth.o","parity","moreprem","cigs","alcohol","ppvt.imp")
covs.catF=c("bwg","female","mlt.birtF","b.marryF","livwhoF","languageF","whenprenF","drugs","othstudy","momed4F","siteF","momraceF","workdur.imp")
covsF=c(covs.cont,covs.catF)
ncovs.cont=length(covs.cont)

usek = na.omit(ihdp[!(ihdp$treat==1 & ihdp$dose400==0),c("iqsb.36","dose400",covsF)])

########################  run methods, record balance, tes, ses  ########################

formyF = as.formula(usek[,c("iqsb.36","dose400",covsF)])
formy0 = as.formula(usek[,c(1:2)])

# estimate the pscore with an additive model
formzF = as.formula(usek[,c("dose400",covsF)])
modqx = glm(formzF,data=usek,family="binomial",x=TRUE)
qx = modqx$fitted
##  now figure total # covs out from the design matrix because of the complications caused by the factors
ncovs=ncol(modqx$x)-1 	
 	
library(mgcv)
# took the smooths off of parity, birth.o and moreprem because # unique cats was so small
# i was getting error messages
form.gamF = as.formula("dose400 ~ s(bw) + s(momage) + s(nnhealth) + birth.o + parity + moreprem + s(cigs) + s(alcohol) + s(ppvt.imp) + bwg + female + mlt.birtF + b.marryF + livwhoF + languageF + whenprenF + drugs + othstudy + momed4F + siteF + momraceF + workdur.imp")

## now pscore estimated using GAM (easier to do outside of matchit)
qxg = gam(form.gamF,data=usek,family="binomial")$fitted

### now create pscore model to try to emphasize the terms i'm going to be judging balance on
form.quadz <- as.formula("dose400 ~ (bw + momage + nnhealth + birth.o + parity + moreprem + cigs + alcohol + ppvt.imp + bwg + female + mlt.birtF + b.marryF + livwhoF + languageF + whenprenF + drugs + othstudy + momed4F + siteF + momraceF + workdur.imp)^2 + I(bw^2) + I(momage^2) + I(nnhealth^2) + I(birth.o^2) + I(parity^2) + I(moreprem^2) + I(cigs^2) + I(alcohol^2) + I(ppvt.imp^2)")
modq0=glm(formula=form.quadz,data=usek,family="binomial",x=TRUE)
qxq=modq0$fitted

# pulling out most important terms from this regression and then adding quadratic
# terms only for those
form.red <- as.formula("dose400 ~ (bw + nnhealth + parity + moreprem + ppvt.imp + bwg + b.marryF + languageF + momed4F + siteF + momraceF)")
mod.red=glm(formula=form.red,data=usek,family="binomial",x=TRUE)
qx.red=mod.red$fitted
# dimnames(mod.red$x)[[2]]
rm(mod.red)

### now create model to set up all the terms i'm going to be judging balance on
form.bal <- as.formula("dose400 ~ (bw + momage + nnhealth + birth.o + parity + moreprem + cigs + alcohol + ppvt.imp + bwg + female + mlt.birtF + b.marryF + livwhoF + languageF + whenprenF + drugs + othstudy + momed4F + siteF + momraceF + workdur.imp)^2")
# balance in dist for continuous will be examined using QQ stats rather than squared terms

## specify the combination of options for matching
methods = c(rep("nearest",24),rep("optimal",12),rep("full",4),rep("iptw",4))
ratios  = c(rep(c(1,2,3),8),rep(c(1,2,3),4),c(1,1,1))
which.link = c(rep(1,3),rep(2,3),rep(3,3),rep(4,3),rep(1,3),rep(2,3),rep(3,3),rep(4,3),rep(1,3),rep(2,3),rep(3,3),rep(4,3),c(1,2,3,4))
links = list(qx=qx,quad=qxq,gam=qxg,reduced=qx.red)
# to differentiate in the plot
# nearest: 
labels = c(rep("N",24),rep("O",12),rep("F",4),rep("W",4))

n.methods = length(methods)+1
tes.400.dm = matrix(0,n.methods,2)
tes.400.ols = matrix(0,n.methods,2)
balb.400.mat = matrix(0,n.methods,3)
balc1.400.mat = matrix(0,n.methods,7)
balc2.400.mat = matrix(0,n.methods,7)
nms = c("no match",rep(NA,n.methods-1))

num.cont=9
num.cov=41
nc=sum(usek$dose400==0)

### FIRST OLS
tes.400.dm[1, ] = summary(lm(formula=formy0, data=usek))$coef[2,1:2]
tes.400.ols[1, ] = summary(lm(formula=formyF, data=usek))$coef[2,1:2]

#######################################################
#############  now all the matching craziness

# default for nearest neighbor matching is no replacement
# only way to run optimal and full matching is with no replacement

###  FIGURE OUT GAM OPTION

library(MatchIt)
## here we can iterate through the following options
for(k in 1:length(methods)){
#
if(methods[k]=="iptw"){
ps=links[[k-40]]
wts = rep(1,nrow(usek))
wts[usek$dose400==0] = ps[usek$dose400==0]/(1-ps[usek$dose400==0])
wts[usek$dose400==0] = wts[usek$dose400==0]*nc/sum(wts[usek$dose400==0])
# to trick matchit into calculating the balance statistics for us
# we first create matchit output for matching with appropriate
# distance measure then replace their weights with ours
m.out = matchit(form=form.bal, data=usek, distance=ps)
m.out$weights=wts
bal = balance.sum(m.out,num.cont=num.cont,num.cov=num.cov)
balb.400.mat[k+1,] = unlist(bal$bal.bin)
balc1.400.mat[k+1,] = unlist(bal$bal.cont1)
balc2.400.mat[k+1,] = unlist(bal$bal.cont2)
tes.400.dm[k+1,]  = summary(lm(formula=formy0,data=usek,weights=wts))$coef[2,1:2]
tes.400.ols[k+1,] = summary(lm(formula=formyF,data=usek,weights=wts))$coef[2,1:2]
cat(k,tes.400.ols[k+1,],"\n")
nms[k+1] = paste(methods[k],k-40,sep=".")
}
else{
if(methods[k]=="full"){
m.out = matchit(form=form.bal, data=usek, method=methods[k], distance=links[[which.link[k]]])
nms[k+1] = paste(methods[k],ratios[k],which.link[k],sep=".")
}
if(methods[k]!="full"){
if(k>0 & k<13){
m.out = matchit(form=form.bal, data=usek, method=methods[k], ratio=ratios[k], distance=links[[which.link[k]]],replace=TRUE)
nms[k+1] = paste(methods[k],ratios[k],which.link[k],"R",sep=".")
}
else{
m.out = matchit(form=form.bal, data=usek, method=methods[k], ratio=ratios[k], distance=links[[which.link[k]]])
nms[k+1] = paste(methods[k],ratios[k],which.link[k],sep=".")
}
}
if(k==1){
bal = balance.sum(m.out,num.cont=num.cont,num.cov=num.cov,matched=FALSE)
balb.400.mat[1,] = unlist(bal$bal.bin)
balc1.400.mat[1,] = unlist(bal$bal.cont1)
balc2.400.mat[1,] = unlist(bal$bal.cont2)
}
bal = balance.sum(m.out,num.cont=num.cont,num.cov=num.cov)
balb.400.mat[k+1,] = unlist(bal$bal.bin)
balc1.400.mat[k+1,] = unlist(bal$bal.cont1)
balc2.400.mat[k+1,] = unlist(bal$bal.cont2)
mdat = match.data(m.out, weights="wts")
tes.400.dm[k+1,]  = summary(lm(formula=formy0,data=mdat,weights=wts))$coef[2,1:2]
tes.400.ols[k+1,] = summary(lm(formula=formyF,data=mdat,weights=wts))$coef[2,1:2]
cat(k,tes.400.ols[k+1,],"\n")
}
}

cbind.data.frame(nms,round(balb.400.mat,2))
cbind.data.frame(nms,round(balc1.400.mat,2))
cbind.data.frame(nms,round(balc2.400.mat,2))

# look only at the treatment effect estimates with the best balance
# on the univariate statistics


## plot of all
## a small amount of jitter is added to make it easier to distinguish multiple
## observations at the same point
ind = balb.400.mat[,2]<.2 & balc1.400.mat[,1]<.2
res=tes.400.ols[ind,1]
plot(y=rep(1,length(res)),x=res,ylim=c(.5,6.5),xlab="treatment effect estimates",ylab="",yaxt="n",xlim=c(6.5,max(res)),mgp=c(2,.5,0),cex=1,pch=labels[ind])
axis(side=2,at=c(1,2,3,4,5),labels=c("Row 1","Row 2","Row 3","Row 4","Row 5"),las=1,mgp=c(4,.5,0),xlab="balance criteria")
text(x=c(6.6,7.6,8.6,9.6),y=rep(6.4,4),labels=c("max","max", "max","pct"),cex=.8)
text(x=c(6.6,7.6,8.6,9.6),y=rep(6,4),labels=c("STD(B)","STD(C)", "EQQmax(C)","STD(C2)>.1"),cex=.8)
text(x=c(6.6,7.6),y=rep(1,2),labels=c(".2",".2"),cex=.8)
#
ind = balb.400.mat[,2]<.15 & balc1.400.mat[,2]<.15
text(x=c(6.6,7.6),y=rep(2,2),labels=c(".15",".15"),cex=.8)
res=tes.400.ols[ind,1]
points(y=rep(2,sum(ind)),x=res,ylim=c(0,4),cex=1,pch=labels[ind])
ind = balb.400.mat[,2]<.13 & balc1.400.mat[,2]<.13 & balc1.400.mat[,7]<.15
text(x=c(6.6,7.6,8.6),y=rep(3,3),labels=c(".13",".13", ".15"),cex=.8)
res=tes.400.ols[ind,1]
points(y=rep(3,sum(ind)),x=res,ylim=c(0,4),cex=1,pch=labels[ind])
ind = balb.400.mat[,2]<.13 & balc1.400.mat[,2]<.13 & balc1.400.mat[,7]<.15 & balc2.400.mat[,3]/742<.2
text(x=c(6.6,7.6,8.6,9.6),y=rep(4,4),labels=c(".13",".13",".15",".2"),cex=.8)
res=tes.400.ols[ind,1]
points(y=rep(4,sum(ind)),x=res,ylim=c(0,4),cex=1,pch=labels[ind])
ind = balb.400.mat[,2]<.11 & balc1.400.mat[,2]<.11 & balc1.400.mat[,7]<.11 & balc2.400.mat[,3]/742<.2
text(x=c(6.6,7.6,8.6,9.6),y=rep(5,4),labels=c(".11",".11",".11",".2"),cex=.8)
res=tes.400.ols[ind,1]
points(y=rep(5,sum(ind)),x=res,ylim=c(0,4),cex=1,pch=labels[ind])

## add in line for BART estimate
abline(v=12.9,col="red")

## add in the genmatch points (for these results see example1.code.gm.R)
points(y=c(1,1),x=c(13.2,13.6),cex=1,pch="G")
points(y=c(2,2),x=c(13.2,13.6),cex=1,pch="G")
points(y=c(3,3),x=c(13.2,13.6),cex=1,pch="G")
points(y=c(4,4),x=c(13.2,13.6),cex=1,pch="G")
points(y=c(5),x=c(13.2),cex=1,pch="G")
## 13.2 for intelligent, 13.6 for default

###########################################################################
## now BART

covs.cont=c("bw","momage","nnhealth","birth.o","parity","moreprem","cigs","alcohol","ppvt.imp")
covs.cat=c("bwg","female","mlt.birt","b.marry","livwho","language","whenpren","drugs","othstudy","mom.lths","mom.hs","mom.coll","mom.scoll","site1","site2","site3","site4","site5","site6","site7","site8","momblack","momhisp","momwhite","workdur.imp")
covs=c(covs.cont,covs.cat)
ncovs=length(covs)

usek = na.omit(ihdp[!(ihdp$treat==1 & ihdp$dose400==0),c("iqsb.36","dose400",covs)])

# important to excluse the treated with "low" dose since "highdose=0"
# for them means something very different than it does for those in
# then control group

xt=as.matrix(usek[,-1])
xp=as.matrix(usek[usek$dose400==1,-1])
xp[,1]=0

y=as.numeric(usek[,1])

library(BayesTree)
bart.tot <- bart(x.train=xt,   y.train=y,  x.test=xp)
save.image()

# check convergence
plot(bart.tot$sigma)

#### results
# first just effect of treatment on treated
diffs=bart.tot$yhat.train[,usek$dose400==1]-bart.tot$yhat.test
mndiffs=apply(diffs,1,mean)
mean(mndiffs)
# 12.9

sd(mndiffs)
# 1.96

# get a sense of t.e. heterogeneity
hist(apply(diffs,2,mean))

