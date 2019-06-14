## Genetic Matching
## (some duplication in this file because these were originally
##  separate files, combined to reduce total number of files)
##############################################################################
### first the "default" method

set.seed(3847293)
library(Matching)

############# first load the data
load("example.data")
source("functions.R")

covs.cont=c("bw","momage","nnhealth","birth.o","parity","moreprem","cigs","alcohol","ppvt.imp")
covs.catF=c("bwg","female","mlt.birtF","b.marryF","livwhoF","languageF","whenprenF","drugs","othstudy","momed4F","siteF","momraceF","workdur.imp")
covsF=c(covs.cont,covs.catF)
# 1 will be added to each of these in the balance summary function to account for qx
ncovs.cont=length(covs.cont)
ncovs=40

usek = na.omit(ihdp[!(ihdp$treat==1 & ihdp$dose400==0),c("iqsb.36","dose400",covsF)])

########################  set-up  ########################
formzF = as.formula(usek[,c("dose400",covsF)])
modqx = glm(formzF,data=usek,family="binomial",x=TRUE)
qx = modqx$fitted

## genetic matching should also explicitly control for qx
usek2=cbind.data.frame(usek[,c("iqsb.36","dose400",covsF)],qx=qx)
formzq = as.formula(usek2[,c("dose400","qx",covsF)])
## use this just to get the design matrix needed for GenMatch
modqx2 = glm(formzq,data=usek,family="binomial",x=TRUE)

## need formula for balance statistics to explicitly control what quadratic
## terms are included so we don't waste time with squared binary variables
form.quadz <- as.formula("dose400 ~ qx + (bw + momage + nnhealth + birth.o + parity + moreprem + cigs + alcohol + ppvt.imp + bwg + female + mlt.birtF + b.marryF + livwhoF + languageF + whenprenF + drugs + othstudy + momed4F + siteF + momraceF + workdur.imp)^2")
mod.quad=glm(formula=form.quadz,data=usek2,x=TRUE)
ncol(mod.quad$x)
# should be (ncovs+1) =41 main effects; 
# 742 interactions
# so 783 total
# design matrix has 783 cols (ignoring constant term)

# need to create a version of the formula that will work in the function
form.quadz2 = formula(cbind.data.frame(dose400=usek2$dose400,data.frame(mod.quad$x[,-1])))

#########  now for genetic matching
library(Matching)
	
# deleted the oldstarting values since we've changed the data
##sv from run1d.R

##The following is needed for cluster stuff. Uncomment all to get a cluster going
#source("~/snow/snow1a.R") #this is a complicated cluster version of
##library(snow)
##setDefaultClusterOptions(master="localhost")
##cl <- makeSOCKcluster(c("localhost","localhost"))

#clusterExport(cl, "ncovs.cont") 
#clusterExport(cl, "ncovs")      
#clusterExport(cl, "form.quadz2")  
#clusterExport(cl, "fix")

#from hill.default.new1.ps10000.Rout, gen 24
sv <- c(2.479897e+02, 3.656868e+01, 4.920903e+02, 4.621738e+02,
        3.500506e+02, 2.839971e+02, 1.943023e+02, 3.391185e+02, 8.441648e+00,
        3.202975e+02, 1.627430e+02, 3.615579e+02, 7.967758e+02, 8.759419e+00,
        2.566829e+02, 7.492769e+02, 7.015629e+02, 5.379703e+02, 4.629878e+02,
        8.281951e+02, 3.194794e+02, 7.710584e+02, 8.747445e+00, 6.982671e+02,
        2.579425e+02, 4.841003e+01, 6.237803e+02, 3.018265e+02, 4.949530e+02,
        9.357231e+02, 7.023428e+02, 5.413681e+02, 5.202681e+02, 2.207877e+01,
        7.122965e+02, 2.508598e+02, 1.521980e+02, 9.563593e+01, 8.290899e+02,
        7.054497e+02, 5.171170e+02)

gm2 <- GenMatch(usek2$dose400,X=modqx2$x[,2:42],M=1,pop.size=1, 
                hard.generation.limit=TRUE,
                max.generations= 1,  
                wait.generations= 0,
#               cluster=cl,                                
                starting.values=sv)

cat("With BiasAdjustment\n")
mm2=Match(usek2[,1], Tr=usek2$dose400, X=modqx2$x[,2:42], Z=modqx2$x[,3:42], Weight.matrix=gm2, replace=TRUE, M=1,BiasAdjust=TRUE)
summary(mm2)
# 13.6 (2.53)

cat("Without BiasAdjustment\n")
mm2=Match(usek2[,1], Tr=usek2$dose400, X=modqx2$x[,2:42], Z=modqx2$x[,3:42], Weight.matrix=gm2, replace=TRUE, M=1,BiasAdjust=FALSE)
summary(mm2)

cat("Print Jennifer's Balance Stats\n")
bal.gm2 <- balance.gm3.sum(mat.out=mm2,matched=TRUE,dat=cbind.data.frame(dose400=usek2$dose400,data.frame(mod.quad$x)),
                          num.cont=ncovs.cont,num.cov=ncovs,form.q=form.quadz2)
print(bal.gm2)

##############################################################################
### now the "intelligent" method

### try fitting various matching methods to the data and compare to bart
### save treatment effect estimates and s.e.'s (approx for matching) for all
### save balance summaries for unmatched and all matched

set.seed(3847293)
library(Matching)

############# first load the data
load("example.data")
source("functions.R")

covs.cont=c("bw","momage","nnhealth","birth.o","parity","moreprem","cigs","alcohol","ppvt.imp")
covs.catF=c("bwg","female","mlt.birtF","b.marryF","livwhoF","languageF","whenprenF","drugs","othstudy","momed4F","siteF","momraceF","workdur.imp")
covsF=c(covs.cont,covs.catF)
# 1 will be added to each of these in the balance summary function to account for qx
ncovs.cont=length(covs.cont)
ncovs=40

usek = na.omit(ihdp[!(ihdp$treat==1 & ihdp$dose400==0),c("iqsb.36","dose400",covsF)])

########################  set-up  ########################
formzF = as.formula(usek[,c("dose400",covsF)])
modqx = glm(formzF,data=usek,family="binomial",x=TRUE)
qx = modqx$fitted

## genetic matching should also explicitly control for qx
usek2=cbind.data.frame(usek[,c("iqsb.36","dose400",covsF)],qx=qx)
formzq = as.formula(usek2[,c("dose400","qx",covsF)])
## use this just to get the design matrix needed for GenMatch
modqx2 = glm(formzq,data=usek,family="binomial",x=TRUE)

## need formula for balance statistics to explicitly control what quadratic
## terms are included so we don't waste time with squared binary variables
form.quadz <- as.formula("dose400 ~ qx + (bw + momage + nnhealth + birth.o + parity + moreprem + cigs + alcohol + ppvt.imp + bwg + female + mlt.birtF + b.marryF + livwhoF + languageF + whenprenF + drugs + othstudy + momed4F + siteF + momraceF + workdur.imp)^2")
mod.quad=glm(formula=form.quadz,data=usek2,x=TRUE)
ncol(mod.quad$x)
# should be (ncovs+1) =41 main effects; 
# 742 interactions
# so 783 total
# design matrix has 783 cols (ignoring constant term)

# need to create a version of the formula that will work in the function
form.quadz2 = formula(cbind.data.frame(dose400=usek2$dose400,data.frame(mod.quad$x[,-1])))

#########  now for genetic matching
library(Matching)

# deleted the oldstarting values since we've changed the data
##sv from run1d.R

myfit <- function(in.matches, BM)
  {
    # note: "matches" has three columns:
    # column 1: index of treated obs
    # column 2: index of control obs
    # column 3: weights for matched-pairs

#    form.quadz2 = formula(cbind.data.frame(dose400=usek2$dose400,data.frame(mod.quad$x[,-1])))
    
    matches <- list()
    matches$index.treated <- in.matches[,1]
    matches$index.control <- in.matches[,2]
    matches$weightsl      <- in.matches[,3]

    class(matches) <- "Match"

    rb <- balance.gm3.sum(mat.out=matches,matched=TRUE,dat=BM,
                         num.cont=ncovs.cont,
                         num.cov=ncovs,
                         form.q=form.quadz2,
                         verbose=FALSE,
                         xdata=in.xdata,
                         Tr=in.Tr)

    #order the values. Note we are minimizing!
    
    #reorder to counts are first
    total.over <- rb$bal.bin$bin.std.over.1+rb$bal.cont1$cont.std.over.1
    ret <- c(total.over,
             rb$bal.bin$bin.std.over.1,
             rb$bal.cont1$cont.std.over.1,
             rb$bal.cont1$cont.maxQQ.max,
             rb$bal.cont2$cont2.std.over.1,

             rb$bal.bin$bin.std.max,     
             rb$bal.cont1$cont.std.max,
             rb$bal.cont1$cont.maxQQ.mn,
             rb$bal.cont2$cont2.std.max,
             rb$bal.cont2$cont2.maxQQ.max,
             
             rb$bal.cont1$cont.medQQ.max,           
             rb$bal.cont2$cont2.medQQ.max,
             
             rb$bal.bin$bin.std.mn,
             rb$bal.cont1$cont.std.mn,
             rb$bal.cont1$cont.medQQ.mn,
             rb$bal.cont2$cont2.std.mn,
             rb$bal.cont2$cont2.medQQ.mn,
             rb$bal.cont2$cont2.maxQQ.mn
             )
             
  }#end of myfit

datafr <- cbind.data.frame(dose400=usek2$dose400,data.frame(mod.quad$x))
formul <- form.quadz2
in.xdata  <- as.data.frame(model.matrix(formul, data=datafr, drop.unused.levels = TRUE))
in.Tr  <- as.real(model.response(model.frame(formul, data=datafr)))

#if we were using a cluster setup we would need the following

##The following is needed for cluster stuff. Uncomment all to get a cluster going
#source("~/snow/snow1a.R") #this is a complicated cluster version of
##library(snow)
##setDefaultClusterOptions(master="localhost")
##cl <- makeSOCKcluster(c("localhost","localhost"))

#clusterExport(cl, "ncovs.cont") 
#clusterExport(cl, "ncovs")      
#clusterExport(cl, "form.quadz2")  
#clusterExport(cl, "fix")
#clusterExport(cl, "balance.gm3.sum")
#clusterExport(cl, "myfit")
#clusterExport(cl, "MatchBalance1")
#clusterExport(cl, "in.xdata")
#clusterExport(cl, "in.Tr")

#solution from hill.myfit1.new3c.ps5000.Rout
sv <- c(3.865330e+02, 5.204360e+01, 5.039569e+02, 5.220178e+02,
        5.061516e+01, 2.608618e+02, 4.279134e+02, 3.030325e+02, 8.374402e+00,
        5.624370e+02, 4.255614e+01, 2.215599e+02, 7.967758e+02, 1.158201e+01,
        2.266262e+02, 1.450229e+02, 6.611446e+02, 5.746659e+02, 2.703830e+01,
        7.718229e+02, 1.644440e+02, 8.985572e+02, 7.563544e+02, 5.451207e+02,
        2.579425e+02, 1.720334e+02, 5.898611e+02, 1.462342e+02, 7.641596e+02,
        9.068453e+02, 2.546357e+02, 5.413681e+02, 6.473745e+02, 3.966466e+01,
        8.220546e+02, 1.199092e+02, 3.303571e+02, 1.624907e+02, 9.325370e+02,
        7.082779e+02, 9.401388e+02)

gm2 <- GenMatch(usek2$dose400,X=modqx2$x[,2:42],M=1,pop.size=1, 
                fit.func=myfit,
                hard.generation.limit=TRUE,
                max.generations=1,                
                wait.generations=0,
#               cluster=cl,                
                starting.values=sv)

cat("With BiasAdjustment\n")
mm2=Match(usek2[,1], Tr=usek2$dose400, X=modqx2$x[,2:42], Z=modqx2$x[,3:42], Weight.matrix=gm2, replace=TRUE, M=1,BiasAdjust=TRUE)
summary(mm2)
# 13.2 (2.96)

cat("Without BiasAdjustment\n")
mm2=Match(usek2[,1], Tr=usek2$dose400, X=modqx2$x[,2:42], Z=modqx2$x[,3:42], Weight.matrix=gm2, replace=TRUE, M=1,BiasAdjust=FALSE)
summary(mm2)

cat("Print Jennifer's Balance Stats\n")
bal.gm2 <- balance.gm3.sum(mat.out=mm2,matched=TRUE,dat=cbind.data.frame(dose400=usek2$dose400,data.frame(mod.quad$x)),
                          num.cont=ncovs.cont,num.cov=ncovs,form.q=form.quadz2)
print(bal.gm2)


