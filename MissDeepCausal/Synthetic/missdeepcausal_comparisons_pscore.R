#!/usr/bin/env Rscript
################################################################
# Compare different ATE estimators that handle missing values
################################################################
# >= 13/08: "dlvm" = {confounders are X}, "dlvm2"/"dlvm3" = {confounders are Z}, 
#           "linear1"/"linear2" =  = {confounders are X}, "mf" = {confounders are Z}

# script.dir <- dirname(sys.frame(1)$ofile)
script.dir <- "~/CausalInference/Simulations/causal-inference-missing/"
setwd(script.dir)
source("./load.R")

library("optparse")
 
option_list = list(
  make_option(c("-l", "--link"), type="character", default="nonlinear3", 
              help="link for outcome and propensity models", metavar="character"),
  make_option(c("-s", "--setting"), type="character", default="dlvm", 
              help="setting", metavar="character"),
  make_option(c("-n", "--nsample"), type="integer", default=1000, 
              help="number of samples [default= %default]", metavar="number"),
  make_option(c("-p", "--ncov"), type="integer", default=10, 
              help="number of covariates [default= %default]", metavar="number"),
  make_option(c("-d", "--nlat"), type="integer", default=3, 
              help="number of latent variables [default= %default]", metavar="number"),
  make_option(c("--mech"), type="character", default="MCAR", 
              help="missing values mechanism [default= %default]", metavar="character"),
  make_option(c("--percmiss"), type="double", default=0.1, 
              help="proportion of missing values [default= %default]", metavar="number"),
  make_option(c("-m", "--mask"), type="logical", default=TRUE, 
              help="add the mask in miwae [default= %default]", metavar="boolean"),
  make_option(c("--hmiwae"), type="integer", default=128, 
              help="number of hidden units used in miwae [default= %default]", metavar="number"),
  make_option(c("--nzmult"), type="integer", default=200, 
              help="number of draws from z|x* [default= %default]", metavar="number"),
  make_option(c("--cit"), type="logical", default=FALSE, 
              help="cit [default= %default]", metavar="boolean"),
  make_option(c("--cit2"), type="logical", default=FALSE, 
              help="cit2 [default= %default]", metavar="boolean"),
  make_option(c("--ci2_imp"), type="character", default="mice", #alternative "imputePCA"
              help="ci2_imp [default= %default]", metavar="character"),
  make_option(c("--out_date"), type="character", default=as.character(Sys.Date()), 
              help="date put in output file [default= %default]", metavar="character"),
  make_option(c("-i", "--id"), type="integer", default=0, 
              help="data set id or seed [default= %default]", metavar="number")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# #######
fit_outcome_glm <- function(X, y, treat){
  df.treated <- data.frame(y = y[which(treat==1)], X = X[which(treat==1),])
  colnames(df.treated) <- c("y", paste("X", 1:dim(X)[2], sep=""))
  df.control <- data.frame(y = y[which(treat==0)], X = X[which(treat==0),])
  colnames(df.control) <- c("y", paste("X", 1:dim(X)[2], sep=""))
  
  df.treated <- df.treated[!duplicated(df.treated[,2:ncol(df.treated)]),]
  one.level.treated <- sapply(df.treated, FUN = function(x) length(unique(x))==1)
  df.treated <- df.treated[,!one.level.treated]
  fmla <- as.formula(paste0("y ~ ", paste( colnames(df.treated[,2:ncol(df.treated)]), collapse= "+")))
  lm.treated <- lm(fmla, data = df.treated)
  
  df.control <- df.control[!duplicated(df.control[,2:ncol(df.control)]),]
  one.level.control <- sapply(df.control, FUN = function(x) length(unique(x))==1)
  df.control <- df.control[,!one.level.control]
  fmla <- as.formula(paste0("y ~ ", paste( colnames(df.control[,2:ncol(df.control)]), collapse= "+")))
  lm.control <- lm(fmla, data = df.control)
  
  colnames(X) <- paste("X", 1:dim(X)[2], sep="")
  y_1.hat <- as.numeric(predict(lm.treated, X[,!one.level.treated[-1]]))
  y_0.hat <- as.numeric(predict(lm.control, X[,!one.level.control[-1]]))
  
  return(data.frame(y_1.hat = y_1.hat, y_0.hat = y_0.hat))
}



data.dir <- "./MissDeepCausal/Synthetic/Data/"
rdata.dir <- "./MissDeepCausal/Synthetic/RData/"

h <- 5 # number of hidden units in DLVM data generating model
sd <- 0.1 # noise sd used in data generation
mi.m <- 10 # number of multiple imputations

# Fetch all arguments
if (opt$mask){
  add_mask <- "_mask"
} else {
  add_mask <- ""
}

link <- opt$link
setting <- opt$setting
n <- opt$nsample
p <- opt$ncov
d <- opt$nlat
mech <- opt$mech
perc.miss <- opt$percmiss
h.miwae <- opt$hmiwae
num.samples.sir <- opt$nzmult
cit <- opt$cit
cit2 <- opt$cit2
ci2_imp <- opt$ci2_imp
id <- opt$id

set.seed(0)
V <- make_V(d, p)


e.logistic.miwae <- c()
e.logistic.full <- c()
e.logistic.cc <- c()
e.logistic.mice.woy <- c()
e.logistic.mice.wy <- c()
e.logistic.mf <- c()
estar.logistic.miwae.zmul <- c()
estar.logistic.miwae.zmul <- c()


e.grf.miwae <- c()
e.grf.full <- c()
e.grf.cc <- c()
e.grf.mice.woy <- c()
e.grf.mice.wy <- c()
e.grf.mf <- c()
estar.grf.mia <- c()
estar.grf.miwae.zmul <- c()
estar.grf.miwae.zmul <- c()


print(paste0("Propensity score estimation for data set ", id))
Z <- NULL
if (setting %in% c("mf", "dlvm2", "dlvm3")){
  Z <- read.csv(file = paste0(data.dir, setting, 
                              "_", n, "_", p, "_", d, "_seed", id, 
                              "_propNA", format(round(perc.miss, 2), nsmall=2),
                              "_z.csv"))
  
  if (dim(Z)[2] != d){
    Z <- Z[,-1]
  }
}

# Generate treatment and outcome using the codes Z
if (setting %in% c("linear1", "linear2", "mf")){
  setting_lin <- setting
  if (setting == "mf") {setting_lin <- "linear3"}
  if (link %in% c("linear3", "linear1", "linear2")){
    tmp <-gen_linear(n = n, p = p, r = d, setting = setting_lin,
                      seed = id, ps.dependence = "strong", sd = sd,
                      mechanism = "MCAR", prop.missing = perc.miss,
                      cit = cit, cio = FALSE,
                      cit2 = cit2, cio2 = FALSE,
                      ci2_imp = ci2_imp,
                      link = "log-lin",
                      V = V)
  } else {
    if (link == "nonlinear3") {
      tmp <-gen_linear(n = n, p = p, r = d, setting = setting_lin,
                       seed = id, ps.dependence = "strong", sd = sd,
                       mechanism = mech, prop.missing = perc.miss,
                       cit = cit, cio = FALSE,
                       cit2 = cit2, cio2 = FALSE,
                       ci2_imp = ci2_imp,
                       link = "nonlinear3",
                       V = V)
    }
  }
}
if (setting == "dlvm") {
    tmp <-gen_dlvm(n=n, p=p, d=d, h = h,
                   sd=sd, 
                   seed = id,
                   mechanism=mech, prop.missing = perc.miss, 
                   cit = cit, cio = FALSE,
                   cit2 = cit2, cio2 = FALSE,
                   ci2_imp = ci2_imp,
                   link = link,
                   sigma.structure = "diagonal")
} 
if (setting == "dlvm2") {
  tmp <-gen_dlvm2(n=n, p=p, d=d, h = h,
                  sd=sd, 
                  seed = id,
                  mechanism=mech, prop.missing = perc.miss, 
                  cit = cit, cio = FALSE,
                  link = link,
                  sigma.structure = "diagonal")
} 
if (setting == "dlvm3") {
  tmp <-gen_dlvm3(n=n, p=p, d=d, h = h,
                  sd=sd, 
                  seed = id,
                  mechanism=mech, prop.missing = perc.miss, 
                  cit = cit, cio = FALSE,
                  link = link,
                  sigma.structure = "diagonal")
}
y <- tmp$y
treat <- tmp$treat
e.true <- tmp$ps

data_miss <- read.csv(file = paste0(data.dir, setting, "_", 
                                    n, "_", p, "_", d,
                                    "_seed", id, 
                                    "_propNA", format(round(perc.miss, 2), nsmall=2),
                                    "_xmiss.csv"))

if (dim(data_miss)[2] > p){
  data_miss <- data_miss[,-1]
}
data_miss <- data.frame(apply(data_miss, c(1,2), FUN = function(x) ifelse(is.nan(x), NA, x)))

data_comp <- read.csv(file = paste0(data.dir, setting, "_", 
                                    n, "_", p, "_", d,
                                    "_seed", id, 
                                    "_propNA", format(round(perc.miss, 2), nsmall=2),
                                    "_xcomp.csv"))
if (dim(data_comp)[2] > p){
  data_comp <- data_comp[,-1]
}

data_zhat_miwae <- read.csv(file = paste0(data.dir, setting, "_", n, "_", p, "_", d, 
                                          "_hmiwae", h.miwae, add_mask, 
                                          "_seed", id, 
                                          "_propNA", format(round(perc.miss, 2), nsmall=2),
                                          "_zhat.csv"), 
                            sep = ";", header=F)
if (dim(data_zhat_miwae)[2] > d){
  data_zhat_miwae <- data_zhat_miwae[,-1]
}

## LOG-LIN on single imputed or full data
tt <- predict_glm(data_zhat_miwae, treat=as.factor(treat), seed=id)
e.logistic.miwae <- rbind(e.logistic.miwae, cbind(id, perc.miss, t(tt$pscore)))


if (setting %in% c("mf", "dlvm2", "dlvm3")){
  tt <- predict_glm(data.frame(Z), treat=as.factor(treat), seed=id)
  e.logistic.full <- rbind(e.logistic.full, cbind(id, perc.miss, t(tt$pscore)))
} else if (setting %in% c("linear1", "linear2", "dlvm")){
  tt <- predict_glm(data.frame(data_comp), treat=as.factor(treat), seed=id)
  e.logistic.full <- rbind(e.logistic.full, cbind(id, perc.miss, t(tt$pscore)))
}

rows.cc <- apply(data_miss, 1, FUN = function(x) sum(is.na(x))==0)
try(e.logistic.cc <- rbind(e.logistic.cc, 
                          cbind(id, perc.miss,
                                predict_glm(data_miss[rows.cc, ], treat=as.factor(treat[rows.cc ]), seed=id)$pscore)))

## FOREST on single imputed or full data
Z.m <- model.matrix(~. , data=data_zhat_miwae)
forest.W <- grf::regression_forest(Z.m, treat, tune.parameters = TRUE)
e.grf.miwae <- rbind(e.grf.miwae, cbind(id, perc.miss, t(predict(forest.W)$predictions)))

if (opt$setting %in% c("mf", "dlvm2", "dlvm3")){
  Z.m <- model.matrix(~. , data=Z)
  forest.W <- grf::regression_forest(Z.m, treat, tune.parameters = TRUE)
  e.grf.full <- rbind(e.grf.full, cbind(id, perc.miss, t(predict(forest.W)$predictions)))
} else if (opt$setting %in% c("linear1", "linear2", "dlvm")){
  Z.m <- model.matrix(~. , data=data_comp)
  forest.W <- grf::regression_forest(Z.m, treat, tune.parameters = TRUE)
  e.grf.full <- rbind(e.grf.full, cbind(id, perc.miss, t(predict(forest.W)$predictions)))
}

try(Z.m <- model.matrix(~. , data=data_miss[rows.cc, ]))  
try(forest.W  <- grf::regression_forest(Z.m, treat[rows.cc ], tune.parameters = TRUE))
try(e.grf.cc <- rbind(e.grf.cc, 
                      cbind(id, perc.miss, t(predict(forest.W)$predictions))))

## Impute the data X with mice
mice.imp.woy <- prepare_data_ate(data_miss, treat, y, 
                                 imputation.method = "mice", mi.m=mi.m, use.outcome = F)
mice.imp.wy <- prepare_data_ate(data_miss, treat, y, 
                                imputation.method = "mice", mi.m=mi.m, use.outcome = T)

res.log <- c()
res.grf <- c()
for (k in 1:mi.m){
  try(res.log <- rbind(res.log, predict_glm(mice.imp.woy$df.imp[[k]], treat=as.factor(treat), seed=id)$pscore))
  try(Z.m <- model.matrix(~. , data=mice.imp.woy$df.imp[[k]]))
  try(forest.W <- grf::regression_forest(Z.m, treat, tune.parameters = TRUE))
  try(res.grf <- rbind(res.grf, predict(forest.W)$predictions))
}
e.logistic.mice.woy <- rbind(e.logistic.mice.woy, 
                           cbind(id, perc.miss, 
                                 t(apply(data.frame(res.log), FUN = mean, MARGIN = 2))))
e.grf.mice.woy <- rbind(e.grf.mice.woy, 
                        cbind(id, perc.miss, 
                              t(apply(data.frame(res.grf), FUN = mean, MARGIN = 2))))

res.log <- c()
res.grf <- c()
for (k in 1:mi.m){
  try(res.log <- rbind(res.log, predict_glm(mice.imp.wy$df.imp[[k]], treat=as.factor(treat), seed=id)$pscore))
  try(Z.m <- model.matrix(~. , data=mice.imp.wy$df.imp[[k]]))
  try(forest.W <- grf::regression_forest(Z.m, treat, tune.parameters = TRUE))
  try(res.grf <- rbind(res.grf, predict(forest.W)$predictions))
}
e.logistic.mice.wy <- rbind(e.logistic.mice.wy, 
                          cbind(id, perc.miss, 
                                t(apply(data.frame(res.log), FUN = mean, MARGIN = 2))))
e.grf.mice.wy <- rbind(e.grf.mice.wy, 
                       cbind(id, perc.miss, 
                             t(apply(data.frame(res.grf), FUN = mean, MARGIN = 2))))


## Estimate Z_hat using low-rank matrix factorization (Kallus, Mao & Udell, 2018)
udell <- prepare_data_ate(data_miss, treat, y, 
                          imputation.method = "udell", use.outcome = F)

tt <- predict_glm(udell$df.imp, treat=as.factor(treat), seed = id)
e.logistic.mf <- rbind(e.logistic.mf, cbind(id, perc.miss, t(tt$pscore)))

Z.m <- model.matrix(~. , data=udell$df.imp)
forest.W <- grf::regression_forest(Z.m, treat, tune.parameters = TRUE)
e.grf.mf <- rbind(e.grf.mf, cbind(id, perc.miss, t(predict(forest.W)$predictions)))


## MIA
mia <- prepare_data_ate(data_miss, treat, y, 
                        imputation.method = "mia.grf.ate", use.outcome = F)
Z.m <- model.matrix(~. , data=mia$df.imp)
forest.W <- grf::regression_forest(Z.m, treat,tune.parameters = TRUE)
estar.grf.mia <- rbind(estar.grf.mia, cbind(id, perc.miss, t(predict(forest.W)$predictions)))


# ####### Add estimation on sampled Z.hat
w.hat.grf.mul <- c()
w.hat.log.mul <- c()
for (k in 0:(num.samples.sir-1)){
  cat(paste(k," "))
  data_zhat_miwae <- read.csv(file = paste0(data.dir, setting, "_", n, "_", p, 
                                            "_", d, "_hmiwae", h.miwae, add_mask, "_seed", id, 
                                            "_propNA", format(round(perc.miss, 2), nsmall=2),
                                            "_zhat_m", k, ".csv"),
                              sep = ";", header=F)
  if (dim(data_zhat_miwae)[2] > d){
    data_zhat_miwae <- data_zhat_miwae[,-1]
  }
  ## grf
  Z.m <- model.matrix(~. , data=data_zhat_miwae)
  forest.W <- regression_forest(Z.m, treat, tune.parameters = TRUE)
  w.hat <- predict(forest.W)$predictions
  w.hat.grf.mul <- rbind(w.hat.grf.mul, w.hat)
  
  ## log-lin
  w.hat <- predict_glm(data_zhat_miwae, as.factor(treat), seed=id)$pscore
  w.hat.log.mul <- rbind(w.hat.log.mul, w.hat)
  
}
## log-lin
estar.logistic.miwae.zmul <- rbind(estar.logistic.miwae.zmul, 
                                   cbind(id, perc.miss, t(apply(w.hat.log.mul, 2, mean))))

## grf
estar.grf.miwae.zmul <- rbind(estar.grf.miwae.zmul, 
                              cbind(id, perc.miss, t(apply(w.hat.grf.mul, 2, mean))))


print(paste0("finished estimation for data set ", id))

######## Merge results
change_colnames <- function(df, col.names){
  colnames(df) <- col.names
  return(df)
}

col.names <- c("id", "prop.missing", paste0("e",1:n))

e.true <- change_colnames(data.frame(cbind(id, perc.miss, t(e.true))), col.names)

e.logistic.mice.woy <- change_colnames(data.frame(e.logistic.mice.woy), col.names)
e.logistic.mice.wy <- change_colnames(data.frame(e.logistic.mice.wy), col.names)
e.logistic.mf <- change_colnames(data.frame(e.logistic.mf), col.names)
e.logistic.miwae <- change_colnames(data.frame(e.logistic.miwae), col.names)
e.logistic.full <- change_colnames(data.frame(e.logistic.full), col.names)
estar.logistic.miwae.zmul <- change_colnames(data.frame(estar.logistic.miwae.zmul), col.names)
try(e.logistic.cc <- change_colnames(data.frame(e.logistic.cc), col.names))

e.grf.mice.woy <- change_colnames(data.frame(e.grf.mice.woy), col.names)
e.grf.mice.wy <- change_colnames(data.frame(e.grf.mice.wy), col.names)
e.grf.mf <- change_colnames(data.frame(e.grf.mf), col.names)
e.grf.miwae <- change_colnames(data.frame(e.grf.miwae), col.names)
e.grf.full <- change_colnames(data.frame(e.grf.full), col.names)
estar.grf.mia <- change_colnames(data.frame(estar.grf.mia), col.names)
estar.grf.miwae.zmul <- change_colnames(data.frame(estar.grf.miwae.zmul), col.names)

try(e.grf.cc <- change_colnames(data.frame(e.grf.cc), col.names))


res <- rbind(cbind(method = rep("true", dim(e.true)[1]), e.true),
             cbind(method = rep("logistic.mice.woy",dim(e.logistic.mice.woy)[1]), e.logistic.mice.woy),
             cbind(method = rep("logistic.mice.wy",dim(e.logistic.mice.wy)[1]), e.logistic.mice.wy),
             cbind(method = rep("logistic.mf",dim(e.logistic.mf)[1]), e.logistic.mf),
             cbind(method = rep("logistic.miwae",dim(e.logistic.miwae)[1]), e.logistic.miwae),
             cbind(method = rep("logistic.full",dim(e.logistic.full)[1]), e.logistic.full),
             cbind(method = rep("logistic.miwae.zmul",dim(estar.logistic.miwae.zmul)[1]), estar.logistic.miwae.zmul),
             cbind(method = rep("grf.mice.woy",dim(e.grf.mice.woy)[1]), e.grf.mice.woy),
             cbind(method = rep("grf.mice.wy",dim(e.grf.mice.wy)[1]), e.grf.mice.wy),
             cbind(method = rep("grf.mf",dim(e.grf.mf)[1]), e.grf.mf),
             cbind(method = rep("grf.miwae",dim(e.grf.miwae)[1]), e.grf.miwae),
             cbind(method = rep("grf.full",dim(e.grf.full)[1]), e.grf.full),
             cbind(method = rep("grf.mia",dim(estar.grf.mia)[1]), estar.grf.mia),
             cbind(method = rep("grf.miwae.zmul",dim(estar.grf.miwae.zmul)[1]), estar.grf.miwae.zmul))

print(paste0("first rbind passed for data set ", id))

try(res <- rbind(res, 
                 cbind(method = rep("grf.cc",dim(e.grf.cc)[1]), e.grf.cc),
                 cbind(method = rep("logistic.cc",dim(e.logistic.cc)[1]), e.logistic.cc)))

print(paste0("second rbind passed for data set ", id))



save(res, file = paste0(rdata.dir, opt$out_date, 
                        "_",link,"_",setting, 
                        "_n", n, "_p", p, "_d", d, "_hmiwae", h.miwae,
                        add_mask,
                        "_sir", num.samples.sir,
                        "_propNA", format(round(perc.miss, 2), nsmall=2),
                        "_cit2", as.logical(cit2),
                        "_results_pscore_seed",id,
                        ".RData"))

print(paste0("saved results for data set ", id))

