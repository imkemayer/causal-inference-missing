#!/usr/bin/env Rscript
################################################################
# Compare different ATE estimators that handle missing values
################################################################
# >= 13/08: "dlvm" = {confounders are X}, "dlvm2"/"dlvm3" = {confounders are Z}, 
#           "linear1"/"linear2" =  = {confounders are X}, "mf" = {confounders are Z}

# script.dir <- dirname(sys.frame(1)$ofile)
script.dir <- "~/CausalInference/Simulations/causal-inference-missing/"
setwd(script.dir)
source("./MissDeepCausal/Synthetic/Scripts/load.R")

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
  make_option(c("-w", "--w"), type="logical", default=FALSE, 
              help="used W in miwae [default= %default]", metavar="boolean"),
  make_option(c("--wy"), type="logical", default=FALSE, 
              help="w and y used in miwae [default= %default]", metavar="boolean"),
  make_option(c("--dmiwae"), type="integer", default=3, 
              help="number of latent variables estimated in miwae [default= %default]", metavar="number"),
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
  make_option(c("--offset"), type="logical", default=TRUE, 
              help="offset in regression [default= %default]", metavar="boolean"),
  make_option(c("--out_date"), type="character", default=as.character(Sys.Date()), 
              help="date put in output file [default= %default]", metavar="character"),
  make_option(c("-i", "--id"), type="integer", default=0, 
              help="data set id or seed [default= %default]", metavar="number")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# #######
regression_lin <- function(X, y, treat, offset=TRUE){
  df <- data.frame(cbind(X,"y"=y, "w" = treat))
  if (offset){
    model <- lm(y ~ . + w, data = df)
  } else {
    model <- lm(y ~ . + w + 0, data = df)
  }
  return(coefficients(model)["w"])
}


prepare_data_ate <- function(df, w, y, imputation.method, mi.m, 
                             mask,
                             use.outcome, use.interaction,
                             r_udell = NULL){ 
  df.imp <- fitted <- NULL                                      
  if (tolower(imputation.method) %in% c("mean", "mean.grf", "mean.grf.ate")){
    df.imp <- get_imputeMean(data.frame(df))
  } 
  if (tolower(imputation.method) =="pca"){
    if (use.outcome){
      df.imp <- data.frame(df)
      if (sum(is.na(df))>0){
        df.y <- data.frame(cbind(df, y = y))
        ncomp <- estim_ncpPCA(df.y, ncp.min=1, ncp.max=dim(df)[2]+1)
        df.imp <- get_imputePC(df.y, ncp = ncomp$ncp, seed=i)
        df.imp <- df.imp[,1:dim(df)[2]]
      }
      
    } else {
      df.imp <- data.frame(df)
      if (sum(is.na(df))>0){
        ncomp <- try(estim_ncpPCA(data.frame(df), ncp.min=1, ncp.max=dim(df)[2]))
        df.imp <- get_imputePC(data.frame(df), ncp = ncomp$ncp, seed=i)
      }
    }
  } 
  if (tolower(imputation.method) =="mice"){
    if(use.outcome){
      df.y <- data.frame(cbind(df, y = y))
      df.imp <- get_MICE(df.y, m = mi.m, idx = 1:dim(df)[2])
    } else {
      df.imp <- get_MICE(data.frame(df), m = mi.m)
    }
  } 
  if (tolower(imputation.method) == "missforest"){
    if (use.outcome){
      df.y <- data.frame(cbind(df, y = y))
      df.imp <- get_MISSFOREST(df.y)
      df.imp <- df.imp[,1:dim(df)[2]]
      
    } else {
      df.imp <- get_MISSFOREST(data.frame(df))
    }
  } 
  if (tolower(imputation.method) %in% c("mia.grf", "mia.grf.ate")){
    df.imp <- get_imputeInf(data.frame(df))
  } 
  if (tolower(imputation.method) == "saem"){
    df.imp <- data.frame(df)
    df.num <- data.frame(sapply(data.frame(df), as.numeric))
    if (max(as.integer(w))!= 1){
      w <- as.integer(w)
      idx.ones <- which(w==max(w))
      idx.zeros <- which(w!=max(w))
      w[idx.ones] <- 1
      w[idx.zeros] <- 0
    }
    Z.misaem <- predict_misaem(df.num, 
                               w, 
                               pattern = mask,
                               use.interaction = use.interaction)
    fitted <- Z.misaem$fitted
  } 
  if (tolower(imputation.method) == "udell") {
    #pca_soft <- recover_pca_gaussian_cv(as.matrix(df), 1:dim(df)[2], nfolds = 3)
    #df.imp <- data.frame(pca_soft$Uhat)
    if (all(sapply(data.frame(df), is.numeric))){
      if (is.null(r_udell)){
        pca_soft <- recover_pca_gaussian_cv(as.matrix(df), 1:dim(df)[2], nfolds = 3)
      } else {
        pca_soft <- recover_pca_gaussian(as.matrix(df), r_udell)
      }
      df.imp <- data.frame(pca_soft$Uhat)
    } else {
      #var.type <- sapply(data.frame(df), FUN = function(x) {if_else(is.numeric(x), "gaussian", "binomial")})
      #df[,var.type=="binomial"] <- sapply(df[,var.type=="binomial"], as.integer)
      #glrm_mod <- mimi::mimi(y=data.frame(df), model = "low-rank",
      #                          var.type = var.type, max.rank = dim(df), algo="bcgd", lambda1=1)
      #ncomp <- missMDA::estim_ncpFAMD(df, ncp.max = dim(df)[2]-1)
    }
  }
  return(list(df.imp = df.imp, fitted = fitted))
}


data.dir <- "./MissDeepCausal/Synthetic/Data/"
rdata.dir <- "./MissDeepCausal/Synthetic/RData/"

h <- 5 # number of hidden units in DLVM data generating model
sd <- 0.1 # noise sd used in data generation
mi.m <- 10 # number of multiple imputations

link <- opt$link
setting <- opt$setting
n <- opt$nsample
p <- opt$ncov
d <- opt$nlat
mech <- opt$mech
perc.miss <- opt$percmiss
d.miwae <- opt$dmiwae
h.miwae <- opt$hmiwae
num.samples.sir <- opt$nzmult
cit <- opt$cit
cit2 <- opt$cit2
ci2_imp <- opt$ci2_imp
offset <- opt$offset

for (id in 0:opt$id){
# opt <- data.frame("link" = "linear3", "setting" = "mf", "nsample" = 10000, "ncov" = 10, "nlat" = 3,
#                   "mech" = "MCAR", "percmiss" = 0, "cit"=FALSE, "cit2"=FALSE, "ci2_imp" = "mice",
#                   "mask"= FALSE, "w"=FALSE, "dmiwae"=10, "hmiwae"=128, "nzmult"=200,
#                   "offset"=TRUE,
#                   "out_date"=as.character(Sys.Date()), "id"=id)
# opt[,c("link", "setting","mech","ci2_imp", "out_date")] <- sapply(opt[,c("link", "setting","mech","ci2_imp", "out_date")], as.character)

# Fetch all arguments
if (opt$mask){
  add_mask <- "_mask"
} else {
  add_mask <- ""
}

if (opt$w){
  use.w <- "_W"
} else {
  use.w <- ""
}

if (opt$wy){
  use.wy <- "_WY"
} else {
  use.wy <- ""
}




set.seed(0)
V <- make_V(d, p)


lin.miwae <- c()
lin.full <- c()
lin.cc <- c()
lin.mice.woy <- c()
lin.mice.wy <- c()
lin.mi.miwae <- c()
lin.mf <- c()
lin.mf.cv <- c()
lin.mean <- c()

lin.miwae.agg <- c()




print(paste0("ATE estimation for data set ", id, " (link: ", link, ", setting: ", setting,", p: ", p, ", propNA: ",perc.miss, ")"))
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

if (opt$wy){
  data_miss <- read.csv(file = paste0(data.dir, setting, "_", link, "_",
                                    n, "_", p, "_", d,
                                    "_seed", id, 
                                    "_propNA", format(round(perc.miss, 2), nsmall=2),
                                    "_xmiss.csv"))
} else {
  data_miss <- read.csv(file = paste0(data.dir, setting, "_", 
                                    n, "_", p, "_", d,
                                    "_seed", id, 
                                    "_propNA", format(round(perc.miss, 2), nsmall=2),
                                    "_xmiss.csv"))
}
if (dim(data_miss)[2] > p){
  data_miss <- data_miss[,-1]
}
data_miss <- data.frame(apply(data_miss, c(1,2), FUN = function(x) ifelse(is.nan(x), NA, x)))


if (opt$wy){
  data_comp <- read.csv(file = paste0(data.dir, setting,  "_", link, "_", 
                                    n, "_", p, "_", d,
                                    "_seed", id, 
                                    "_propNA", format(round(perc.miss, 2), nsmall=2),
                                    "_xcomp.csv"))
} else {
  data_comp <- read.csv(file = paste0(data.dir, setting, "_", 
                                    n, "_", p, "_", d,
                                    "_seed", id, 
                                    "_propNA", format(round(perc.miss, 2), nsmall=2),
                                    "_xcomp.csv"))
}
if (dim(data_comp)[2] > p){
  data_comp <- data_comp[,-1]
}
data_zhat_miwae <- NULL
if (opt$w){
  data_zhat_miwae <- read.csv(file = paste0(data.dir, setting,  "_", link,"_", n, "_", p, "_", d, 
                                              "_W_dmiwae", d.miwae, "_hmiwae", h.miwae, add_mask, 
                                              "_seed", id, 
                                              "_propNA", format(round(perc.miss, 2), nsmall=2),
                                              "_zhat.csv"), 
                                sep = ";", header=F)
} else if (opt$wy){
  data_zhat_miwae <- read.csv(file = paste0(data.dir, setting,  "_", link,"_", n, "_", p, "_", d, 
                                            "_WY_dmiwae", d.miwae, "_hmiwae", h.miwae, add_mask, 
                                            "_seed", id, 
                                            "_propNA", format(round(perc.miss, 2), nsmall=2),
                                            "_zhat.csv"), 
                              sep = ";", header=F)
} else {
  
  try(data_zhat_miwae <- read.csv(file = paste0(data.dir, setting, "_", n, "_", p, "_", d, 
                                            "_dmiwae", d.miwae, "_hmiwae", h.miwae, add_mask, 
                                            "_seed", id, 
                                            "_propNA", format(round(perc.miss, 2), nsmall=2),
                                            "_zhat.csv"), 
                              sep = ";", header=F))
  if (is.null(data_zhat_miwae)){
    try(data_zhat_miwae <- read.csv(file = paste0(data.dir, setting,  "_", link,"_", n, "_", p, "_", d, 
                                            "_dmiwae", d.miwae, "_hmiwae", h.miwae, add_mask, 
                                            "_seed", id, 
                                            "_propNA", format(round(perc.miss, 2), nsmall=2),
                                            "_zhat.csv"), 
                              sep = ";", header=F))
  }
}
if (dim(data_zhat_miwae)[2] > d.miwae){
  data_zhat_miwae <- data_zhat_miwae[,-1]
}

## Propensity scores
mean.imp <- prepare_data_ate(data_miss, treat, y, imputation.method = "mean")
e.logistic.mean <- predict_glm(mean.imp$df.imp, treat=as.factor(treat), seed=id)$pscore

e.logistic.miwae <- predict_glm(data_zhat_miwae, treat=as.factor(treat), seed=id)$pscore
if (setting %in% c("mf", "dlvm2", "dlvm3")){
  e.logistic.full <- predict_glm(data.frame(Z), treat=as.factor(treat), seed=id)$pscore
} else if (setting %in% c("linear1", "linear2", "dlvm")){
  e.logistic.full <- predict_glm(data.frame(data_comp), treat=as.factor(treat), seed=id)$pscore
}
rows.cc <- apply(data_miss, 1, FUN = function(x) sum(is.na(x))==0)
try(e.logistic.cc <- predict_glm(data_miss[rows.cc, ], treat=as.factor(treat[rows.cc ]), seed=id)$pscore)

## Impute the data X with mice
mice.imp.woy <- prepare_data_ate(data_miss, treat, y, 
                                 imputation.method = "mice", mi.m=mi.m, use.outcome = F)
mice.imp.wy <- prepare_data_ate(data_miss, treat, y, 
                                imputation.method = "mice", mi.m=mi.m, use.outcome = T)

res.log <- c()
for (k in 1:mi.m){
  try(res.log <- rbind(res.log, predict_glm(mice.imp.woy$df.imp[[k]], treat=as.factor(treat), seed=id)$pscore))
}
e.logistic.mice.woy <- apply(data.frame(res.log), FUN = mean, MARGIN = 2)

res.log <- c()
for (k in 1:mi.m){
  try(res.log <- rbind(res.log, predict_glm(mice.imp.wy$df.imp[[k]], treat=as.factor(treat), seed=id)$pscore))
}
e.logistic.mice.wy <- apply(data.frame(res.log), FUN = mean, MARGIN = 2)

## Estimate Z_hat using low-rank matrix factorization (Kallus, Mao & Udell, 2018)
udell.cv <- prepare_data_ate(data_miss, treat, y, 
                            imputation.method = "udell", use.outcome = F)
udell <- prepare_data_ate(data_miss, treat, y, 
                            imputation.method = "udell", use.outcome = F, r_udell = d)

e.logistic.mf.cv <- predict_glm(udell.cv$df.imp, treat=as.factor(treat), seed = id)$pscore
e.logistic.mf <- predict_glm(udell$df.imp, treat=as.factor(treat), seed = id)$pscore


## Linear regressions
tt <- regression_lin(data_zhat_miwae, y=y, treat=treat, offset = offset)
lin.miwae <- rbind(lin.miwae, cbind(tt, perc.miss, id, "glm"))
tt <- regression_lin(cbind(data_zhat_miwae, e.logistic.miwae), y=y, treat=treat, offset = offset)
lin.miwae <- rbind(lin.miwae, cbind(tt, perc.miss, id, "glm.wPS"))


tt <- regression_lin(mean.imp$df.imp, y=y, treat=treat, offset = offset)
lin.mean <- rbind(lin.mean, cbind(tt, perc.miss, id, "glm"))
tt <- regression_lin(cbind(mean.imp$df.imp, e.logistic.mean), y=y, treat=treat, offset = offset)
lin.mean <- rbind(lin.mean, cbind(tt, perc.miss, id, "glm.wPS"))


if (setting %in% c("mf", "dlvm2", "dlvm3")){
  tt <- regression_lin(data.frame(Z), y=y, treat=treat, offset = offset)
  lin.full <- rbind(lin.full, cbind(tt, perc.miss, id, "glm"))
  tt <- regression_lin(data.frame(cbind(Z, e.logistic.full)), y=y, treat=treat, offset = offset)
  lin.full <- rbind(lin.full, cbind(tt, perc.miss, id, "glm.wPS"))
  
} else if (setting %in% c("linear1", "linear2", "dlvm")){
  tt <- regression_lin(data.frame(data_comp), y=y, treat=treat, offset = offset)
  lin.full <- rbind(lin.full, cbind(tt, perc.miss, id, "glm"))
  tt <- regression_lin(data.frame(cbind(data_comp, e.logistic.full)), y=y, treat=treat, offset = offset)
  lin.full <- rbind(lin.full, cbind(tt, perc.miss, id, "glm.wPS"))
}

rows.cc <- apply(data_miss, 1, FUN = function(x) sum(is.na(x))==0)
try(lin.cc <- rbind(lin.cc, 
                    cbind(regression_lin(data_miss[rows.cc, ], y=y[rows.cc ], treat=treat[rows.cc ], offset = offset),
                          perc.miss, id, "glm")))
try(lin.cc <- rbind(lin.cc, 
                    cbind(regression_lin(cbind(data_miss[rows.cc, ], e.logistic.cc), y=y[rows.cc ], treat=treat[rows.cc ], offset = offset),
                          perc.miss, id, "glm.wPS")))



res <- c()
res.ps <- c()
for (k in 1:mi.m){
  try(res <- rbind(res,regression_lin(X = mice.imp.woy$df.imp[[k]], y = y, treat = treat, offset = offset)))
  try(res.ps <- rbind(res.ps,regression_lin(X = cbind(mice.imp.woy$df.imp[[k]], e.logistic.mice.woy), y = y, treat = treat, offset = offset)))
}
lin.mice.woy <- rbind(lin.mice.woy, cbind(t(apply(data.frame(res), FUN = mean, MARGIN = 2)),
                                                      perc.miss, id, "glm"))
lin.mice.woy <- rbind(lin.mice.woy, cbind(t(apply(data.frame(res.ps), FUN = mean, MARGIN = 2)),
                                                      perc.miss, id, "glm.wPS"))

res <- c()
res.ps <- c()
for (k in 1:mi.m){
  try(res <- rbind(res,regression_lin(X = mice.imp.wy$df.imp[[k]], y = y, treat = treat, offset = offset)))
  try(res.ps <- rbind(res.ps,regression_lin(X = cbind(mice.imp.wy$df.imp[[k]], e.logistic.mice.wy), y = y, treat = treat, offset = offset)))
}
lin.mice.wy <- rbind(lin.mice.wy, cbind(t(apply(data.frame(res), FUN = mean, MARGIN = 2)),
                                                      perc.miss, id, "glm"))
lin.mice.wy <- rbind(lin.mice.wy, cbind(t(apply(data.frame(res.ps), FUN = mean, MARGIN = 2)),
                                                      perc.miss, id, "glm.wPS"))


## Estimate Z_hat using low-rank matrix factorization (Kallus, Mao & Udell, 2018)
tt <- regression_lin(udell.cv$df.imp, y=y, treat=treat, offset = offset)
lin.mf.cv <- rbind(lin.mf.cv, cbind(tt, perc.miss, id, "glm"))
tt <- regression_lin(cbind(udell.cv$df.imp, e.logistic.mf.cv), y=y, treat=treat, offset = offset)
lin.mf.cv <- rbind(lin.mf.cv, cbind(tt, perc.miss, id, "glm.wPS"))

tt <- regression_lin(udell$df.imp, y=y, treat=treat, offset = offset)
lin.mf <- rbind(lin.mf, cbind(tt, perc.miss, id, "glm"))
tt <- regression_lin(cbind(udell$df.imp, e.logistic.mf), y=y, treat=treat, offset = offset)
lin.mf <- rbind(lin.mf, cbind(tt, perc.miss, id, "glm.wPS"))



change_colnames <- function(df, col.names){
  colnames(df) <- col.names
  return(df)
}



col.names <- c("tau", "prop.missing", "seed", "regression")


lin.miwae <- change_colnames(data.frame(lin.miwae), col.names)
lin.full <-  change_colnames(data.frame(lin.full), col.names)
try(lin.cc <- change_colnames(data.frame(lin.cc), col.names))
lin.mice.woy <-  change_colnames(data.frame(lin.mice.woy), col.names)
lin.mice.wy <- change_colnames(data.frame(lin.mice.wy), col.names)
lin.mf <- change_colnames(data.frame(lin.mf), col.names)
lin.mf.cv <- change_colnames(data.frame(lin.mf.cv), col.names)
lin.mean <- change_colnames(data.frame(lin.mean), col.names)



print(paste0('renamed first set of results for data set ', id))

# ####### Add estimation on sampled Z.hat

lin.agg <- c()
lin.agg.ps <- c()
w.hat.log.mul <- c()

for (k in 0:(num.samples.sir-1)){
  cat(paste(k," "))
  data_zhat_miwae <- NULL
  if (opt$w){
      try(data_zhat_miwae <- read.csv(file = paste0(data.dir, setting, "_", link, "_", n, "_", p, "_", d, 
                                                "_W_dmiwae", d.miwae, "_hmiwae", h.miwae, add_mask, 
                                                "_seed", id, 
                                                "_propNA", format(round(perc.miss, 2), nsmall=2),
                                                "_zhat_m",k,".csv"), 
                                  sep = ";", header=F))
  } else if (opt$wy) {
    try(data_zhat_miwae <- read.csv(file = paste0(data.dir, setting, "_", link, "_", n, "_", p, "_", d, 
                                                "_WY_dmiwae", d.miwae, "_hmiwae", h.miwae, add_mask, 
                                                "_seed", id, 
                                                "_propNA", format(round(perc.miss, 2), nsmall=2),
                                                "_zhat_m",k,".csv"), 
                                  sep = ";", header=F))
  } else {
      
    try(data_zhat_miwae <- read.csv(file = paste0(data.dir, setting, "_", n, "_", p, "_", d, 
                                                "_dmiwae", d.miwae, "_hmiwae", h.miwae, add_mask, 
                                                "_seed", id, 
                                                "_propNA", format(round(perc.miss, 2), nsmall=2),
                                                "_zhat_m",k,".csv"), 
                                  sep = ";", header=F))
    if (is.null(data_zhat_miwae)){
        try(data_zhat_miwae <- read.csv(file = paste0(data.dir, setting, "_", link, "_", n, "_", p, "_", d, 
                                                "_dmiwae", d.miwae, "_hmiwae", h.miwae, add_mask, 
                                                "_seed", id, 
                                                "_propNA", format(round(perc.miss, 2), nsmall=2),
                                                "_zhat_m",k,".csv"), 
                                  sep = ";", header=F))
    }
  }

  if (!is.null(data_zhat_miwae) & dim(data_zhat_miwae)[1] == length(treat)){  
    if (dim(data_zhat_miwae)[2] > d.miwae){
      data_zhat_miwae <- data_zhat_miwae[,-1]
    }
    w.hat <- predict_glm(data_zhat_miwae, as.factor(treat), seed=id)$pscore
    w.hat.log.mul <- rbind(w.hat.log.mul, w.hat)

    tt <- regression_lin(data_zhat_miwae, y, treat, offset = offset)
    lin.agg <- rbind(lin.agg, cbind(tt, perc.miss, id, "glm"))
    tt <- regression_lin(cbind(data_zhat_miwae, w.hat), y, treat, offset = offset)
    lin.agg.ps <- rbind(lin.agg.ps, cbind(tt, perc.miss, id, "glm.wPS"))
  }
}

## aggregate tau estimates
lin.miwae.agg <- rbind(lin.miwae.agg, cbind(apply(data.frame(as.numeric(lin.agg[,1])),2,mean), perc.miss, id, "glm"))
lin.miwae.agg <- rbind(lin.miwae.agg, cbind(apply(data.frame(as.numeric(lin.agg.ps[,1])),2,mean), perc.miss, id, "glm"))

lin.miwae.agg <- change_colnames(lin.miwae.agg, col.names)


print(paste0("finished estimation for data set ", id))

######## Aggregate results

tau.lin <- rbind(cbind(lin.mice.woy, method = rep("lin.mice.woy",dim(lin.mice.woy)[1])),
                 cbind(lin.mice.wy, method = rep("lin.mice.wy",dim(lin.mice.wy)[1])),
                 cbind(lin.mf, method = rep("lin.mf",dim(lin.mf)[1])),
                 cbind(lin.mean, method = rep("lin.mean",dim(lin.mean)[1])),
                 cbind(lin.mf.cv, method = rep("lin.mf.cv",dim(lin.mf.cv)[1])),
                 cbind(lin.miwae, method = rep("lin.miwae",dim(lin.miwae)[1])),
                 cbind(lin.full, method = rep("lin.Z_true",dim(lin.full)[1])))
print(paste0("first lin rbind passed for data set ", id))

try(tau.lin <- rbind(tau.lin, 
                     cbind(lin.cc, method = rep("lin.cc",dim(lin.cc)[1]))))
print(paste0("second lin rbind passed for data set ", id))



try(tau.lin <- rbind(tau.lin,
                     cbind(data.frame(lin.miwae.agg), method = rep("lin.miwae.agg",dim(lin.miwae.agg)[1]))))
print(paste0("third lin rbind passed for data set ", id))


print(paste0(rdata.dir, opt$out_date, 
                            "_", link, "_", setting, 
                            "_n", n, "_p", p, "_d", d, 
                            use.w, use.wy, "_dmiwae", d.miwae, "_hmiwae", h.miwae,
                            add_mask,
                            "_sir", num.samples.sir,
                            "_propNA", format(round(perc.miss, 2), nsmall=2),
                            "_cit2", as.logical(cit2),
                            "_ci2Imp", ci2_imp,
                            "_offset", as.logical(offset),
                            "_results_regression_seed", id,
                            ".RData"))

save(tau.lin, file = paste0(rdata.dir, opt$out_date, 
                            "_", link, "_", setting, 
                            "_n", n, "_p", p, "_d", d, 
                            use.w, use.wy, "_dmiwae", d.miwae, "_hmiwae", h.miwae,
                            add_mask,
                            "_sir", num.samples.sir,
                            "_propNA", format(round(perc.miss, 2), nsmall=2),
                            "_cit2", as.logical(cit2),
                            "_ci2Imp", ci2_imp,
                            "_offset", as.logical(offset),
                            "_results_regression_seed", id,
                            ".RData"))

}


