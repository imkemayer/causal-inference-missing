#!/usr/bin/env Rscript
################################################################
# Compare different ATE estimators that handle missing values
################################################################

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
              help="w used in miwae [default= %default]", metavar="boolean"),
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
  make_option(c("--forest"), type="character", default="grf", #alternative "ranger"
              help="forest [default= %default]", metavar="character"),
  make_option(c("--out_date"), type="character", default=as.character(Sys.Date()), 
              help="date put in output file [default= %default]", metavar="character"),
  make_option(c("-i", "--id"), type="integer", default=0, 
              help="data set id or seed [default= %default]", metavar="number")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

################################
# # Define some helper functions
################################

### Function: fit conditional response surfaces models with linear regression
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

### Function: prepare data by imputing missing values or estimating latent factors
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

### Function: rename entry labels in result vectors
change_colnames <- function(df, col.names){
  colnames(df) <- col.names
  return(df)
}

dr.col.names <- c("dr", "se", "prop.missing", "seed", "regression")
ipw.col.names <- c("ipw1", "ipw2", "se", "prop.missing", "seed", "regression")


########################
# # Set global variables
########################


data.dir <- "./MissDeepCausal/Synthetic/Data/" # location of csv files to be loaded for analysis
rdata.dir <- "./MissDeepCausal/Synthetic/RData/" # location where to save analysis results (RData files)

h <- 5 # number of hidden units in DLVM data generating model
sd <- 0.1 # noise sd used in data generation
mi.m <- 10 # number of multiple imputations

### Fetch all arguments from script call
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
use.w <- opt$w
d.miwae <- opt$dmiwae
h.miwae <- opt$hmiwae
num.samples.sir <- opt$nzmult
cit <- opt$cit
cit2 <- opt$cit2
ci2_imp <- opt$ci2_imp
forest <- opt$forest
id <- opt$id

### Set temporary variables
set.seed(0)
V <- make_V(d, p) # only relevant for `mf` setting (this is the loadings matrix)


### Initialize the vectors that will contain analysis results
dr.loglin.miwae <- ipw.loglin.miwae <- c()
dr.loglin.full <- ipw.loglin.full <- c()
dr.loglin.cc <- ipw.loglin.cc <- c()
dr.loglin.mice.woy <- ipw.loglin.mice.woy <- c()
dr.loglin.mice.wy <- ipw.loglin.mice.wy <- c()
dr.loglin.mf <- ipw.loglin.mf <- c()
dr.loglin.mf.cv <- ipw.loglin.mf.cv <- c()


dr.grf.miwae <- ipw.grf.miwae <- c()
dr.grf.full <- ipw.grf.full <- c()
dr.grf.cc <- ipw.grf.cc <- c()
dr.grf.mice.woy <- ipw.grf.mice.woy <- c()
dr.grf.mice.wy <- ipw.grf.mice.wy <- c()
dr.grf.mf <- ipw.grf.mf <- c()
dr.grf.mf.cv <- ipw.grf.mf.cv <- c()
dr.grf.mia <- ipw.grf.mia <- c()


dr.grf.miwae.zmul <- ipw.grf.miwae.zmul <- c()
dr.loglin.miwae.zmul <- ipw.loglin.miwae.zmul <- c()
dr.grf.miwae.agg <- ipw.grf.miwae.agg <- c()
dr.loglin.miwae.agg <- ipw.loglin.miwae.agg <- c()



print(paste0("ATE estimation for data set ", id))
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

# if (d.miwae == 3){
#   if (use.w){
#     data_zhat_miwae <- read.csv(file = paste0(data.dir, setting, "_", n, "_", p, "_", d, 
#                                               "_W_hmiwae", h.miwae, add_mask, 
#                                               "_seed", id, 
#                                               "_propNA", format(round(perc.miss, 2), nsmall=2),
#                                               "_zhat.csv"), 
#                                 sep = ";", header=F)
#   } else {
#     data_zhat_miwae <- read.csv(file = paste0(data.dir, setting, "_", n, "_", p, "_", d, 
#                                               "_hmiwae", h.miwae, add_mask, 
#                                               "_seed", id, 
#                                               "_propNA", format(round(perc.miss, 2), nsmall=2),
#                                               "_zhat.csv"), 
#                                 sep = ";", header=F)
#   }
# } else {
  if (use.w){
    data_zhat_miwae <- read.csv(file = paste0(data.dir, setting, "_", n, "_", p, "_", d, 
                                              "_W_dmiwae", d.miwae, "_hmiwae", h.miwae, add_mask, 
                                              "_seed", id, 
                                              "_propNA", format(round(perc.miss, 2), nsmall=2),
                                              "_zhat.csv"), 
                                sep = ";", header=F)
  } else {
    data_zhat_miwae <- read.csv(file = paste0(data.dir, setting, "_", n, "_", p, "_", d, 
                                              "_dmiwae", d.miwae, "_hmiwae", h.miwae, add_mask, 
                                              "_seed", id, 
                                              "_propNA", format(round(perc.miss, 2), nsmall=2),
                                              "_zhat.csv"), 
                                sep = ";", header=F)
  }
# }
if (dim(data_zhat_miwae)[2] > d.miwae){
  data_zhat_miwae <- data_zhat_miwae[,-1]
}

## LOG-LIN
tt <- dr(data_zhat_miwae, outcome=y, treat=treat, 
         ps.method="glm",  
         out.method = "glm")
dr.loglin.miwae <- rbind(dr.loglin.miwae, cbind(tt, perc.miss, id, "glm"))

tt <- ipw(data_zhat_miwae, outcome=y, treat=treat,
          ps.method="glm")
ipw.loglin.miwae <- rbind(ipw.loglin.miwae, cbind(tt, perc.miss, id, "glm"))

if (setting %in% c("mf", "dlvm2", "dlvm3")){
  tt <- dr(data.frame(Z), outcome=y, treat=treat,
           ps.method="glm",
           out.method = "glm")
  dr.loglin.full <- rbind(dr.loglin.full, cbind(tt, perc.miss, id, "glm"))
  
  tt <- ipw(data.frame(Z), outcome=y, treat=treat,
            ps.method="glm")
  ipw.loglin.full <- rbind(ipw.loglin.full, cbind(tt, perc.miss, id, "glm"))
} else if (setting %in% c("linear1", "linear2", "dlvm")){
  tt <- dr(data.frame(data_comp), outcome=y, treat=treat,
           ps.method="glm",
           out.method = "glm")
  dr.loglin.full <- rbind(dr.loglin.full, cbind(tt, perc.miss, id, "glm"))
  
  tt <- ipw(data.frame(data_comp), outcome=y, treat=treat,
            ps.method="glm")
  ipw.loglin.full <- rbind(ipw.loglin.full, cbind(tt, perc.miss, id, "glm"))
}

rows.cc <- apply(data_miss, 1, FUN = function(x) sum(is.na(x))==0)
try(dr.loglin.cc <- rbind(dr.loglin.cc, 
                          cbind(dr(data_miss[rows.cc, ], outcome=y[rows.cc ], treat=treat[rows.cc ],
                                   ps.method="glm",
                                   out.method = "glm"),
                                perc.miss, id, "glm")))

try(ipw.loglin.cc <- rbind(ipw.loglin.cc,
                           cbind(ipw(data_miss[rows.cc, ], outcome=y[rows.cc ], treat=treat[rows.cc ],
                                     ps.method="glm"),
                                 perc.miss, id, "glm")))

## FOREST
if (forest == "grf"){
  c.forest.miwae <- grf::causal_forest(data_zhat_miwae, y, treat)
  ate.miwae <- average_treatment_effect(c.forest.miwae, target.sample = "all")
  dr.grf.miwae <- rbind(dr.grf.miwae, cbind(t(ate.miwae), perc.miss, id, forest))
} else {
  df <- data.frame(cbind(data_zhat_miwae, "y" = y))
  rg.y <- ranger::ranger(y ~., data = df)
  df <- data.frame(cbind(data_zhat_miwae, "treat" = as.factor(treat)))
  rg.w <- ranger::ranger(treat ~., data = df, probability = TRUE)
  rg.y.hat <- predictions(rg.y)
  rg.w.hat <- predictions(rg.w)[,2]
  c.forest.miwae <- grf::causal_forest(data_zhat_miwae, y, treat, Y.hat = rg.y.hat, W.hat = rg.w.hat)
  ate.miwae <- average_treatment_effect(c.forest.miwae, target.sample = "all")
  dr.grf.miwae <- rbind(dr.grf.miwae, cbind(t(ate.miwae), perc.miss, id, forest))
}

tt <- ipw(data_zhat_miwae, outcome=y, treat=treat,
          ps.method=forest, trimming_weight = 0.999)
ipw.grf.miwae <- rbind(ipw.grf.miwae, cbind(tt, perc.miss, id, forest))

if (setting %in% c("mf", "dlvm2", "dlvm3")){
  if (forest == "grf"){
    c.forest.full = grf::causal_forest(Z, y, treat)
  } else {
    df <- data.frame(cbind(Z, "y" = y))
    rg.y <- ranger::ranger(y ~., data = df)
    df <- data.frame(cbind(Z, "treat" = as.factor(treat)))
    rg.w <- ranger::ranger(treat ~., data = df, probability = TRUE)
    rg.y.hat <- predictions(rg.y)
    rg.w.hat <- predictions(rg.w)[,2]
    c.forest.full <- grf::causal_forest(Z, y, treat, Y.hat = rg.y.hat, W.hat = rg.w.hat)
  }
  ate.full <- average_treatment_effect(c.forest.full, target.sample = "all")
  dr.grf.full <- rbind(dr.grf.full, cbind(t(ate.full), perc.miss, id, forest))
  
  tt <- ipw(data.frame(Z), outcome=y, treat=treat,
            ps.method=forest, trimming_weight = 0.999)
  ipw.grf.full <- rbind(ipw.grf.full, cbind(tt, perc.miss, id, forest))
} else if (setting %in% c("linear1", "linear2", "dlvm")){
  if (forest == "grf"){
    c.forest.full = grf::causal_forest(data_comp, y, treat)
  } else {
    df <- data.frame(cbind(data_comp, "y" = y))
    rg.y <- ranger::ranger(y ~., data = df)
    df <- data.frame(cbind(data_comp, "treat" = as.factor(treat)))
    rg.w <- ranger::ranger(treat ~., data = df, probability = TRUE)
    rg.y.hat <- predictions(rg.y)
    rg.w.hat <- predictions(rg.w)[,2]
    c.forest.full <- grf::causal_forest(data_comp, y, treat, Y.hat = rg.y.hat, W.hat = rg.w.hat)
  }
  ate.full <- average_treatment_effect(c.forest.full, target.sample = "all")
  dr.grf.full <- rbind(dr.grf.full, cbind(t(ate.full), perc.miss, id, forest))
  
  tt <- ipw(data.frame(data_comp), outcome=y, treat=treat,
            ps.method=forest, trimming_weight = 0.999)
  ipw.grf.full <- rbind(ipw.grf.full, cbind(tt, perc.miss, id, forest))
}

## Complete cases
if (forest == "grf"){
  try(c.forest.cc <- grf::causal_forest(data_miss[rows.cc, ], y[rows.cc ], treat[rows.cc ]))
} else {
  try(df <- data.frame(cbind(data_miss[rows.cc, ], "y" = y[rows.cc ])))
  try(rg.y <- ranger::ranger(y ~., data = df))
  try(df <- data.frame(cbind(data_miss[rows.cc, ], "w" = as.factor(treat[rows.cc ]))))
  try(rg.w <- ranger::ranger(w ~., data = df, probability=TRUE))
  try(rg.y.hat <- predictions(rg.y))
  try(rg.w.hat <- predictions(rg.w)[,2])
  try(c.forest.cc <- grf::causal_forest(data_miss[rows.cc, ], y[rows.cc ], treat[rows.cc ], Y.hat = rg.y.hat, W.hat = rg.w.hat))
}  
try(dr.grf.cc <- rbind(dr.grf.cc, 
                       cbind(t(average_treatment_effect(c.forest.cc, target.sample = "all")),
                             perc.miss, id, forest)))

try(ipw.grf.cc <- rbind(ipw.grf.cc,
                        cbind(ipw(data_miss[rows.cc, ], outcome=y[rows.cc ], treat=treat[rows.cc ],
                                  ps.method=forest, trimming_weight = 0.999),
                              perc.miss, id, forest)))

## Impute the data X with mice
mice.imp.woy <- prepare_data_ate(data_miss, treat, y, 
                                 imputation.method = "mice", mi.m=mi.m, use.outcome = F)

mice.imp.wy <- prepare_data_ate(data_miss, treat, y, 
                                imputation.method = "mice", mi.m=mi.m, use.outcome = T)

res <- c()
res.grf <- c()
for (k in 1:mi.m){
  try(res <- rbind(res,dr(X = mice.imp.woy$df.imp[[k]],
                          outcome = y, 
                          treat = treat, 
                          ps.method="glm", 
                          target= "all",
                          out.method = "glm")))
  if (forest == "grf"){
    try(c.forest <- grf::causal_forest(mice.imp.woy$df.imp[[k]], y, treat))
  } else {
    try(df <- data.frame(cbind(mice.imp.woy$df.imp[[k]], "y" = y)))
    try(rg.y <- ranger::ranger(y ~., data = df))
    try(df <- data.frame(cbind(mice.imp.woy$df.imp[[k]], "treat" = treat)))
    try(rg.w <- ranger::ranger(treat ~., data = df, probability = TRUE))
    try(rg.y.hat <- predictions(rg.y))
    try(rg.w.hat <- predictions(rg.w)[,2])
    try(c.forest <- grf::causal_forest(mice.imp.woy$df.imp[[k]], y, treat, Y.hat = rg.y.hat, W.hat = rg.w.hat))
  }  
  try(res.grf <- rbind(res.grf, t(average_treatment_effect(c.forest, target.sample = "all"))))
}
dr.loglin.mice.woy <- rbind(dr.loglin.mice.woy, cbind(t(apply(data.frame(res), FUN = mean, MARGIN = 2)),
                                                      perc.miss, id, "glm"))
dr.grf.mice.woy <- rbind(dr.grf.mice.woy, cbind(t(apply(data.frame(res.grf), FUN = mean, MARGIN = 2)),
                                                      perc.miss, id, forest))
res <- c()
res.grf <- c()
for (k in 1:mi.m){
  try(res <- rbind(res, ipw(X = mice.imp.woy$df.imp[[k]],
                            outcome = y, 
                            treat = treat, 
                            ps.method = "glm", 
                            target = "all")))
  try(res.grf <- rbind(res.grf, ipw(X = mice.imp.woy$df.imp[[k]],
                                    outcome = y, 
                                    treat = treat, 
                                    ps.method = forest, 
                                    trimming_weight = 0.999,
                                    target = "all")))
}
ipw.loglin.mice.woy <- rbind(ipw.loglin.mice.woy, cbind(t(apply(data.frame(res), FUN = mean, MARGIN = 2)),
                                                        perc.miss, id, "glm"))
ipw.grf.mice.woy <- rbind(ipw.loglin.mice.woy, cbind(t(apply(data.frame(res.grf), FUN = mean, MARGIN = 2)),
                                                        perc.miss, id, forest))

res <- c()
res.grf <- c()
for (k in 1:mi.m){
  try(res <- rbind(res,dr(X = mice.imp.wy$df.imp[[k]],
                          outcome = y, 
                          treat = treat, 
                          ps.method="glm", 
                          target= "all",
                          out.method = "glm")))
  if (forest == "grf"){
    try(c.forest <- grf::causal_forest(mice.imp.wy$df.imp[[k]], y, treat))
  } else {
    try(df <- data.frame(cbind(mice.imp.wy$df.imp[[k]], "y" = y)))
    try(rg.y <- ranger::ranger(y ~., data = df))
    try(df <- data.frame(cbind(mice.imp.wy$df.imp[[k]], "treat" = treat)))
    try(rg.w <- ranger::ranger(treat ~., data = df, probability = TRUE))
    try(rg.y.hat <- predictions(rg.y))
    try(rg.w.hat <- predictions(rg.w)[,2])
    try(c.forest <- grf::causal_forest(mice.imp.wy$df.imp[[k]], y, treat, Y.hat = rg.y.hat, W.hat = rg.w.hat))
  }  
  try(res.grf <- rbind(res.grf, t(average_treatment_effect(c.forest, target.sample = "all"))))
}
dr.loglin.mice.wy <- rbind(dr.loglin.mice.wy, cbind(t(apply(data.frame(res), FUN = mean, MARGIN = 2)),
                                                    perc.miss, id, "glm"))
dr.grf.mice.wy <- rbind(dr.grf.mice.wy, cbind(t(apply(data.frame(res.grf), FUN = mean, MARGIN = 2)),
                                                    perc.miss, id, forest))

res <- c()
res.grf <- c()
for (k in 1:mi.m){
  try(res <- rbind(res, ipw(X = mice.imp.wy$df.imp[[k]],
                            outcome = y, 
                            treat = treat, 
                            ps.method = "glm", 
                            target = "all")))
  try(res.grf <- rbind(res.grf, ipw(X = mice.imp.wy$df.imp[[k]],
                                    outcome = y, 
                                    treat = treat, 
                                    ps.method = forest,
                                    trimming_weight = 0.999, 
                                    target = "all")))
}
ipw.loglin.mice.wy <- rbind(ipw.loglin.mice.wy, cbind(t(apply(data.frame(res), FUN = mean, MARGIN = 2)),
                                                      perc.miss, id, "glm"))
ipw.grf.mice.wy <- rbind(ipw.grf.mice.wy, cbind(t(apply(data.frame(res.grf), FUN = mean, MARGIN = 2)),
                                                      perc.miss, id, forest))

## Estimate Z_hat using low-rank matrix factorization (Kallus, Mao & Udell, 2018)
udell.cv <- prepare_data_ate(data_miss, treat, y, 
                            imputation.method = "udell", use.outcome = F)
udell <- prepare_data_ate(data_miss, treat, y, 
                          imputation.method = "udell", use.outcome = F, r_udell = d)
# with cross-validation
tt <- dr(udell.cv$df.imp, outcome=y, treat=treat, 
         ps.method="glm",  
         out.method = "glm")
dr.loglin.mf.cv <- rbind(dr.loglin.mf.cv, cbind(tt, perc.miss, id, "glm"))

tt <- ipw(udell.cv$df.imp, outcome=y, treat=treat,
          ps.method="glm")
ipw.loglin.mf.cv <- rbind(ipw.loglin.mf.cv, cbind(tt, perc.miss, id, "glm"))

if (forest == "grf"){
  c.forest.mf = grf::causal_forest(udell.cv$df.imp, y, treat)
} else {
  df <- data.frame(cbind(udell.cv$df.imp, "y" = y))
  rg.y <- ranger::ranger(y ~., data = df)
  df <- data.frame(cbind(udell.cv$df.imp, "treat" = as.factor(treat)))
  rg.w <- ranger::ranger(treat ~., data = df, probability = TRUE)
  rg.y.hat <- predictions(rg.y)
  rg.w.hat <- predictions(rg.w)[,2]
  c.forest.mf <- grf::causal_forest(udell.cv$df.imp, y, treat, Y.hat = rg.y.hat, W.hat = rg.w.hat)
}
ate.mf <- average_treatment_effect(c.forest.mf, target.sample = "all")
dr.grf.mf.cv <- rbind(dr.grf.mf.cv, cbind(t(ate.mf), perc.miss, id, forest))

tt <- ipw(udell.cv$df.imp, outcome=y, treat=treat,
          ps.method=forest, trimming_weight = 0.999)
ipw.grf.mf.cv <- rbind(ipw.grf.mf.cv, cbind(tt, perc.miss, id, forest))

# without cross-validation
tt <- dr(udell$df.imp, outcome=y, treat=treat, 
         ps.method="glm",  
         out.method = "glm")
dr.loglin.mf <- rbind(dr.loglin.mf, cbind(tt, perc.miss, id, "glm"))

tt <- ipw(udell$df.imp, outcome=y, treat=treat,
          ps.method="glm")
ipw.loglin.mf <- rbind(ipw.loglin.mf, cbind(tt, perc.miss, id, "glm"))

if (forest == "grf"){
  c.forest.mf = grf::causal_forest(udell$df.imp, y, treat)
} else {
  df <- data.frame(cbind(udell$df.imp, "y" = y))
  rg.y <- ranger::ranger(y ~., data = df)
  df <- data.frame(cbind(udell$df.imp, "treat" = as.factor(treat)))
  rg.w <- ranger::ranger(treat ~., data = df, probability = TRUE)
  rg.y.hat <- predictions(rg.y)
  rg.w.hat <- predictions(rg.w)[,2]
  c.forest.mf <- grf::causal_forest(udell$df.imp, y, treat, Y.hat = rg.y.hat, W.hat = rg.w.hat)
}
ate.mf <- average_treatment_effect(c.forest.mf, target.sample = "all")
dr.grf.mf <- rbind(dr.grf.mf, cbind(t(ate.mf), perc.miss, id, forest))

tt <- ipw(udell$df.imp, outcome=y, treat=treat,
          ps.method=forest, trimming_weight = 0.999)
ipw.grf.mf <- rbind(ipw.grf.mf, cbind(tt, perc.miss, id, forest))


## MIA
mia <- prepare_data_ate(data_miss, treat, y, 
                        imputation.method = "mia.grf.ate", use.outcome = F)

tt <- ipw(mia$df.imp, outcome=y, treat=treat,
          ps.method=forest, trimming_weight = 0.999)
ipw.grf.mia <- rbind(ipw.grf.mia, cbind(tt, perc.miss, id, forest))

if (forest == "grf"){
  c.forest.mia = grf::causal_forest(mia$df.imp, y, treat)
} else {
  df <- data.frame(cbind(mia$df.imp, "y" = y))
  rg.y <- ranger::ranger(y ~., data = df)
  df <- data.frame(cbind(mia$df.imp, "treat" = as.factor(treat)))
  rg.w <- ranger::ranger(treat ~., data = df, probability = TRUE)
  rg.y.hat <- predictions(rg.y)
  rg.w.hat <- predictions(rg.w)[,2]
  c.forest.mia <- grf::causal_forest(mia$df.imp, y, treat, Y.hat = rg.y.hat, W.hat = rg.w.hat)
}
ate.mia <- average_treatment_effect(c.forest.mia, target.sample = "all")
dr.grf.mia <- rbind(dr.grf.mia, cbind(t(ate.mia), perc.miss, id, forest))




###########################################
# # Prepare and save first batch of results
###########################################

### Rename results
dr.loglin.miwae <- change_colnames(data.frame(dr.loglin.miwae), dr.col.names)
dr.loglin.miwae$dr <- as.double(as.character(dr.loglin.miwae$dr))
dr.loglin.miwae$se <- as.double(as.character(dr.loglin.miwae$se))

dr.loglin.full <-  change_colnames(data.frame(dr.loglin.full), dr.col.names)
dr.loglin.full$dr <- as.double(as.character(dr.loglin.full$dr))
dr.loglin.full$se <- as.double(as.character(dr.loglin.full$se))

try(dr.loglin.cc <- change_colnames(data.frame(dr.loglin.cc), dr.col.names))
try(dr.loglin.cc$dr <- as.double(as.character(dr.loglin.cc$dr)))
try(dr.loglin.cc$se <- as.double(as.character(dr.loglin.cc$se)))

dr.loglin.mice.woy <-  change_colnames(data.frame(dr.loglin.mice.woy), dr.col.names)
dr.loglin.mice.woy$dr <- as.double(as.character(dr.loglin.mice.woy$dr))
dr.loglin.mice.woy$se <- as.double(as.character(dr.loglin.mice.woy$se))

dr.loglin.mice.wy <- change_colnames(data.frame(dr.loglin.mice.wy), dr.col.names)
dr.loglin.mice.wy$dr <- as.double(as.character(dr.loglin.mice.wy$dr))
dr.loglin.mice.wy$se <- as.double(as.character(dr.loglin.mice.wy$se))

dr.loglin.mf <- change_colnames(data.frame(dr.loglin.mf), dr.col.names)
dr.loglin.mf$dr <- as.double(as.character(dr.loglin.mf$dr))
dr.loglin.mf$se <- as.double(as.character(dr.loglin.mf$se))
dr.loglin.mf.cv <- change_colnames(data.frame(dr.loglin.mf.cv), dr.col.names)
dr.loglin.mf.cv$dr <- as.double(as.character(dr.loglin.mf.cv$dr))
dr.loglin.mf.cv$se <- as.double(as.character(dr.loglin.mf.cv$se))

dr.grf.miwae <- change_colnames(data.frame(dr.grf.miwae), dr.col.names)
dr.grf.full <- change_colnames(data.frame(dr.grf.full), dr.col.names)
try(dr.grf.cc <- change_colnames(data.frame(dr.grf.cc), dr.col.names))
dr.grf.mice.woy <-  change_colnames(data.frame(dr.grf.mice.woy), dr.col.names)
dr.grf.mice.wy <- change_colnames(data.frame(dr.grf.mice.wy), dr.col.names)
dr.grf.mf <- change_colnames(data.frame(dr.grf.mf), dr.col.names)
dr.grf.mf.cv <- change_colnames(data.frame(dr.grf.mf.cv), dr.col.names)
dr.grf.mia <- change_colnames(data.frame(dr.grf.mia), dr.col.names)


ipw.loglin.miwae <- change_colnames(data.frame(ipw.loglin.miwae), ipw.col.names)
ipw.loglin.full <- change_colnames(data.frame(ipw.loglin.full), ipw.col.names)
try(ipw.loglin.cc <- change_colnames(data.frame(ipw.loglin.cc), ipw.col.names))
ipw.loglin.mice.woy <- change_colnames(data.frame(ipw.loglin.mice.woy), ipw.col.names)
ipw.loglin.mice.wy <- change_colnames(data.frame(ipw.loglin.mice.wy), ipw.col.names)
ipw.loglin.mf <- change_colnames(data.frame(ipw.loglin.mf), ipw.col.names)
ipw.loglin.mf.cv <- change_colnames(data.frame(ipw.loglin.mf.cv), ipw.col.names)

ipw.grf.full <- change_colnames(data.frame(ipw.grf.full), ipw.col.names)
ipw.grf.miwae <- change_colnames(data.frame(ipw.grf.miwae), ipw.col.names)
try(ipw.grf.cc <- change_colnames(data.frame(ipw.grf.cc), ipw.col.names))
ipw.grf.mice.woy <- change_colnames(data.frame(ipw.grf.mice.woy), ipw.col.names)
ipw.grf.mice.wy <- change_colnames(data.frame(ipw.grf.mice.wy), ipw.col.names)
ipw.grf.mia <- change_colnames(data.frame(ipw.grf.mia), ipw.col.names)
ipw.grf.mf <- change_colnames(data.frame(ipw.grf.mf), ipw.col.names)
ipw.grf.mf.cv <- change_colnames(data.frame(ipw.grf.mf.cv), ipw.col.names)

print(paste0('renamed first set of results for data set ', id))

### Aggregate results

dr.loglin <- rbind(cbind(dr.loglin.mice.woy, method = rep("loglin.mice.woy",dim(dr.loglin.mice.woy)[1])),
                   cbind(dr.loglin.mice.wy, method = rep("loglin.mice.wy",dim(dr.loglin.mice.wy)[1])),
                   cbind(dr.loglin.mf, method = rep("loglin.mf",dim(dr.loglin.mf)[1])),
                   cbind(dr.loglin.mf.cv, method = rep("loglin.mf.cv",dim(dr.loglin.mf.cv)[1])),
                   cbind(dr.loglin.miwae, method = rep("loglin.miwae",dim(dr.loglin.miwae)[1])),
                   cbind(dr.loglin.full, method = rep("loglin.Z_true",dim(dr.loglin.full)[1])),
                   cbind(dr.grf.mice.woy, method = rep("grf.mice.woy",dim(dr.grf.mice.woy)[1])),
                   cbind(dr.grf.mice.wy, method = rep("grf.mice.wy",dim(dr.grf.mice.wy)[1])),
                   cbind(dr.grf.miwae, method = rep("grf.miwae",dim(dr.grf.miwae)[1])),
                   cbind(dr.grf.full, method = rep("grf.Z_true",dim(dr.grf.full)[1])),
                   cbind(dr.grf.mf, method = rep("grf.mf",dim(dr.grf.mf)[1])),
                   cbind(dr.grf.mf.cv, method = rep("grf.mf.cv",dim(dr.grf.mf.cv)[1])),
                   cbind(dr.grf.mia, method = rep("grf.mia",dim(dr.grf.mia)[1])))
print(paste0("first dr rbind passed for data set ", id))

try(dr.loglin <- rbind(dr.loglin, 
                       cbind(dr.grf.cc, method = rep("grf.cc",dim(dr.grf.cc)[1])),
                       cbind(dr.loglin.cc, method = rep("loglin.cc",dim(dr.loglin.cc)[1]))))

print(paste0("second dr rbind passed for data set ", id))

ipw.loglin <- rbind(cbind(data.frame(ipw.loglin.mice.woy), method = rep("loglin.mice.woy",dim(ipw.loglin.mice.woy)[1])),
                    cbind(data.frame(ipw.loglin.mice.wy), method = rep("loglin.mice.wy",dim(ipw.loglin.mice.wy)[1])),
                    cbind(data.frame(ipw.loglin.mf), method = rep("loglin.mf",dim(ipw.loglin.mf)[1])),
                    cbind(data.frame(ipw.loglin.mf.cv), method = rep("loglin.mf.cv",dim(ipw.loglin.mf.cv)[1])),
                    cbind(data.frame(ipw.loglin.miwae), method = rep("loglin.miwae",dim(ipw.loglin.miwae)[1])),
                    cbind(data.frame(ipw.loglin.full), method = rep("loglin.Z_true",dim(ipw.loglin.full)[1])),
                    cbind(data.frame(ipw.grf.mice.woy), method = rep("grf.mice.woy",dim(ipw.grf.mice.woy)[1])),
                    cbind(data.frame(ipw.grf.mice.wy), method = rep("grf.mice.wy",dim(ipw.grf.mice.wy)[1])),
                    cbind(data.frame(ipw.grf.full), method = rep("grf.Z_true",dim(ipw.grf.full)[1])),
                    cbind(data.frame(ipw.grf.mia), method = rep("grf.mia",dim(ipw.grf.mia)[1])),
                    cbind(data.frame(ipw.grf.miwae), method = rep("grf.miwae",dim(ipw.grf.miwae)[1])),
                    cbind(data.frame(ipw.grf.mf), method = rep("grf.mf",dim(ipw.grf.mf)[1])),
                    cbind(data.frame(ipw.grf.mf.cv), method = rep("grf.mf.cv",dim(ipw.grf.mf.cv)[1])))
print(paste0("first ipw rbind passed for data set ", id))

try(ipw.loglin <- rbind(ipw.loglin,
                        cbind(data.frame(ipw.loglin.cc), method = rep("loglin.cc",dim(ipw.loglin.cc)[1])),
                        cbind(data.frame(ipw.grf.cc), method = rep("grf.cc",dim(ipw.grf.cc)[1]))))

print(paste0("second ipw rbind passed for data set ", id))


### Save results
if (use.w){
  save(dr.loglin, ipw.loglin, file = paste0(rdata.dir, opt$out_date, 
                                          "_",link,"_",setting, 
                                          "_n", n, "_p", p, "_d", d, 
                                          "_W_dmiwae", d.miwae, "_hmiwae", h.miwae,
                                          "_",add_mask,
                                          "_sir", num.samples.sir,
                                          "_propNA", format(round(perc.miss, 2), nsmall=2),
                                          "_cit2", as.logical(cit2),
                                          "_ci2Imp", ci2_imp,
                                          "_forest", forest,
                                          "_results_dr_ipw_seed",id,
                                          ".RData"))
} else {
  save(dr.loglin, ipw.loglin, file = paste0(rdata.dir, opt$out_date, 
                                            "_",link,"_",setting, 
                                            "_n", n, "_p", p, "_d", d, 
                                            "_dmiwae", d.miwae, "_hmiwae", h.miwae,
                                            "_",add_mask,
                                            "_sir", num.samples.sir,
                                            "_propNA", format(round(perc.miss, 2), nsmall=2),
                                            "_cit2", as.logical(cit2),
                                            "_ci2Imp", ci2_imp,
                                            "_forest", forest,
                                            "_results_dr_ipw_seed",id,
                                            ".RData"))
}

############################################
# # Add estimation on sampled Z.hat
# # (multiple imputation and e^* estimation)
############################################

w.hat.grf.agg <- y.hat.grf.agg <- c()
ipw.grf.agg <- c()
dr.grf.agg <- c()
w.hat.log.agg <- y_1.hat.lin.agg <- y_0.hat.lin.agg <- c()
ipw.loglin.agg <- c()
dr.loglin.agg <- c()


w.hat.grf.mul <- y.hat.grf.mul <- c()
w.hat.log.mul <- y_1.hat.lin.mul <- y_0.hat.lin.mul <- c()
for (k in 0:(num.samples.sir-1)){
  cat(paste(k," "))
  # if (d.miwae == 3){
  #   if (use.w){
  #     data_zhat_miwae <- read.csv(file = paste0(data.dir, setting, "_", n, "_", p, "_", d, 
  #                                               "_W_hmiwae", h.miwae, add_mask, 
  #                                               "_seed", id, 
  #                                               "_propNA", format(round(perc.miss, 2), nsmall=2),
  #                                               "_zhat_m",k,".csv"), 
  #                                 sep = ";", header=F)
  #   } else {
  #     data_zhat_miwae <- read.csv(file = paste0(data.dir, setting, "_", n, "_", p, "_", d, 
  #                                               "_hmiwae", h.miwae, add_mask, 
  #                                               "_seed", id, 
  #                                               "_propNA", format(round(perc.miss, 2), nsmall=2),
  #                                               "_zhat_m",k,".csv"), 
  #                                 sep = ";", header=F)
  #   }
  # } else {
    if (use.w){
      data_zhat_miwae <- read.csv(file = paste0(data.dir, setting, "_", n, "_", p, "_", d, 
                                                "_W_dmiwae", d.miwae, "_hmiwae", h.miwae, add_mask, 
                                                "_seed", id, 
                                                "_propNA", format(round(perc.miss, 2), nsmall=2),
                                                "_zhat_m",k,".csv"), 
                                  sep = ";", header=F)
    } else {
      data_zhat_miwae <- read.csv(file = paste0(data.dir, setting, "_", n, "_", p, "_", d, 
                                                "_dmiwae", d.miwae, "_hmiwae", h.miwae, add_mask, 
                                                "_seed", id, 
                                                "_propNA", format(round(perc.miss, 2), nsmall=2),
                                                "_zhat_m",k,".csv"), 
                                  sep = ";", header=F)
    }
  # }
  if (dim(data_zhat_miwae)[2] > d.miwae){
    data_zhat_miwae <- data_zhat_miwae[,-1]
  }

  ## grf mul
  if (forest == "grf"){
    Z.m = model.matrix(~. , data=data_zhat_miwae)
    forest.W = regression_forest(Z.m, treat, tune.parameters = TRUE)
    w.hat = predict(forest.W)$predictions
    forest.Y = regression_forest(Z.m, y, tune.parameters = TRUE)
    y.hat = predict(forest.Y)$predictions
  } else {
    df <- data.frame(cbind(data_zhat_miwae, "y" = y))
    rg.y <- ranger::ranger(y ~., data = df)
    df <- data.frame(cbind(data_zhat_miwae, "treat" = as.factor(treat)))
    rg.w <- ranger::ranger(treat ~., data = df, probability = TRUE)
    y.hat <- predictions(rg.y)
    w.hat <- predictions(rg.w)[,2]
  }
  w.hat.grf.mul <- rbind(w.hat.grf.mul, w.hat)
  y.hat.grf.mul <- rbind(y.hat.grf.mul, y.hat)

  ## grf agg (estimate ATE for each sample Z|X and aggregate the taus)
  c.forest.miwae = grf::causal_forest(data_zhat_miwae, y, treat,
                                    W.hat = w.hat,
                                    Y.hat = y.hat)
  ate.miwae <- average_treatment_effect(c.forest.miwae, target.sample = "all")
  dr.grf.agg <- rbind(dr.grf.agg, t(ate.miwae))

  w.hat <- (treat)/w.hat+ (1 - treat)/ (1 - w.hat)
  if (any(is.nan(w.hat)) | any(is.infinite(w.hat))){
    weightMax <- quantile(w.hat, 0.999, na.rm=TRUE) 
    w.hat[is.na(w.hat)] <- weightMax # If the weight is NA this means that pscore was 0
    w.hat[w.hat > weightMax] <- weightMax
    weightMin <- quantile(w.hat, 0.001) 
    w.hat[w.hat < weightMin] <- weightMin
  }
  ipw1 <- 1/n*(sum(y[which(treat==1)]*w.hat[which(treat==max(treat))]) - sum(y[which(!(treat==max(treat)))]*w.hat[which(!(treat==max(treat)))]))
  ipw2 <- 1/sum(w.hat[which(treat==max(treat))]) * sum(y[which(treat==max(treat))] * w.hat[which(treat==max(treat))]) -
    1/sum(w.hat[which(!(treat==max(treat)))]) * sum(y[which(!(treat==max(treat)))] * w.hat[which(!(treat==max(treat)))])
  mod <- lm(y~treat, weights = w.hat)
  se.ipw2 <- sqrt(diag(sandwich::vcovHC(mod, type = "HC")))[2]
  tt <- cbind(ipw1, ipw2, se.ipw2)
  ipw.grf.agg <- rbind(ipw.grf.agg, tt)
  
  ## log-lin mul
  w.hat <- predict_glm(data_zhat_miwae, as.factor(treat), seed=id)$pscore
  w.hat.log.mul <- rbind(w.hat.log.mul, t(w.hat))
  tmp <- fit_outcome_glm(data_zhat_miwae, y, treat)
  y_1.hat.lin.mul <- rbind(y_1.hat.lin.mul, tmp[,1])
  y_0.hat.lin.mul <- rbind(y_0.hat.lin.mul, tmp[,2])

  ## log-lin agg (estimate ATE for each sample Z|X and aggregate the taus)
  w.hat <- (treat)/w.hat + (1 - treat)/ (1 - w.hat)
  delta_i <- tmp[,1] -  tmp[,2] + treat*(y-tmp[,1])*w.hat - (1-treat)*(y-tmp[,2])*w.hat
  dr.loglin.agg <- rbind(dr.loglin.agg,
                         cbind(mean(delta_i), sqrt(sum((delta_i-mean(delta_i))^2))/n))

  ipw1 <- 1/n*(sum(y[which(treat==1)]*w.hat[which(treat==max(treat))]) - sum(y[which(!(treat==max(treat)))]*w.hat[which(!(treat==max(treat)))]))
  ipw2 <- 1/sum(w.hat[which(treat==max(treat))]) * sum(y[which(treat==max(treat))] * w.hat[which(treat==max(treat))]) -
  1/sum(w.hat[which(!(treat==max(treat)))]) * sum(y[which(!(treat==max(treat)))] * w.hat[which(!(treat==max(treat)))])
  mod <- lm(y~treat, weights = w.hat)
  se.ipw2 <- sqrt(diag(sandwich::vcovHC(mod, type = "HC")))[2]
  tt <- cbind(ipw1, ipw2, se.ipw2)
  ipw.loglin.agg <- rbind(ipw.loglin.agg, tt)

  
}
## log-lin
w.hat.log.mul <- apply(w.hat.log.mul, 2, mean)
y_1.hat.lin.mul <- apply(y_1.hat.lin.mul, 2, mean)
y_0.hat.lin.mul <- apply(y_0.hat.lin.mul, 2, mean)

w.hat.log.mul <- (treat)/w.hat.log.mul + (1 - treat)/ (1 - w.hat.log.mul)
delta_i <- y_1.hat.lin.mul -  y_0.hat.lin.mul + 
  treat*(y-y_1.hat.lin.mul)*w.hat.log.mul - (1-treat)*(y-y_0.hat.lin.mul)*w.hat.log.mul
dr.loglin.miwae.zmul <- rbind(dr.loglin.miwae.zmul,
                              cbind(mean(delta_i), sqrt(sum((delta_i-mean(delta_i))^2))/n, 
                                    opt$percmiss, opt$id, "glm"))

ipw1 <- 1/n*(sum(y[which(treat==1)]*w.hat.log.mul[which(treat==max(treat))]) - sum(y[which(!(treat==max(treat)))]*w.hat.log.mul[which(!(treat==max(treat)))]))
ipw2 <- 1/sum(w.hat.log.mul[which(treat==max(treat))]) * sum(y[which(treat==max(treat))] * w.hat.log.mul[which(treat==max(treat))]) -
  1/sum(w.hat.log.mul[which(!(treat==max(treat)))]) * sum(y[which(!(treat==max(treat)))] * w.hat.log.mul[which(!(treat==max(treat)))])
mod <- lm(y~treat, weights = w.hat.log.mul)
se.ipw2 <- sqrt(diag(sandwich::vcovHC(mod, type = "HC")))[2]
tt <- cbind(ipw1, ipw2, se.ipw2)
ipw.loglin.miwae.zmul <- rbind(ipw.loglin.miwae.zmul, cbind(tt, perc.miss, id, "glm"))

## grf
w.hat.grf.mul <- apply(w.hat.grf.mul, 2, mean)
y.hat.grf.mul <- apply(y.hat.grf.mul, 2, mean)
c.forest.miwae = grf::causal_forest(data_zhat_miwae, y, treat,
                                    W.hat = w.hat.grf.mul,
                                    Y.hat = y.hat.grf.mul)
ate.miwae <- average_treatment_effect(c.forest.miwae, target.sample = "all")
dr.grf.miwae.zmul <- rbind(dr.grf.miwae.zmul, cbind(t(ate.miwae), perc.miss, id, forest))

w.hat.grf.mul <- (treat)/w.hat.grf.mul + (1 - treat)/ (1 - w.hat.grf.mul)
if (any(is.nan(w.hat.grf.mul)) | any(is.infinite(w.hat.grf.mul))){
  weightMax <- quantile(w.hat.grf.mul, 0.999, na.rm=TRUE) 
  w.hat.grf.mul[is.na(w.hat.grf.mul)] <- weightMax # If the weight is NA this means that pscore was 0
  w.hat.grf.mul[w.hat.grf.mul > weightMax] <- weightMax
  weightMin <- quantile(w.hat.grf.mul, 0.001) 
  w.hat.grf.mul[w.hat.grf.mul < weightMin] <- weightMin
}
ipw1 <- 1/n*(sum(y[which(treat==1)]*w.hat.grf.mul[which(treat==max(treat))]) - sum(y[which(!(treat==max(treat)))]*w.hat.grf.mul[which(!(treat==max(treat)))]))
ipw2 <- 1/sum(w.hat.grf.mul[which(treat==max(treat))]) * sum(y[which(treat==max(treat))] * w.hat.grf.mul[which(treat==max(treat))]) -
  1/sum(w.hat.grf.mul[which(!(treat==max(treat)))]) * sum(y[which(!(treat==max(treat)))] * w.hat.grf.mul[which(!(treat==max(treat)))])
mod <- lm(y~treat, weights = w.hat.grf.mul)
se.ipw2 <- sqrt(diag(sandwich::vcovHC(mod, type = "HC")))[2]
tt <- cbind(ipw1, ipw2, se.ipw2)
ipw.grf.miwae.zmul <- rbind(ipw.grf.miwae.zmul, cbind(tt, perc.miss, id, forest))



### Multiple imputation strategy: aggregate tau estimates of 
dr.grf.miwae.agg <- rbind(dr.grf.miwae.agg, cbind(t(apply(dr.grf.agg[,1:2],2,mean)), perc.miss, id, forest))
dr.loglin.miwae.agg <- rbind(dr.loglin.miwae.agg, cbind(t(apply(dr.loglin.agg[,1:2],2,mean)), perc.miss, id, "glm"))
ipw.grf.miwae.agg <- rbind(ipw.grf.miwae.agg, cbind(t(apply(ipw.grf.agg[,1:3],2,mean)), perc.miss, id, forest))
ipw.loglin.miwae.agg <- rbind(ipw.loglin.miwae.agg, cbind(t(apply(ipw.loglin.agg[,1:3],2,mean)), perc.miss, id, "glm"))


#############################################
# # Prepare, aggregate and save final results
#############################################


dr.grf.miwae.agg <- change_colnames(dr.grf.miwae.agg, dr.col.names)
dr.loglin.miwae.agg <- change_colnames(dr.loglin.miwae.agg, dr.col.names)
ipw.grf.miwae.agg <- change_colnames(ipw.grf.miwae.agg, ipw.col.names)
ipw.loglin.miwae.agg <- change_colnames(ipw.loglin.miwae.agg, ipw.col.names)

dr.grf.miwae.zmul <- change_colnames(dr.grf.miwae.zmul, dr.col.names)
dr.loglin.miwae.zmul <- change_colnames(dr.loglin.miwae.zmul, dr.col.names)
ipw.grf.miwae.zmul <- change_colnames(ipw.grf.miwae.zmul, ipw.col.names)
ipw.loglin.miwae.zmul <- change_colnames(ipw.loglin.miwae.zmul, ipw.col.names)


print(paste0("finished estimation for data set ", id))


try(dr.loglin <- rbind(dr.loglin,
                       cbind(data.frame(dr.grf.miwae.zmul), method = rep("grf.miwae.zmul",dim(dr.grf.miwae.zmul)[1])),
                       cbind(data.frame(dr.loglin.miwae.zmul), method = rep("loglin.miwae.zmul",dim(dr.loglin.miwae.zmul)[1]))))
print(paste0("third dr rbind passed for data set ", id))

try(dr.loglin <- rbind(dr.loglin,
                       cbind(data.frame(dr.grf.miwae.agg), method = rep("grf.miwae.agg",dim(dr.grf.miwae.agg)[1])),
                       cbind(data.frame(dr.loglin.miwae.agg), method = rep("loglin.miwae.agg",dim(dr.loglin.miwae.agg)[1]))))
print(paste0("fourth dr rbind passed for data set ", id))


try(ipw.loglin <- rbind(ipw.loglin,
                        cbind(data.frame(ipw.grf.miwae.zmul), method = rep("grf.miwae.zmul",dim(ipw.grf.miwae.zmul)[1])),
                        cbind(data.frame(ipw.loglin.miwae.zmul), method = rep("loglin.miwae.zmul",dim(ipw.loglin.miwae.zmul)[1]))))
print(paste0("third ipw rbind passed for data set ", id))

try(ipw.loglin <- rbind(ipw.loglin,
                        cbind(data.frame(ipw.grf.miwae.agg), method = rep("grf.miwae.agg",dim(ipw.grf.miwae.agg)[1])),
                        cbind(data.frame(ipw.loglin.miwae.agg), method = rep("loglin.miwae.agg",dim(ipw.loglin.miwae.agg)[1]))))
print(paste0("fourth ipw rbind passed for data set ", id))

### Final save results
if (use.w){
  save(dr.loglin, ipw.loglin, file = paste0(rdata.dir, opt$out_date, 
                                          "_",link,"_",setting, 
                                          "_n", n, "_p", p, "_d", d, 
                                          "_W_dmiwae", d.miwae, "_hmiwae", h.miwae,
                                          "_",add_mask,
                                          "_sir", num.samples.sir,
                                          "_propNA", format(round(perc.miss, 2), nsmall=2),
                                          "_cit2", as.logical(cit2),
                                          "_ci2Imp", ci2_imp,
                                          "_forest", forest,
                                          "_results_dr_ipw_seed",id,
                                          ".RData"))
} else {
  save(dr.loglin, ipw.loglin, file = paste0(rdata.dir, opt$out_date, 
                                            "_",link,"_",setting, 
                                            "_n", n, "_p", p, "_d", d, 
                                            "_dmiwae", d.miwae, "_hmiwae", h.miwae,
                                            "_",add_mask,
                                            "_sir", num.samples.sir,
                                            "_propNA", format(round(perc.miss, 2), nsmall=2),
                                            "_cit2", as.logical(cit2),
                                            "_ci2Imp", ci2_imp,
                                            "_forest", forest,
                                            "_results_dr_ipw_seed",id,
                                            ".RData"))




