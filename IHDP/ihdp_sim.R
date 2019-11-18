library(MASS) # mvrnorm
library(caret) # cross-validated training
library(randomForest) # Random Forests
library(ranger) # Random Forests prediction
library(FactoMineR) # Factorial data analysis
library(grf) # Doubly robust treatment effect estimation
library(cobalt)

library(misaem) # Logistic regression with missing values
library(missMDA) # PCA/MCA with missing values + iterative PC imputation
library(missForest) # Random forests with missing values
library(mice) # Multiple imputation

library(foreach)
library(doParallel)


library(assertthat)
library(norm)
library(softImpute)

library(dplyr)
library(tidyr)
library(directlabels)
min_inf_proxy <- -1e300
max_inf_proxy <- 1e300

setwd("~/Documents/TraumaMatrix/CausalInference/Simulations/causal-inference-missing/IHDP/")

source("./vdorie-npci/data.R")
source("./vdorie-npci/util.R")
source("./vdorie-npci/transform.R")

source("../Helper/helper_causalInference.R")
source("../Helper/helper_simulations.R")
source("../Helper/helper_imputation.R")
source("../Helper/helper_udell.R")
source("../Utils/miss.saem.v2.R")
source("../Utils/amputation.R")

loadDataInCurrentEnvironment(covariates = "select", p.score = "none")

x <- as.matrix(x)
w <- rep(0.5, ncol(x))

#######################################################################################
##### Complete
#######################################################################################

N <- 10

methods <- c("glm", "grf.ate")
ate_ipw <- matrix(NA, nrow = N, ncol = 2*length(methods))
ate_dr <- matrix(NA, nrow = N, ncol = 2*length(methods))

for (seed in 1:N){
  
  generateDataForIterInCurrentEnvironment(seed, x, z, w, overlap = "all")
  
  x.r <- data.frame(x.r)
  vars.binary <- colnames(x.r)[apply(x.r, 2, FUN = function(e) length(unique(e))==2)]
  for (j in vars.binary){
    x.r[,j] <- as.factor(x.r[,j])
  }
  

  for (m in 1:length(methods)){
    if (!(methods[m] %in% c("grf.ate"))){
      ate_ipw[seed,(2*m-1):(2*m)] <- ipw(X=x.r, outcome = y, treat = z,
                                         ps.method = methods[m], seed = seed)
    }
    res <- dr(X=x.r, 
              outcome = y, 
              treat = z,
              ps.method= methods[m],  
              target= "all", 
              seed = seed,
              out.method = methods[m])
    ate_dr[seed,(2*m-1):(2*m)] <-  c(res$dr, res$se)
    
  }
}

colnames(ate_ipw) <- kronecker(methods, c("ate1","ate2"),function(x,r) paste0(r,"_",x))
colnames(ate_dr) <- kronecker(methods, c("ate","sd.ate"),function(x,r) paste0(r,"_",x))


#######################################################################################
##### MISSING VALUES
#######################################################################################
imp.methods <- c("mean", "mean.grf.ate", "mia.grf.ate","saem", "udell")
prop.missing <- 0.7
res_ate_ipw_miss <- c()
res_ate_dr_miss <- c()


mech <- "MNAR"
for (cit in c(F,T)){
  for (cio in c(F,T)){
    for (use.mask in c(F,T)){
      writeLines(paste0("cit: ", cit,", cio: ", cio, ", use.mask: ", use.mask))
      for (seed in 1:10){
        cat(seed, " ")
        
        generateDataForIterInCurrentEnvironment(seed, x, z, w, overlap = "all")
        x.r <- data.frame(x.r)
        
        #patterns <- t(rmultinom(dim(x.r)[1], size=floor(dim(x.r)[2]/3), prob=rep(1/dim(x.r)[2],dim(x.r)[2])))
        #patterns <- apply(patterns, c(1,2), FUN = function(x) as.integer(x>0))
        #x.miss <- mice::ampute(x.r, prop = prop.missing, patterns = patterns, freq = rep(1/dim(x.r)[1], dim(x.r)[1]), mech=mech, std = F)$amp
        x.miss <- produce_NA(x.r, mechanism = mech, perc.missing = prop.missing)$data.incomp
        
        vars.binary <- colnames(x)[apply(x.r, 2, FUN = function(e) length(unique(e))==2)]
        for (j in vars.binary){
          x.miss[,j] <- as.factor(x.miss[,j])
        }
        
        idx_incomp <- which(apply(x.miss[,which(sapply(x.miss, is.numeric))], MARGIN = 1, FUN = function(elt) sum(is.na(elt))==length(elt)))
        p.num <- sum(sapply(x.miss, is.numeric))
        for (i in idx_incomp){
          idx_col <- rbinom(1, p.num-1, 0.5)+1
          x.miss[i, idx_col] <- x.r[i, idx_col]
        }
        
        generateDataForIterInCurrentEnvironment(seed, x, z, w, mask = is.na(x.miss),
                                                cit=cit,cio=cio,overlap = "all")
        
        for (imp in imp.methods){
          if (imp %in% c("saem","udell")){
            x.imp.fitted <- prepare_data_ate(df=x.miss[,which(sapply(x.miss, is.numeric))], w=z.r, y=y, 
                                             imputation.method=imp, mi.m = 10,
                                             mask = NULL,
                                             use.outcome=F, use.interaction=F)
            x.miss <- x.miss[,which(sapply(x.miss, is.numeric))]
          } else{
            x.imp.fitted <- prepare_data_ate(df=x.miss, w=z.r, y=y, 
                                             imputation.method=imp, mi.m = 10,
                                             mask = NULL,
                                             use.outcome=F, use.interaction=F)
          }
          x.imp <- x.imp.fitted$df.imp
          fitted <- x.imp.fitted$fitted
          
          
          if (imp %in% c("mean", "saem", "pc", "udell")){ 
            ps.method <- "glm"
            out.method <- "glm"
          } 
          if (imp %in% c("mia.grf.ate", "mean.grf.ate")) {
            ps.method <- "grf.ate"
            out.method <- "grf.ate"
          }
          
          mask <- NULL
          if (use.mask){
            mask <- is.na(x.miss)
          }
          
          res <- ipw(X=data.frame(x.imp), outcome = y, treat = z.r,
                     ps.method = ps.method, seed = seed,
                     mask = mask,
                     fitted = fitted)
          res_ate_ipw_miss <- rbind(res_ate_ipw_miss,
                                c(res[1], NA,"ipw1", imp, mech, use.mask, paste0(cit,".", cio), cit, cio))
          res_ate_ipw_miss <- rbind(res_ate_ipw_miss,
                                c(res[2], NA, "ipw2", imp, mech, use.mask, paste0(cit,".", cio), cit, cio))
          
          res <- data.frame(dr = c(NA,NA))
          
          try(res <- dr(X=data.frame(x.imp), outcome = y, treat = z.r,
                        ps.method= ps.method,  
                        target= "all", 
                        seed = seed,
                        mask = mask,
                        out.method = out.method,
                        fitted = fitted))
          res_ate_dr_miss <- rbind(res_ate_dr_miss, 
                               c(res$dr, res$se, "dr", imp, mech, use.mask, paste0(cit,".", cio), cit, cio))
          
        }
        cat("\n")
      }
    }
  }
}
res_ate_ipw_miss <- data.frame(res_ate_ipw_miss)
res_ate_dr_miss <- data.frame(res_ate_dr_miss)

colnames(res_ate_ipw_miss) <- c("ate", "se.ate", "type", "method", "mechanism", "use.mask", "citcio", "cit", "cio")
colnames(res_ate_dr_miss) <- c("ate", "se.ate", "type", "method", "mechanism", "use.mask", "citcio", "cit", "cio")

res_ate_dr_miss$ate <- as.numeric(levels(res_ate_dr_miss$ate))[res_ate_dr_miss$ate]
res_ate_dr_miss$se.ate <- as.numeric(levels(res_ate_dr_miss$se.ate))[res_ate_dr_miss$se.ate]
res_ate_ipw_miss$ate <- as.numeric(levels(res_ate_ipw_miss$ate))[res_ate_ipw_miss$ate]


ate_dr_miss <- data.frame(res_ate_dr_miss) %>%
                group_by(mechanism, type, method, use.mask, citcio) %>%
                summarise(ate.mean = mean(ate)) %>%
                ungroup()

ate_ipw_miss <- data.frame(res_ate_ipw_miss) %>%
                  group_by(mechanism, type, method, use.mask, citcio) %>%
                  summarise(ate.mean = mean(ate)) %>%
                  ungroup()

#res_dr_mar02 <- res_ate_dr_miss
#res_ipw_mar02 <- res_ate_ipw_miss
#res_dr_mar05 <- res_ate_dr_miss
#res_ipw_mar05 <- res_ate_ipw_miss
# res_dr_mar07 <- res_ate_dr_miss
# res_ipw_mar07 <- res_ate_ipw_miss
res_mar <- rbind(res_ipw_mar02, res_ipw_mar05, res_ipw_mar07,
                 res_dr_mar02, res_dr_mar05, res_dr_mar07)

#res_dr_mnar02 <- res_ate_dr_miss
#res_ipw_mnar02 <- res_ate_ipw_miss
#res_dr_mnar05 <- res_ate_dr_miss
#res_ipw_mnar05 <- res_ate_ipw_miss
#res_dr_mnar07 <- res_ate_dr_miss
#res_ipw_mnar07 <- res_ate_ipw_miss
res_mnar <- rbind(res_ipw_mnar02, res_ipw_mnar05, res_ipw_mnar07,
                  res_dr_mnar02, res_dr_mnar05, res_dr_mnar07)

save(res_mar, res_mnar,
     file=paste0("./RData/", Sys.date(), "_ihdp_N", N, "_results_compare_miss.RData"))

##################
### PLOTS
##################

load(file="./RData/2019-06-14_ihdp_N10_results_compare_miss.RData")

levels(res_mar$method)[which(levels(res_mar$method)=="mean")] <- "mean.loglin"
levels(res_mnar$method)[which(levels(res_mnar$method)=="mean")] <- "mean.loglin"

levels(res_mar$method)[which(levels(res_mar$method)=="udell")] <- "mf"
levels(res_mnar$method)[which(levels(res_mnar$method)=="udell")] <- "mf"

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2", "#D55E00",  "#F0E442")
names(cbPalette) <- c("mean.loglin", "mia.grf.ate", "mice.woy", "mice.wy", "saem", "mf", "mean.grf.ate", "mice")


ate_mar <- data.frame(res_mar) %>%
  group_by(mechanism, prop.missing, type, method, use.mask, citcio) %>%
  summarise(ate.mean = mean(ate), ate.se = mean(se.ate)) %>%
  ungroup()

ate_mnar <- data.frame(res_mnar) %>%
  group_by(mechanism, prop.missing, type, method, use.mask, citcio) %>%
  summarise(ate.mean = mean(ate), ate.se = mean(se.ate)) %>%
  ungroup()

df <-ate_mar[which(ate_mar$method != "mf"),]
plot_ate_bias <- function(df, mechanism){
  
  plot_ipw <- ggplot(data=df[which(df$type == "ipw2"),],
                     aes(x=prop.missing, y=abs(4-ate.mean), color = method, linetype = use.mask)) +
    geom_line() +
    geom_point() +
    #geom_errorbar(aes(ymin=ate-sd, ymax=ate+sd), width=.1) +
    geom_hline(aes(yintercept=0), colour="black") +
    ylab("Bias (absolute value)") +
    #scale_y_continuous(limits = c(0.99*min_ate, 1.01*max_ate)) +
    ggtitle(paste0(mechanism, ", IPW")) +
    facet_wrap(~citcio, nrow=2, labeller = label_both) +
    scale_color_manual(values=cbPalette,
                       name = "Method", labels = sort(unlist(unique(df[which(df$type == "ipw2"), "method"]))))
  
  
  plot_dr <- ggplot(data=df[which(df$type=="dr"),],
                    aes(x=prop.missing, y=abs(4-ate.mean), color = method, linetype = use.mask)) +
    geom_line() +
    geom_point() +
    #geom_errorbar(aes(ymin=ate-sd, ymax=ate+sd), width=.1) +
    geom_hline(aes(yintercept=0), colour="black") +
    ylab("Bias (absolute value)") +
    #scale_y_continuous(limits = c(0.99*min_ate, 1.01*max_ate)) +
    ggtitle(paste0(mechanism, ", DR")) +
    facet_wrap(~citcio, nrow=2, labeller = label_both) +
    scale_color_manual(values=cbPalette,
                       name = "Method", labels = sort(unlist(unique(df[which(ate_mar$type=="dr"),"method"]))))
  
  return(list(plot_ipw = plot_ipw, plot_dr = plot_dr))
}

df <-ate_mar[which(ate_mar$method != "mf"),]
plots_mar <- plot_ate_bias(df, mechanism = "MAR")
plots_mar$plot_ipw
plots_mar$plot_dr

df <-ate_mnar[which(ate_mnar$method != "mf"),]
plots_mnar <- plot_ate_bias(df, mechanism = "MNAR")
plots_mnar$plot_ipw
plots_mnar$plot_dr
