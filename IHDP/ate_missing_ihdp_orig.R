################################################################
# Example script to run a simple experiment with synthetic data 
# to compare different ATE estimators that handle missing values
################################################################

#script.dir <- dirname(sys.frame(1)$ofile)
#setwd(script.dir)
setwd("~/cluster_files/IHDP/")
#setwd("~/Documents/TraumaMatrix/CausalInference/Simulations/causal-inference-missing/IHDP/")
source("../load.R")

source("./vdorie-npci/data.R")
source("./vdorie-npci/util.R")
source("./vdorie-npci/transform.R")



loadDataInCurrentEnvironment(covariates = "select", p.score = "none")

x <- as.matrix(x)
w <- rep(0.5, ncol(x))

# Number of nodes
nb.cores <- 19

# Choose methods to compare
imp.methods <- c("mice", "mean", "mia.grf.ate", "mean.grf.ate", "saem")#, "udell")# c("mean", "mia.grf.ate", "udell")#

# Number of imputations for multiple imputations
mi.m <- 10

# Number of experiments per setting
N <- 200#100


col.names <- c("ate", "se.ate", "type", "method", "mechanism", "prop.missing", "use.mask", "use.outcome", "cit.cio", "cit", "cio")


cl <- makeCluster(nb.cores)
registerDoParallel(cl)

for (mech in c("MCAR","MNAR","MAR")){#c("MCAR", "MAR", "MNAR")){
  for (prop.missing in c(0.1, 0.3, 0.5, 0.7)){
    writeLines(paste0("Mechanism: ", mech,", p.miss: ", prop.missing))
    results <- data.frame(t(rep(NA, length(col.names))))
    colnames(results) <- col.names
    results <- results[-1,]
    results <- foreach(seed=1:N, .export=c("expit"),
                       .packages = c("MASS", "caret", "grf",
                                     "missMDA", "missForest", "dbarts",
                                     "mice", "misaem", "softImpute", "norm", "assertthat", "pracma"),
                       .combine = 'rbind') %dopar% {
                         # for (seed in 1:N){
                         #   print(seed)
     cat(seed," ")
     res_tmp <- data.frame(t(rep(NA, length(col.names))))
     colnames(res_tmp) <- col.names
     res_tmp <- res_tmp[-1,]
     for (cio in c(FALSE, TRUE)){                 
       for (use.mask in c(FALSE,TRUE)){
       #writeLines(paste0("cit: ", cit,", cio: ", cio, ", use.mask: ", use.mask))
       
       generateDataForIterInCurrentEnvironment(seed, x[,1:6], z, w,  overlap = "all")
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
       
       generateDataForIterInCurrentEnvironment(seed, x.r[,1:6], z, w,  mask = is.na(x.miss), cio = cio, overlap = "all")
       
       x.miss <- x.miss[,1:6]
       for (imp in imp.methods){
         
         if (imp %in% c("saem", "udell")){
           x.imp.fitted <- prepare_data_ate(df=x.miss[,which(sapply(x.miss, is.numeric))], w=z.r, y=y, 
                                            imputation.method=imp, mi.m = mi.m,
                                            mask = NULL,
                                            use.outcome=F, use.interaction=F)
           x.miss <- x.miss[,which(sapply(x.miss, is.numeric))]
           x.imp <- x.imp.fitted$df.imp
           fitted <- x.imp.fitted$fitted
         } 
         else if (imp != "mice") {
           x.imp.fitted <- prepare_data_ate(df=x.miss, w=z.r, y=y, 
                                            imputation.method=imp, mi.m = mi.m,
                                            mask = NULL,
                                            use.outcome=F, use.interaction=F)
           x.imp <- x.imp.fitted$df.imp
           fitted <- x.imp.fitted$fitted
         }
         
         
         
         if (imp %in% c("mean", "saem", "pc", "udell", "mice")){ 
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
         
         if (imp != "mice"){
           res <- ipw(X=data.frame(x.imp), outcome = y, treat = z.r,
                      ps.method = ps.method, seed = seed,
                      mask = mask,
                      fitted = fitted)
           tmp1 <- cbind(res[1], NA,"ipw1", imp, mech, prop.missing, use.mask, FALSE, paste0("-.",cio), "-", cio)
           colnames(tmp1) <- col.names
           
           tmp2 <- cbind(res[2], NA,"ipw2", imp, mech, prop.missing, use.mask, FALSE, paste0("-.",cio), "-", cio)
           colnames(tmp2) <- col.names
           
           res <- data.frame(dr = c(NA,NA))
           
           try(res <- dr(X=data.frame(x.imp), outcome = y, treat = z.r,
                         ps.method= ps.method,  
                         target= "all", 
                         seed = seed,
                         out.method = out.method,
                         mask = mask,
                         fitted = fitted))
           tmp3 <- cbind(res$dr, res$se,"dr", imp, mech, prop.missing, use.mask, FALSE, paste0("-.",cio), "-", cio)
           colnames(tmp3) <- col.names
           
           res_tmp <- rbind(res_tmp, tmp1, tmp2, tmp3)
           
         } else {
           for (use.outcome in c(F,T)){
             x.imp.fitted <- prepare_data_ate(df=x.miss, w=z.r, y=y, 
                                              imputation.method=imp, mi.m = mi.m,
                                              mask = NULL, use.outcome = use.outcome,
                                              use.interaction=F)
             x.imp <- x.imp.fitted$df.imp
             fitted <- x.imp.fitted$fitted
             
             res <- c()
             for (k in 1:mi.m){
               try(res <- rbind(res, ipw(X = data.frame(x.imp[[k]]), outcome = y, treat = z.r, 
                                         ps.method = ps.method, seed = seed,
                                         mask = mask,
                                         fitted = fitted)))
               
             }
             tmp1 <- cbind(mean(res[,1]), NA, "ipw1", imp, mech, prop.missing, use.mask, use.outcome, paste0("-.",cio), "-", cio)
             colnames(tmp1) <- col.names
             
             tmp2 <- cbind(mean(res[,2]), NA, "ipw2", imp, mech, prop.missing, use.mask, use.outcome, paste0("-.",cio), "-", cio)
             colnames(tmp2) <- col.names
             
             res <- c()
             for (k in 1:mi.m){
               try(res <- rbind(res, dr(X=data.frame(x.imp[[k]]), outcome = y, treat = z.r,
                                        ps.method= ps.method,  
                                        target= "all", 
                                        seed = seed,
                                        out.method = out.method,
                                        mask = mask,
                                        fitted = fitted)))
               
             }
             tmp3 <- cbind(mean(res[,1]), mean(res[,2]), "dr", imp, mech, prop.missing, use.mask, use.outcome, paste0("-.",cio), "-", cio)
             colnames(tmp3) <- col.names
             
             res_tmp <- rbind(res_tmp, tmp1, tmp2, tmp3)
           }
           
         }
                                 
                              
         }
       }
       #results <- rbind(results, res_tmp)
     }
     res_tmp
    }
    cat("\n")
    results_compare_dr_miss <- results[which(results$type == "dr"),]
    results_compare_ipw_miss <- results[which(results$type %in% c("ipw1","ipw2")),]
    results_compare_dr_miss$ate <- as.numeric(as.character(results_compare_dr_miss$ate))
    results_compare_dr_miss$se.ate <- as.numeric(as.character(results_compare_dr_miss$se.ate))
    results_compare_ipw_miss$ate <- as.numeric(as.character(results_compare_ipw_miss$ate))
    results_compare_ipw_miss$se.ate <- as.numeric(as.character(results_compare_ipw_miss$se.ate))
    results_compare_dr_miss$prop.missing <- as.numeric(as.character(results_compare_dr_miss$prop.missing))
    results_compare_ipw_miss$prop.missing <- as.numeric(as.character(results_compare_ipw_miss$prop.missing))
    
    save(results_compare_dr_miss,
         results_compare_ipw_miss,
         N,
         prop.missing,
         file = paste0("./RData/", Sys.Date(), "_ihdp_orig_N", N, "_", mech, "_propNA",
                       prop.missing, "_results_compare_miss.RData"))
  }
}


stopCluster(cl)
