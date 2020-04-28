nb_cores <- max(parallel::detectCores(logical=F)-1,1)

###################################################################################
## Define functions for simulations on complete data
###################################################################################
ipw_estimation <- function(N,  n, setting = 1, Sigma = diag(p), link = "strong", trimming_weight = 1, covariates_link = "linear"){
  results <- data.frame("ipw1_true" = rep(NA, N), "ipw2_true" = rep(NA,N),
                        "ipw1_false" = rep(NA, N), "ipw2_false" = rep(NA,N))
  
  if (trimming_weight == 1) {
    trimming <- FALSE
  } else {
    trimming <- TRUE
  }
  
  cl <- makeCluster(nb_cores)
  registerDoParallel(cl)
  
  t <- Sys.time()
  results <- foreach(i=1:N, 
                     .export=c("generate_sim", "ipw", "expit",
                               "p", "r", "p.udell", "alpha.star.strong", "alpha.star.moderate", "alpha.star.low",
                               "beta.star", "tau", "predict_glm", "predict_grf", "make_V", "design_matrix",
                               "perturbation_gaussian"), 
                     .packages = c("MASS", "caret", "assertthat", "grf"), 
                     .combine = rbind) %dopar% {
                                        
       sample <- generate_sim(n=n, seed = i, link=link, setting = setting, Sigma = Sigma, covariates_link = covariates_link)
       results[i,] <- cbind(ipw(X = data.frame(sample$Z),
                                outcome = sample$y, 
                                treat = sample$treat, 
                                ps.method="glm", 
                                target= "all",
                                seed = i,
                                trimming_weight = trimming_weight),
                            ipw(X = data.frame(sample$X),
                                outcome = sample$y, 
                                treat = sample$treat, 
                                ps.method="glm", 
                                target= "all", 
                                seed = i,
                                trimming_weight = trimming_weight))
  }
  stopCluster(cl)
  results <- data.frame(results)
  results <- cbind(results, 
                   n = rep(n, dim(results)[1]), 
                   setting = rep(setting, dim(results)[1]),
                   link = rep(link, dim(results)[1]),
                   trimming = rep(trimming, dim(results)[1]),
                   covariates_link = rep(covariates_link, dim(results)[1]))
  
  colnames(results) <- c("ipw1_true", "ipw2_true", 
                         "ipw1_false", "ipw2_false", 
                         "n", "setting", "link", "trimming", "covariates_link")
  writeLines("Running time:")
  print(Sys.time()-t)
  return(results) 
}

dr_estimation <- function(N, n, 
                          setting = 1, 
                          Sigma = diag(p), 
                          link = "strong", 
                          trimming_weight = 1, covariates_link = "linear"){
  results <- data.frame("dr_true_true" = rep(NA, N),
                        "dr_false_true" = rep(NA, N),
                        "dr_true_false" = rep(NA, N),
                        "dr_false_false" = rep(NA, N))
  
  if (trimming_weight == 1) {
    trimming <- FALSE
  } else {
    trimming <- TRUE
  }
  
  cl <- makeCluster(nb_cores)
  registerDoParallel(cl)
  
  t <- Sys.time()
  
  results <- foreach(i=1:N, 
                     .export=c("generate_sim", "dr", "expit", "min_inf_proxy", "max_inf_proxy",
                               "p", "r", "p.udell", "alpha.star.strong", "alpha.star.moderate", "alpha.star.low",
                               "beta.star", "tau", "predict_glm", "predict_grf", "make_V", "design_matrix",
                               "perturbation_gaussian", "get_imputeEM", "imputeEllP"), 
                     .packages = c("MASS", "caret", "assertthat", "norm","grf"),
                     .combine = rbind) %dopar% {
     sample <- generate_sim(n=n, seed = i, link=link, setting = setting, Sigma = Sigma, covariates_link = covariates_link)
     results[i, ] <- cbind(dr(X = data.frame(sample$Z),
                              outcome = sample$y, 
                              treat = sample$treat, 
                              ps.method="glm", 
                              out.method = "glm",
                              target= "all",
                              seed = i,
                              trimming_weight = trimming_weight),
                           dr(X = data.frame(sample$Z),
                              X.for.prop = data.frame(sample$X),
                              outcome = sample$y, 
                              treat = sample$treat, 
                              ps.method="glm", 
                              out.method = "glm",
                              target= "all",
                              seed = i,
                              trimming_weight = trimming_weight),
                           dr(X = data.frame(sample$Z),
                              X.for.outcome = data.frame(sample$X),
                              outcome = sample$y, 
                              treat = sample$treat, 
                              ps.method="glm", 
                              out.method = "glm",
                              target= "all", 
                              seed = i,
                              trimming_weight = trimming_weight),
                           dr(X = data.frame(sample$X),
                              outcome = sample$y, 
                              treat = sample$treat, 
                              ps.method="glm", 
                              out.method = "glm",
                              target= "all", 
                              seed = i,
                              trimming_weight = trimming_weight))
  }
  results <- data.frame(results)
  
  results <- cbind(results, 
                   n = rep(n, dim(results)[1]), 
                   setting = rep(setting, dim(results)[1]),
                   link = rep(link, dim(results)[1]),
                   trimming = rep(trimming, dim(results)[1]),
                   covariates_link = rep(covariates_link, dim(results)[1]))
  
  colnames(results) <- c("dr_true_true", 
                         "dr_false_true",
                         "dr_true_false",
                         "dr_false_false", 
                         "n", "setting", "link", "trimming", "covariates_link")
  
  stopCluster(cl)
  writeLines("Running time:")
  print(Sys.time()-t)
  return(results) 
}




###################################################################################
## Define functions for simulations on incomplete data
###################################################################################

prepare_data_ate <- function(df, w, y, imputation.method, mi.m, 
                             mask,
                             use.outcome, use.interaction){ 
  df.imp <- fitted <- NULL                                      
  if (tolower(imputation.method) %in% c("mean")){
    df.imp <- get_imputeMean(data.frame(df))
  } 
  if (tolower(imputation.method) %in% c("mean.grf", "mean.grf.ate")){
    df.imp <- get_imputeMeanNA(data.frame(df))
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
      pca_soft <- recover_pca_gaussian_cv(as.matrix(df), 1:dim(df)[2], nfolds = 3)
      df.imp <- data.frame(pca_soft$Uhat)
    } else {
      var.type <- sapply(data.frame(df), FUN = function(x) {if_else(is.numeric(x), "gaussian", "binomial")})
      df[,var.type=="binomial"] <- sapply(df[,var.type=="binomial"], as.integer)
      lambdas <- cv.mimi(y=data.frame(df), model = "low-rank", var.type, 
                        algo = "bcgd", 
                        maxit = 100, max.rank = NULL, trace.it = F, parallel = F,
                        len = 15)
      glrm_mod <- mimi::mimi(y=data.frame(df), model = "low-rank",
                                var.type = var.type, max.rank = dim(df), algo="bcgd", lambda1=lambdas$lambda)
      #ncomp <- missMDA::estim_ncpFAMD(df, ncp.max = dim(df)[2]-1)
    }
  }
  return(list(df.imp = df.imp, fitted = fitted))
}


ipw_estimation_miss <- function(N, n, p, r, sd = 0.1,
                                prob, setting = "linear2", ps.dependence = "strong", 
                                link = "linear", 
                                class.interaction = FALSE, sigma.structure = "diagonal",
                                trimming_weight = 1, 
                                imputation.method = "mean",
                                mi.m = 1,
                                nb.cores = nb_cores,
                                use.outcome = FALSE,
                                mechanism = "MCAR",
                                cit = FALSE, cio = FALSE,
                                use.interaction = FALSE,
                                use.mask = FALSE,
                                lib_path = NULL,
                                local = TRUE,
                                new.mia = FALSE){
  results <- data.frame("ipw1_true" = rep(NA, N), "ipw2_true" = rep(NA,N))
  
  if (trimming_weight == 1) {
    trimming <- FALSE
  } else {
    trimming <- TRUE
  }
  
  if (tolower(imputation.method) %in% c("mia.grf", "mean.grf","mia.grf.ate", "mean.grf.ate")) {
    ps.method <- "grf"
    out.method <- "grf"
  } else {
    ps.method <- "glm"
    out.method <- "glm"
  }
  
  if (imputation.method != "mice"){
    mi.m <- 1
  }
  
  set.seed(0)
  V <- make_V(r, p)
  
  sample <- generate_sim(n=5000, p = p, r = r, ng = ng, sd = sd,
                           setting = setting,
                           mechanism=mechanism, prop.missing = prob,
                           cit = cit, cio = cio,
                           V = V,
                           class.interaction = FALSE,
                           ps.dependence = ps.dependence,
                           link = link,
                           sigma.structure = sigma.structure, 
                           seed = 0)
  tau <- mean(sample$tau)

  t <- Sys.time()

  if (!is.null(lib_path)){
    .libPaths(lib_path)
    #clusterEvalQ(cl, .libPaths(lib_path))
  }

  if (local){
    cl <- makeCluster(nb.cores)
  } else {
    cl <- parallel::makePSOCKcluster(nb.cores, outfile='')
  }
  
  registerDoParallel(cl)
  
  results <- foreach(i=1:N, 
                     .export=c("generate_sim", "gen_linear", "gen_latentclass", "gen_dlvm", "gen_multisvd",
                               "gen_treat_lat", "gen_treat_deep", "gen_treat_svd",
                               "gen_out_lat", "gen_out_deep", "gen_out_svd",
                               "ipw", "expit", 
                               "alpha.star.strong", "alpha.star.moderate", "alpha.star.low",
                               "beta.star",
                               "tau", "predict_glm",
                               "latent_confounders",
                               "get_imputeMean", "get_imputeMeanNA", "get_imputePC",
                               "get_MICE", "get_MISSFOREST",
                               "get_imputeInf", "predict_grf",
                               "predict_misaem", "cast_types",
                               "make_V", "design_matrix", "perturbation_gaussian",
                               "cross_valid", "compute_folds",
                               "recover_pca_gaussian_cv", "get_imputeEM", "imputeEllP",
                               "miss.saem.v2",
                               "prepare_data_ate"), 
                     .packages = c("MASS", "caret", "grf",
                                   "missMDA", "missForest",
                                   "mice", "misaem", "softImpute", "norm", "assertthat", "pracma"),
                     .combine = rbind) %dopar% {
    sample <- generate_sim(n=n, p = p, r = r, ng = ng, sd = sd,
                           setting = setting,
                           mechanism=mechanism, prop.missing = prob,
                           cit = cit, cio = cio,
                           V = V,
                           class.interaction = FALSE,
                           ps.dependence = ps.dependence,
                           link = link,
                           sigma.structure = sigma.structure, 
                           seed = i)
    
    fitted <- NULL
    X.mask <- NULL
    X.imp <- NULL
    
    
  
    # Add the mask
    if (use.mask) {
      X.mask <- data.frame(is.na(sample$X.incomp))
      missing.cov <- sapply(X.mask, FUN = function(x) length(unique(x)) == 2)
      X.mask <- X.mask[,which(missing.cov)]
    } 
    
    try(tmp <- prepare_data_ate(sample$X.incomp, sample$treat, sample$y, 
                            imputation.method, mi.m, X.mask,
                            use.outcome, use.interaction))
    try(X.imp <- tmp$df.imp)
    try(fitted <- tmp$fitted)
    
    if (new.mia & imputation.method %in% c("mia.grf", "mia.grf.ate")){
      X.imp <- sample$X.incomp
      fitted <- NULL
    }
    
    if (mi.m == 1){
      try(results[i,] <- ipw(X = X.imp,
                         outcome = sample$y, 
                         treat = sample$treat, 
                         ps.method = ps.method, 
                         target = "all",
                         seed = i,
                         trimming_weight = trimming_weight,
                         fitted = fitted,
                         mask = X.mask,
                         use.interaction = use.interaction))
    } else {
      res <- c()
      for (k in 1:mi.m){
        try(res <- rbind(res, ipw(X = X.imp[[k]],
                     outcome = sample$y, 
                     treat = sample$treat, 
                     ps.method = ps.method, 
                     target = "all",
                     seed = i,
                     trimming_weight = trimming_weight,
                     fitted = fitted,
                     mask = X.mask,
                     use.interaction = use.interaction)))
        
      }
      results[i,] <- t(apply(data.frame(res), FUN = mean, MARGIN = 2))
    }
  }
  
  stopCluster(cl)
  results <- data.frame(results)
  
  results <- cbind(results, 
                   n = rep(n, N), 
                   p = rep(p, N), 
                   r = rep(r, N), 
                   tau = rep(tau, N),
                   setting = rep(setting, N),
                   link = rep(link, N),
                   trimming = rep(trimming, N),
                   ps.method = rep(ps.method, N),
                   out.method = rep("-", N),
                   use.outcome = rep(use.outcome, N),
                   imputation.method = rep(imputation.method, N),
                   mi.m = rep(mi.m, N),
                   use.mask = rep(use.mask, N),
                   prop.missing = rep(prob, N),
                   mechanism = rep(mechanism, N),
                   cit = rep(cit, N),
                   cio = rep(cio, N),
                   use.interaction = rep(use.interaction, N))
  
  
  
  colnames(results) <- c("ipw1_true", "ipw2_true", 
                         "n", "p", "r", "tau",
                         "setting", "link", "trimming", 
                         "ps.method", "out.method", "use.outcome",
                         "imputation.method", "mi.m",
                         "use.mask",
                         "prop.missing",
                         "mechanism",
                         "cit", "cio",
                         "use.interaction") 
  
  writeLines("Running time:")
  print(Sys.time()-t)
  return(results) 
}



dr_estimation_miss <- function(N, n, p, r, sd = 0.1,
                               prob, setting = "linear2", ps.dependence = "strong", 
                               link = "linear", 
                               class.interaction = FALSE, sigma.structure = "diagonal",
                               trimming_weight = 1, imputation.method="mean", 
                               mi.m = 1, 
                               nb.cores = nb_cores,
                               use.outcome = FALSE,
                               mechanism = "MCAR",
                               cit = FALSE, cio = FALSE,
                               use.interaction = FALSE,
                               use.mask = FALSE,
                               lib_path = NULL,
                               local = TRUE,
                               new.mia = FALSE){
  results <- data.frame("dr_true_true" = rep(NA, N), "se" = rep(NA,N))
  
  if (trimming_weight == 1) {
    trimming <- FALSE
  } else {
    trimming <- TRUE
  }
  
  if (tolower(imputation.method) %in% c("mia.grf", "mean.grf")) {
    ps.method <- "grf"
    out.method <- "grf"
  } else if (tolower(imputation.method) %in% c("mia.grf.ate", "mean.grf.ate")) {
    ps.method <- "grf.ate"
    out.method <- "grf.ate"
  } else {
    ps.method <- "glm"
    out.method <- "glm"
  }
  
  
  
  set.seed(0)
  V <- make_V(r, p)
  
  sample <- generate_sim(n=5000, p = p, r = r, ng = ng, sd = sd,
                           setting = setting,
                           mechanism=mechanism, prop.missing = prob,
                           cit = cit, cio = cio,
                           V = V,
                           class.interaction = FALSE,
                           ps.dependence = ps.dependence,
                           link = link,
                           sigma.structure = sigma.structure, 
                           seed = 0)
  tau <- mean(sample$tau)

  t <- Sys.time()
  
  
  if (!is.null(lib_path)){
    .libPaths(lib_path)
    #clusterEvalQ(cl, .libPaths(lib_path))
  }

  if (local){
    cl <- makeCluster(nb.cores)
  } else {
    cl <- parallel::makePSOCKcluster(nb.cores, outfile='')
  }

  registerDoParallel(cl)
  
  results <- foreach(i=1:N, 
                     .export=c("generate_sim", "gen_linear", "gen_latentclass", "gen_dlvm", "gen_multisvd",
                               "gen_treat_lat", "gen_treat_deep", "gen_treat_svd",
                               "gen_out_lat", "gen_out_deep", "gen_out_svd",
                               "dr", "expit", 
                               "alpha.star.strong", "alpha.star.moderate", "alpha.star.low",
                               "beta.star",
                               "tau", "predict_glm",
                               "latent_confounders",
                               "get_imputeMean", "get_imputeMeanNA", "get_imputePC",
                               "get_MICE", "get_MISSFOREST",
                               "get_imputeInf", "predict_grf",
                               "predict_misaem", "cast_types",
                               "make_V", "design_matrix", "perturbation_gaussian",
                               "cross_valid", "compute_folds",
                               "recover_pca_gaussian_cv", "get_imputeEM", "imputeEllP",
                               "miss.saem.v2",
                               "prepare_data_ate"), 
                     .packages = c("MASS", "caret", "grf",
                                   "missMDA", "missForest",
                                   "mice", "misaem", "softImpute", "norm", "assertthat", "pracma"),
                     .combine = rbind) %dopar% {
                  
    sample <- generate_sim(n=n, p = p, r = r, ng = ng, sd = sd,
                           setting = setting,
                           mechanism=mechanism, prop.missing = prob,
                           cit = cit, cio = cio,
                           V = V,
                           class.interaction = FALSE,
                           ps.dependence = ps.dependence,
                           link = link,
                           sigma.structure = sigma.structure, 
                           seed = i)

    
    fitted <- NULL
    X.mask <- NULL
    
    if (use.mask){
      X.mask <- data.frame(is.na(sample$X.incomp))
      missing.cov <- sapply(X.mask, FUN = function(x) length(unique(x)) ==2)
      X.mask <- X.mask[,which(missing.cov)] 
    } 
    
    
    try(tmp <- prepare_data_ate(sample$X.incomp, sample$treat, sample$y, 
                            imputation.method, mi.m, X.mask,
                            use.outcome, use.interaction))
    try(X.imp <- tmp$df.imp)
    try(fitted <- tmp$fitted)
    
    if (new.mia & imputation.method %in% c("mia.grf", "mia.grf.ate")){
      X.imp <- sample$X.incomp
      fitted <- NULL
    }
    
  
    if (mi.m == 1){
      try(results[i, ] <- dr(X = X.imp,
                         outcome = sample$y, 
                         treat = sample$treat, 
                         ps.method=ps.method, 
                         target= "all",
                         seed = i,
                         fitted = fitted,
                         trimming_weight = trimming_weight,
                         out.method = out.method,
                         mask = X.mask,
                         use.interaction = use.interaction))
    } else {
      res <- c()
      for (k in 1:mi.m){
        try(res <- rbind(res,dr(X = X.imp[[k]],
                        outcome = sample$y, 
                        treat = sample$treat, 
                        ps.method=ps.method, 
                        target= "all",
                        seed = i,
                        fitted = fitted,
                        trimming_weight = trimming_weight,
                        out.method = out.method,
                        mask = X.mask,
                        use.interaction = use.interaction)))
      }
      results[i,] <- t(apply(data.frame(res), FUN = mean, MARGIN = 2))
    }
  }
  stopCluster(cl)
  
  results <- data.frame(results)
  
  results <- cbind(results, 
                   n = rep(n, N),
                   p = rep(p, N),
                   r = rep(r, N),
                   tau = rep(tau, N),
                   setting = rep(setting, N),
                   link = rep(link, N),
                   trimming = rep(trimming, N),
                   ps.method = rep(ps.method, N),
                   out.method = rep(out.method, N),
                   use.outcome = rep(use.outcome, N),
                   imputation.method = rep(imputation.method, N),
                   mi.m = rep(mi.m, N),
                   use.mask = rep(use.mask, N),
                   prop.missing = rep(prob, N),
                   mechanism = rep(mechanism, N),
                   cit = rep(cit, N),
                   cio = rep(cio, N),
                   use.interaction = rep(use.interaction, N))
  
  colnames(results) <- c("dr_true_true", "se",
                         "n", "p", "r", "tau",
                         "setting", "link", "trimming", 
                         "ps.method", "out.method", "use.outcome",
                         "imputation.method", "mi.m",
                         "use.mask",
                         "prop.missing",
                         "mechanism",
                         "cit", "cio",
                         "use.interaction") 
  
  
  writeLines("Running time:")
  print(Sys.time()-t)
  return(results) 
}



ate_estimation_miss <- function(N, n, p, r, 
                                prob, setting = "linear2", ps.dependence = "strong", 
                                link = "linear", 
                                class.interaction = FALSE, sigma.structure = "diagonal",
                                trimming_weight = 1, imputation.methods = c("mean", "saem"), 
                                mi.m = 1,
                                nb.cores = nb_cores,
                                use.outcome = FALSE, use.mask = FALSE,
                                mechanism = "MCAR",
                                cit = FALSE, cio = FALSE,
                                use.interaction = FALSE,
                                lib_path = NULL,
                                local = TRUE,
                                new.mia = FALSE){
  results_dr_miss <- c()
  results_ipw_miss <- c()
  for (imp in imputation.methods){
    print(imp)
    if (imp == "mice"){
      m <- mi.m
    } else {
      m <- 1
    }
    
    if (imp == "mice"){
      for (use.y in c(FALSE, TRUE)){
        try(results_dr_miss <- rbind(results_dr_miss,
                                 dr_estimation_miss(N, n, p, r, 
                                                    prob= prob, setting = setting, ps.dependence = ps.dependence, 
                                                    link = link, 
                                                    class.interaction = class.interaction, sigma.structure = sigma.structure,
                                                    trimming_weight = trimming_weight, 
                                                    imputation.method = imp, mi.m = m,
                                                    nb.cores = nb.cores,
                                                    use.outcome = use.y,
                                                    use.mask = use.mask,
                                                    mechanism = mechanism,
                                                    cit = cit, cio = cio, 
                                                    use.interaction = use.interaction,
                                                    lib_path = lib_path,
                                                    local = local,
                                                    new.mia = new.mia)))
        
        try(results_ipw_miss <- rbind(results_ipw_miss,
                                  ipw_estimation_miss(N, n, p, r, 
                                                      prob= prob, setting = setting, ps.dependence = ps.dependence, 
                                                      link = link, 
                                                      class.interaction = class.interaction, sigma.structure = sigma.structure,
                                                      trimming_weight = trimming_weight, 
                                                      imputation.method = imp, mi.m = m,
                                                      nb.cores = nb.cores,
                                                      use.outcome = use.y,
                                                      use.mask = use.mask,
                                                      mechanism = mechanism, 
                                                      cit = cit, cio = cio, 
                                                      use.interaction = use.interaction,
                                                      lib_path = lib_path,
                                                      local = local,
                                                      new.mia = new.mia)))
      }
    }
    
    else {
      try(results_dr_miss <- rbind(results_dr_miss,
                               dr_estimation_miss(N, n, p, r, 
                                                  prob= prob, setting = setting, ps.dependence = ps.dependence, 
                                                  link = link, 
                                                  class.interaction = class.interaction, sigma.structure = sigma.structure,
                                                  trimming_weight = trimming_weight, 
                                                  imputation.method = imp, mi.m = m,
                                                  nb.cores = nb.cores,
                                                  use.outcome = FALSE,
                                                  use.mask = use.mask,
                                                  mechanism = mechanism,
                                                  cit = cit, cio = cio, 
                                                  use.interaction = use.interaction,
                                                  lib_path = lib_path,
                                                  local = local,
                                                  new.mia = new.mia)))
      
      try(results_ipw_miss <- rbind(results_ipw_miss,
                                ipw_estimation_miss(N, n, p, r, 
                                                    prob= prob, setting = setting, ps.dependence = ps.dependence, 
                                                    link = link, 
                                                    class.interaction = class.interaction, sigma.structure = sigma.structure,
                                                    trimming_weight = trimming_weight, 
                                                    imputation.method = imp, mi.m = m,
                                                    nb.cores = nb.cores,
                                                    use.outcome = FALSE,
                                                    use.mask = use.mask,
                                                    mechanism = mechanism, 
                                                    cit = cit, cio = cio, 
                                                    use.interaction = use.interaction,
                                                    lib_path = lib_path,
                                                    local = local,
                                                    new.mia = new.mia)))
    }
  }
  
  return(list(results_dr_miss = results_dr_miss, results_ipw_miss = results_ipw_miss))
}