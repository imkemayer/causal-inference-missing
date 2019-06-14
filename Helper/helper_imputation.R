# Mean/mode imputation
get_imputeMean <- function(df, args_list = NULL){
  df_imp <- df
  vars_real <- colnames(df)[sapply(df, is.numeric) & !sapply(df, is.integer)]
  vars_int <- colnames(df)[sapply(df, is.integer)]
  vars_factor <- colnames(df)[!sapply(df, is.numeric)]
  
  if (length(vars_real) >0){
    mean_real <- sapply(df[,vars_real], mean, na.rm=T)
    df_imp[ , vars_real] <- sapply(1:length(vars_real), 
                                   function(x) ifelse(is.na(df_imp[,vars_real[x]]), mean_real[x], df_imp[,vars_real[x]]))
  }
  
  if (length(vars_int) > 0 ){
    mode_int <- as.vector(apply(data.frame(df[, vars_int]), MARGIN=2,
                                FUN = function(x) as.integer(as.character(as.data.frame(table(x))[which.max(table(x)),1]))))
    df_imp[ , vars_int] <- sapply(1:length(vars_int), 
                                  function(x) ifelse(is.na(df_imp[,vars_int[x]]), mode_int[x], df_imp[,vars_int[x]]))
  }
  
  if (length(vars_factor) > 0 ){
    fact.levels <- sapply(vars_factor, FUN=function(x) levels(df[,x]))
    mode_factor <- as.vector(apply(data.frame(df[, vars_factor]), MARGIN=2, 
                                   FUN = function(x) as.data.frame(table(x))[which.max(table(x)),1]))
    df_imp[ , vars_factor] <- sapply(1:length(vars_factor), 
                                     function(x) ifelse(is.na(df_imp[,vars_factor[x]]), mode_factor[x], as.character(df_imp[,vars_factor[x]])))
    for (i in 1:length(vars_factor)){
      df_imp[,vars_factor[i]] <- as.factor(df_imp[,vars_factor[i]])
      levels(df_imp[,vars_factor[i]]) <- levels(df[,vars_factor[i]])
    }
  }
  return(df_imp)
}

# Mean imputation with noise
get_imputeNoisyMean <- function(df){
  imputedData <- df
  vars_real <- colnames(df)[sapply(df, is.numeric) & !sapply(df, is.integer)]
  vars_int <- colnames(df)[sapply(df, is.integer)]
  vars_factor <- colnames(df)[!sapply(df, is.numeric)]
  
  if (length(vars_real) >0){
    mean_real <- sapply(df[,vars_real], mean, na.rm=TRUE)
    sd_real <- sapply(df[,vars_real], sd, na.rm=TRUE)
    for (i in 1:length(vars_real)){
      rows <- which(is.na(imputedData[ , vars_real[i] ]))
      imputedData[ rows , vars_real[i] ] <- rnorm(length(rows), mean_real[i], sd_real[i])
    }
  }
  
  if (length(vars_int) > 0 ){
    mode_int <- as.vector(sapply(df[, vars_int], 
                                 FUN = function(x) as.integer(as.character(as.data.frame(table(x))[which.max(table(x)),1]))))
    for (i in 1:length(vars_int)){
      imputedData[ which(is.na(imputedData[ , vars_int[i] ])), vars_int[i] ] <- mode_int[i]
    }
  }
  
  if (length(vars_factor) > 0 ){
    mode_factor <- as.vector(sapply(df[, vars_factor], 
                                    FUN = function(x) as.character(as.data.frame(table(x))[which.max(table(x)),1])))
    for (i in 1:length(vars_factor)){
      imputedData[ which(is.na(imputedData[ , vars_factor[i] ])), vars_factor[i] ] <- mode_factor[i]
    }
  }
  return(imputedData)
}

# Imputation using PCA/MCA/FAMD
get_imputePC <- function(df, seed, 
                         ncp=3, Method = "Regularized", scale = TRUE, threshold = 1e-06) {
  data.num <- sapply(df, is.numeric)
  if (sum(data.num) == dim(df)[2]){
    imputed <- imputePCA(df, ncp = ncp, seed = seed)
  } else if (sum(data.num) == 0) {
    imputed <- imputeMCA(df, ncp = ncp, seed = seed)
  } else {
    imputed <- imputeFAMD(df, ncp = ncp, seed = seed)
  }
  imputedData <- data.frame(imputed$completeObs)
  
  
  data.num <- colnames(imputedData)[which(sapply(df, FUN=function(x) (is.numeric(x) & !(is.integer(x)))))]
  
  for (j in colnames(imputedData)){
    imputedData[,j] <- cast_types(j, imputedData, data.num)
  }
  return(imputedData)
}

# Recode values of imputed categorical variables and recast some numericals variables into integers 
cast_types = function(i,df, data.num){
  if (is.factor(df[,i])){
    df[,i] = plyr::mapvalues(df[,i], from = levels(df[,i]), to = gsub(paste(i,"_",sep=''), "", levels(df[,i])))
  } else {
    if(i %in% data.num){
      df[,i] <- round(df[,i],digits=1)
    }  else{
      df[,i] <- as.integer(round(df[,i],digits=0))
    }
  }
  return(df[,i])
}

# Imputation using MICE
get_MICE <- function(df, seed=0, m=1, maxit=5, idx = NULL) {
  if (is.null(idx)){
    idx <- 1:dim(df)[2]
  }
  
  imputed <- mice(df, m=m, maxit=maxit, seed=seed, nnet.MaxNWts=5500, printFlag = FALSE)
  if (m > 1){
    imputedData <- list()
    for (k in 1:m){
      imputedData[[k]] <- mice::complete(imputed, k)
      imputedData[[k]] <- imputedData[[k]][,idx]
    }
  } else {
    imputedData <- mice::complete(imputed)
    imputedData <- imputedData[,idx]
  }
  return(imputedData)
}

# Imputation using MissForest
get_MISSFOREST <- function(df, maxit=2) {
  imputed <- missForest(df,  maxiter=maxit, ntree=100, variablewise=TRUE, replace=TRUE)
  imputedData <- imputed$ximp
  return(imputedData)
}

get_imputeInf <- function(df, max_inf_proxy = 1e300, min_inf_proxy = -1e300) {
  if (sum(is.na(df))==0){
    return(df)
  } else {
    na.cols <- colnames(df)[sapply(df, function(x) sum(is.na(x)) > 0)]
    
    # Size before duplication of variables with missing values
    p1 <- dim(df)[2]
    
    k <- 1
    for (j in na.cols){
      colnames.old <- colnames(df)
      if (sum(is.na(df[,j]))>0){
        df  <- cbind(df, replicate(1, df[,j]))
        colnames(df) <- c(colnames.old, paste0("X",p1+k))
        k <- k+1
      }
    }
    
    # Remember which variables are numeric
    vars.num <- colnames(df)[sapply(df, is.numeric)]
    
    # Size after duplication of variables with missing values
    p2 <- dim(df)[2]
    
    
    
    replace_NA_minusInf <- function(x){
      if (is.na(x)){
        x <- min_inf_proxy
      }
      return(x)
    }
    
    replace_NA_plusInf <- function(x){
      if (is.na(x)){
        x <- max_inf_proxy
      }
      return(x)
    }
    df[, na.cols] <- data.frame(apply(df[,na.cols], replace_NA_minusInf, MARGIN = c(1,2)))
    
    df[, (p1+1):p2] <- as.data.frame(apply(df[,(p1+1):p2], replace_NA_plusInf, MARGIN = c(1,2)))
    
    
    for (j in vars.num) {
      df[,j] <- as.numeric(df[,j])
    }
    return(df)
  }
}


get_imputeEM <- function(X){
  s <- prelim.norm(X)
  thetahat <- em.norm(s, showits= FALSE, criterion = sqrt(.Machine$double.eps))
  params <- getparam.norm(s, thetahat)
  X.prep <- t(t(X) - params$mu)
  Inv.Sigma.tmp <- solve(params$sigma)
  miss.rowi = which(rowSums(is.na(X.prep)) > 0.5)
  X.new <- X.prep
  X.new <- X.new[miss.rowi,,drop=FALSE]
  X.prep[miss.rowi,] <- t(apply(X.new, 1, imputeEllP, Inv.Sigma.tmp))
  Ximp= (t(t(X.prep) + params$mu))
  return(list(Ximp = Ximp, mu =  params$mu, sigma = params$sigma))
}

imputeEllP <- function(point, Sigma.inv){
  point.new <- point
  d <- length(point)
  index <- which(is.na(point))
  if (length(index) == 1){
    point.new[index] <- -point.new[-index] %*% Sigma.inv[index,-index] /
      Sigma.inv[index,index]
  }else{
    index <- which(is.na(point))
    A <- Sigma.inv[index,index]
    b <- -Sigma.inv[index,(1:d)[-index], drop = FALSE] %*% point[-index]
    point.new[index] <- solve(A) %*% b
  }
  return (point.new)
}
