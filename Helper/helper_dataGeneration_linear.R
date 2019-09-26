###################################################################################
## Define variables and functions for data generation
###################################################################################

alpha.star.strong <- c(0.6, -0.6, 0.6)
alpha.star.moderate <- c(0.3, -0.3, 0.3)
alpha.star.low <- c(0.1, -0.1, 0.1)
alpha.star.strong.s3 <- c(0.6, -0.6, 0.6, -0.6, 0.6, -0.6, -0.6, -0.6, 0.6, 0.6)
alpha.star.moderate.s3 <- c(0.3, -0.3, 0.3, -0.3, 0.3, -0.3, -0.3, -0.3, 0.3, 0.3)
alpha.star.low.s3 <- c(0.1, -0.1, 0.1, -0.1, 0.1, -0.1, -0.1, -0.1, 0.1, 0.1)

gamma.1 <- c(0.1, 0.2, 0.1, -0.1, -0.1)
gamma.2 <- c(1, -0.1, -0.1, -0.5, 0.5)

tau <- 1
beta.star <- c(-1, 1, -1, 2)
beta.star.s3 <- c(-1, -1, 1, -1, 2, 2, 2, 1, -1, -2, 1)


rho <- 0.6

expit <- function(x){
  return(1/(1+exp(-x)))
}

make_V <- function(r, p){
  matrix(rnorm(r*p), nrow = p, ncol = r)
}

design_matrix <- function(V, n, r, p){
  # this function builds the design list with U, V, X components 
  #   U: n * r matrix with entries from standard gaussian distribution
  #   V: n * p matrix given in the input 
  #   X: X = UV^T
  
  design = vector("list")
  design$U = matrix(rnorm(n*r), nrow = n, ncol = r)
  design$V = V
  design$X = design$U %*% t(design$V)
  
  assert_that(are_equal(dim(design$U), c(n, r)))
  assert_that(are_equal(dim(design$V), c(p, r)))
  assert_that(are_equal(dim(design$X), c(n, p)))
  
  design
}

perturbation_gaussian <- function(design, noise_sd = 5){
  # add Gaussian noise to the UV^T matrix, which creates Gaussian noisy proxies 
  n = nrow(design$X)
  p = ncol(design$X)
  design$X + matrix(rnorm(n*p, 0, noise_sd), nrow = n, ncol = p)
}

#Data matrix and parameter matrix simulations 
latent_confounders <- function(V, n=100, p=10,r=3,sig=0.1){ # sig=0.25
  
  design <- design_matrix(V, n, r, p)
  X.noisy <- design$X
  if (sig>0){
    X.noisy <- perturbation_gaussian(design, noise_sd = sig)
  }
  return(list(Z = design$U, X = design$X, Xnoisy = X.noisy))
}


gen_linear <- function(n, p=10, r = 3, setting = "linear1",
                       seed = 0, ps.dependence = "strong", sd = 1, 
                       mechanism=FALSE, prop.missing = 0,
                       cit = FALSE, cio = FALSE, 
                       cit2 = FALSE, cio2 = FALSE,
                       ci2_imp = "mice",
                       link = "log-lin",
                       V = NULL){
  # setting="linear1" 
  #   - case 1: mechanism="MCAR": p orthogonal covariates with NA in all covariates
  #   - case 2: mechanism="MAR": p covariates with small correlation, NA in half of the covariates (with same Pr of missing in both)
  #
  # setting="linear2" 
  #   - case 1: mechanism="MCAR": p covariates with strong correlation, NA in all covariates
  #   - case 2: mechanism="MAR": p covariates with strong correlation, NA in half of the covariates (with same Pr of missing in both)
  #
  # setting="linear3" 
  #   - case 1: mechanism="MCAR": 3 latent confounders, p covariates with NA in all covariates
  #   - case 2: mechanism="MAR": 3 latent confounders, p covariates with NA in half of the covariates with NA depending on p/2 remaining covariates
  #
  # setting="linear4" 
  #   - case 1: mechanism="MCAR": k latent confounders, p covariates with NA in all covariates
  #   - case 2: mechanism="MAR": k latent confounders, p covariates with NA in half of the covariates with NA depending on p/2 remaining covariates
  #
  # cit/cio: one treatment and outcome model per observed response pattern
  # cit2/cio2: missing values in confounders are imputed (multiple imputation with mice and take elementwise average over imputations)
  #            and treatment and outcome are defined using the imputed confounders 
  #            (only applicable for linear1, linear2)
  
  if (setting=="linear1") setting <- 1
  if (setting=="linear2") setting <- 2
  if (setting=="linear3") setting <- 3
  if (setting=="linear4") setting <- 4
  
  set.seed(seed)
  # Generate coefficients of logistic regression for generating MAR
  gamma.mar.z <- (10/p)*cbind(0.5*runif(n=floor(p/2))-0.25, 0.25*runif(n=floor(p/2))-0.125)
  gamma.mar.x <- (20/p)*cbind(runif(n=floor(p/2))-0.5,0.5*runif(n=floor(p/2))-0.25)
  
  if (setting == 1){
    rho <- 0.3
    Sigma <- diag(p) + rho*upper.tri(diag(p)) + rho*lower.tri(diag(p))
  } else if (setting == 2) {
    rho <- 0.6
    Sigma <- diag(p) + rho*upper.tri(diag(p)) + rho*lower.tri(diag(p))
  }
  
  Z <- NULL
  
  # Covariates
  if (setting %in% c(1,2)){
    
      X <- mvrnorm(n=n, mu=rep(1, p), Sigma=Sigma)
      
      # Modified covariates
      X.proxy <- X
      for (j in 1:p){
        if (mod(j,10)==1) X.proxy[,j] <- exp(0.5*X[,j])
        if (mod(j,10)==2) X.proxy[,j] <- X[,j]/(1+exp(X[,j-1])) + 10
        if (mod(j,10)==3) X.proxy[,j] <- X[,j] <- (X[,j-2]*X[,j]/25 + 0.6)^3
        if (mod(j,10)==4) X.proxy[,j] <- 0.5*exp(X[,j])
        if (mod(j,10)==5) X.proxy[,j] <- X[,j]/(1+exp(X[,j-4])) - 5
        if (mod(j,10)==6) X.proxy[,j] <- 0.015*(X[,j]*X[,j-1] + 0.6)^3
        if (mod(j,10)==7) X.proxy[,j] <- 0.2*sin(X[,j]^2)*cos(X[,j-1])
        if (mod(j,10)==8) X.proxy[,j] <- X[,j]*X[,j-6]/(1+exp(X[,j-3]))
        if (mod(j,10)==9) X.proxy[,j] <- 0.015*(X[,j]*X[,j-6])
        if (mod(j,10)==0) X.proxy[,j] <-  X[,j]/(1+exp(X[,j-9])) + 1
      }
      
      V <- NULL
    
  } else if (setting == 3){
    dat <- latent_confounders(V, n=n, p=p, r=r)
    X <- dat$Z
    X.proxy <- dat$Xnoisy
  }
  

  # Missing covariates
  if (setting == 3){
    X.incomp <- X.proxy
  } else {
    X.incomp <- X
  }

  if (!is.null(mechanism)){
    if (setting %in% c(1,2)){
      if (mechanism == "MCAR"){
        
        idx_NA <- matrix(runif(n*(dim(X)[2])), nrow=n, ncol=dim(X)[2]) <= prop.missing
        idx_NA[rowSums(idx_NA)==ncol(X), sample(ncol(X),1)] = FALSE # avoid having empty observations
        X.incomp[idx_NA] <- NA
        
      } else if (mechanism == "MAR"){
        idx_NA <- c()
        for (j in 1:ceiling(p/2)) {
          m <- sapply(expit(1*(setting==1) + 0.09*(setting==2) + X[,(ceiling(p/2)+1):p]%*%gamma.mar.x[,setting]), FUN= function(pr) rbinom(n=1, size=1, prob = pr))
          X.incomp[which(m==1), j] <- NA
          idx_NA <- cbind(idx_NA, m==1)
        }
        idx_NA <- cbind(idx_NA, matrix(FALSE, nrow = n, ncol = floor(p/2)))
        idx_NA <- matrix(idx_NA, nrow=n, ncol = p)
        
      } else if (mechanism == "MNAR"){
        idx_NA <- c()
        for (j in 1:ceiling(p/2)) {
          if (prop.missing > 0.4 & prop.missing < 0.6){
            m <- (mod(j,2)==0)*(X[,j]<median(X[,j])) + (mod(j,2)==1)*(X[,j]>median(X[,j]))
          } else {
            m <- (mod(j,2)==0)*(X[,j]<quantile(X[,j], prop.missing)) + (mod(j,2)==1)*(X[,j]>quantile(X[,j], 1-prop.missing))
          }
          X.incomp[which(m==1), j] <- NA
          idx_NA <- cbind(idx_NA, m==1)
        }
        idx_NA <- cbind(idx_NA, matrix(FALSE, nrow = n, ncol = floor(p/2)))
        idx_NA <- matrix(idx_NA, nrow=n, ncol = p)
      }
    } else if (setting == 3) {
      if (mechanism == "MCAR"){
        idx_NA = matrix(runif(n*dim(X)[2]), nrow=n, ncol=dim(X)[2]) <= prop.missing
        idx_NA[rowSums(idx_NA)==ncol(X), sample(ncol(X),1)] = FALSE # avoid having empty observations
        X.incomp[idx_NA] <- NA

      } 
      if (mechanism == "MAR"){
        idx_NA <- c()
        for (j in 1:ceiling(p/2)) {
          m <- sapply(expit(X[,(ceiling(p/2)+1):p]%*%gamma.mar.x[,1]), FUN= function(p) rbinom(n=1, size=1, prob = p))
          X.incomp[which(m==1), j] <- NA
          idx_NA <- cbind(idx_NA, m==1)
        }
        idx_NA <- cbind(idx_NA, matrix(FALSE, nrow = n, ncol = floor(p/2)))
        idx_NA <- matrix(idx_NA, nrow=n, ncol = p)

      } 
      if (mechanism == "MNAR"){
        idx_NA <- c()
        for (j in 1:ceiling(p/2)) {
          if (prop.missing > 0.4 & prop.missing < 0.6){
            m <- (mod(j,2)==0)*(X.incomp[,j]<median(X.incomp[,j])) + (mod(j,2)==1)*(X.incomp[,j]>median(X.incomp[,j]))
          } else {
            m <- (mod(j,2)==0)*(X.incomp[,j]<quantile(X.incomp[,j], prop.missing)) + (mod(j,2)==1)*(X.incomp[,j]>quantile(X.incomp[,j], 1-prop.missing))
          }
          X.incomp[which(m==1), j] <- NA
          idx_NA <- cbind(idx_NA, m==1)
        }
        idx_NA <- cbind(idx_NA, matrix(FALSE, nrow = n, ncol = floor(p/2)))
        idx_NA <- matrix(idx_NA, nrow=n, ncol = p)
      }
    }
  }
  
  
  X.tmp <- X
  
  # Propensities 
  if (cit & setting %in% c(1,2)) {
    X.tmp[idx_NA] <- 0
  } 
  
  if (cit2 & setting %in% c(1,2)) {
    if (ci2_imp == "mice"){
      X.imp <- get_MICE(X.incomp, seed=0, m=5, maxit=5)
    } else {
      X.imp <- missMDA::MIPCA(X.incomp, ncp = d, nboot = 5)$res.MI
    }
    tmp <-Reduce("+", X.imp) / length(X.imp) # take elementwise average over imputations
    X.tmp[idx_NA] <- tmp[idx_NA]
  } 

  if (ps.dependence=="strong"){
      alpha <- array(alpha.star.strong, dim = dim(X.tmp)[2])
    } else if (ps.dependence == "moderate") {
      alpha <- array(alpha.star.moderate, dim = dim(X.tmp)[2])
    } else {
      alpha <- array(alpha.star.low, dim = dim(X.tmp)[2])
    }
  
  
  offsets <- seq(-100, 100, length.out = 50)
  balanced <- FALSE
  i <- 1
  best_idx <- 1
  min_diff <- 1
  if (link == "nonlinear3"){
    for (j in 1:dim(X.tmp)[2]){
      X.tmp[,j] <- (mod(j,5)==1)*((X.tmp[,j]<quantile(X.tmp[,j],0.7)) + (X.tmp[,j]> quantile(X.tmp[,j],0.2))) +
                   (mod(j,5)==2)*(1/(0.001+exp(X.tmp[,j]*X.tmp[,1]))) +
                   (mod(j,5)==3)*(-(X.tmp[,j])*(X.tmp[,2]>0)) +
                   (mod(j,5)==4)*(-2.5*sqrt(abs(X.tmp[,j]))) +
                   (mod(j,5)==0)*(X.tmp[,3]*X.tmp[,j])
    } 
  }
  while ((i <= length(offsets)) & !(balanced) ){
    if (link == "log-lin"){
      prop_scores <- apply(data.frame(X.tmp), MARGIN=1, 
                          FUN = function(z) expit(offsets[i]+2*(setting==3)+ z%*%alpha))
    } else if (link == "nonlinear3"){
      prop_scores <- apply(data.frame(X.tmp), MARGIN=1, 
                          FUN = function(z) expit(offsets[i]+2*(setting==3)+ z%*%alpha))
    }
    else {
      prop_scores <- rep(0, dim(X.tmp)[1])
      for (j in 1:dim(X.tmp)[2]){
        prop_scores <- prop_scores + (mod(j,5)==1)*(X.tmp[,j]*0.03*cos(5*X.tmp[,j]) - 0.2*X.tmp[,j]^2) +
                                     (mod(j,5)==2)*(-0.5*sin(X.tmp[,j])) +
                                     (mod(j,5)==3)*(0.78*(X.tmp[,j]>0)) +
                                     (mod(j,5)==4)*(-2.5*sqrt(abs(X.tmp[,j]))) +
                                     (mod(j,5)==0)*(X.tmp[,j-2]*X.tmp[,j])
      }
      prop_scores <- expit(offsets[i] + prop_scores)
    }
      
    # Treatment assignment
    treat <- sapply(prop_scores, FUN= function(p) rbinom(n=1, size=1, prob = p)) # ifelse(runif(n, 0, 1) <= c(prob_scores), 1, 0)
    
    # test whether there are at least 30% of the observations in each group
    balanced <- (sum(treat==1)/length(treat) > 0.3 & sum(treat==1)/length(treat)<0.7)
    
    diff <- abs(sum(treat==1)/length(treat) - sum(treat==0)/length(treat))
    if (diff < min_diff){
      best_idx <- i
      min_diff <- diff
    }
    i <- i+1
  }
  
  if (i > length(offsets)){
    prop_scores <- apply(data.frame(X.tmp), MARGIN=1, 
                         FUN = function(z) expit(offsets[best_idx]+2*(setting==3)+ z%*%alpha))
    # Treatment assignment
    treat <- sapply(prop_scores, FUN= function(p) rbinom(n=1, size=1, prob = p))
  }
  
  
  # Outcome
  epsilons <- rnorm(n, sd = sd)
  y <- rep(0,n)
  
  X.tmp <- X
  
  if (cio & setting %in% c(1,2)) {
    X.tmp[idx_NA] <- 0
  }
  if (cio2 & setting %in% c(1,2)) {
    if (ci2_imp == "mice"){
      X.imp <- get_MICE(X.incomp, seed=0, m=5, maxit=5)
    } else {
      X.imp <- missMDA::MIPCA(X.incomp, ncp = d, nboot = 5)$res.MI
    }
    tmp <-Reduce("+", X.imp) / length(X.imp) # take elementwise average over imputations
    X.tmp[idx_NA] <- tmp[idx_NA]
  }  
  if (link == "nonlinear3") {
    for (j in 1:dim(X.tmp)[2]){
      X.tmp[,j] <- (mod(j,5)==1)*((X.tmp[,j]<quantile(X.tmp[,j],0.7)) + (X.tmp[,j]> quantile(X.tmp[,j],0.2))) +
        (mod(j,5)==2)*(1/(0.001+exp(X.tmp[,j]*X.tmp[,1]))) +
        (mod(j,5)==3)*(-(X.tmp[,j])*(X.tmp[,2]>0)) +
        (mod(j,5)==4)*(-2.5*sqrt(abs(X.tmp[,j]))) +
        (mod(j,5)==0)*(X.tmp[,3]*X.tmp[,j])
    }
  }
   


  
  beta <- array(beta.star[2:length(beta.star)], dim = dim(X.tmp)[2])

  
  if (link %in% c("log-lin", "nonlinear3")){
    y[which(treat==1)] <- X.tmp[which(treat==1),]%*%beta + beta.star[1] + tau + epsilons[which(treat==1)]
    y[which(!(treat==1))] <- X.tmp[which(!(treat==1)),]%*%beta + beta.star[1] + epsilons[which(!(treat==1))]
  }
  else {
    y[which(treat==1)] <- exp(X.tmp[which(treat==1),]%*%beta + beta.star[1]) + tau + epsilons[which(treat==1)]
    y[which(!(treat==1))] <- exp(X.tmp[which(!(treat==1)),]%*%beta + beta.star[1]) + epsilons[which(!(treat==1))]
  }
  if (setting == 3){
    X <- X.proxy
    Z <- dat$Z
  }
  
  sample <- list("X" = X,
                 "ps" = prop_scores,
                 "treat" = treat,
                 "y" = y,
                 "X.incomp" = X.incomp,
                 "tau" = rep(tau, n),
                 "class" = Z,
                 "X.proxy" = X.proxy)
  return(sample)
}

