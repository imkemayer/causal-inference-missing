
# corresponding to latent class setting
gen_treat_lat <- function(x, class.interaction = FALSE, link = "nonlinear"){
  if (link == "nonlinear"){
    if (class.interaction){
      class <- x[length(x)]
      x <- x[1:length(x)-1]
      if (length(x) ==3){
        f.x <- as.numeric(-class + 0.03*cos(5*x[1]) - 0.5*sin(x[2]+class) + x[3]*x[2] - 0.2*x[1]^2 + 0.8*exp(-class))
      } else {
        f.x <- as.numeric(-((class-1)^2) + 0.03*cos(5*x[1]) - 0.5*sin(x[2]+class) + 
                     x[3]*x[5] - 0.2*x[1]^2 + 0.8*exp(-class) +
                     0.78*(x[4]>0) + 0.5*class*sqrt(abs(x[10])))
      }
    } else {
      if (length(x)==3){
        f.x <- as.numeric(-2 + 0.03*cos(5*x[1]) - 0.5*sin(x[2]) + x[3]*x[2] - 0.2*x[1]^2)
      } else {
        f.x <- as.numeric(-2 + 0.03*cos(5*x[1]) - 0.5*sin(x[2]) + x[3]*x[5] - 0.2*x[1]^2 +
                            0.78*(x[4]>0) + 2.5*sqrt(abs(x[10])))
      }
    }
    ps <- expit(f.x)
    treat <- sapply(ps, FUN= function(pr) rbinom(n=1, size=1, prob = pr))
  } 
  else if (link == "nonlinear2"){
    if (class.interaction){
      class <- x[length(x)]
      x <- x[1:length(x)-1]
      if (length(x) == 3){
        f.x <- as.numeric((class==1)*(-2 + 0.03*cos(5*x[1]) - 0.5*sin(x[2]) + cos(x[3]/2)*cos(2*x[2])) +
                            (class==2)*(1 - x[1]*x[3] + 0.3*x[2]^3 - 1/(abs(x[2]+x[3])+1)) +
                            (class==3)*(-3 + 0.7*(x[1]>0) + 2.5*sqrt(abs(x[3])) - 0.5*(x[1]+x[2]<1) + (x[3]/(abs(x[2])+0.01)<1)))
      } else {
        f.x <- as.numeric((class==1)*(-2 + 0.03*cos(5*x[1]) - 0.5*sin(x[2]) + cos(x[3]/2)*cos(2*x[5]) - 0.1*sin(x[10])) +
                            (class==2)*(1 - x[1]*x[4] + 0.3*x[5]^3 - 1/(abs(x[2]+x[10])+1) + x[3]*x[2]) +
                            (class==3)*(-3 + 0.7*(x[4]>0) + 2.5*sqrt(abs(x[10])) - 0.5*(x[1]+x[2]<1) + (x[3]/(abs(x[5])+0.01)<1)))
      }
    } else {
      if (length(x)==3){
        f.x <- as.numeric(-2 + 0.03*cos(5*x[1]) - 0.5*sin(x[2]) + x[3]*x[2] - 0.2*x[1]^2)
      } else {
        f.x <- as.numeric(-2 + 0.03*cos(5*x[1]) - 0.5*sin(x[2]) + x[3]*x[5] - 0.2*x[1]^2 +
                            0.78*(x[4]>0) + 2.5*sqrt(abs(x[10])))
      }
    }
    ps <- expit(f.x)
    treat <- sapply(ps, FUN= function(pr) rbinom(n=1, size=1, prob = pr))
  } 
  else if (link == "nonlinear3"){
    f.x <- 0
    for (j in 1:length(x)){
      f.x <- f.x + (mod(j,5)==1)*(x[j]*0.03*cos(5*x[j]) - 0.2*x[j]^2) +
                   (mod(j,5)==2)*(-0.5*sin(x[j])) +
                   (mod(j,5)==3)*(0.78*(x[j]>0)) +
                   (mod(j,5)==4)*(-2.5*sqrt(abs(x[j])))
      if (mod(j,5)==0) f.x <- f.x + (x[j-2]*x[j])
    }
    offsets <- c(0.5, -2, seq(-100, 100, length.out = 50))
    balanced <- FALSE
    i <- 1
    best_idx <- 1
    min_diff <- 1
    while ((i <= length(offsets)) & !(balanced) ){
      ps <- expit(offsets[i] + f.x)
      treat <- sapply(ps, FUN= function(pr) rbinom(n=1, size=1, prob = pr))
      
      balanced <- (sum(treat==1)/length(treat) > 0.3 & sum(treat==1)/length(treat)<0.7)
      
      diff <- abs(sum(treat==1)/length(treat) - sum(treat==0)/length(treat))
      if (diff < min_diff){
        best_idx <- i
        min_diff <- diff
      }
      i <- i+1
    }
    if (i > length(offsets)){
      ps <- expit(offsets[best_idx] + f.x)
      treat <- sapply(ps, FUN= function(p) rbinom(n=1, size=1, prob = p)) 
    }
  }
  else if (link=="linear") {
    if (class.interaction){
      class <- as.integer(x[length(x)])
      beta <- array(c(class*0.15, -(class-1.5)*0.5, -0.1*class^2, 0.1, -0.53, 2/(class+1)), dim = length(x))
      x <- x[1:length(x)-1]
      f.x <- as.numeric(sum(x*beta))
    } else {
      beta <- array(c(0.15, -0.5, -0.1, 1, -0.53, 2), dim = length(x))
      f.x <- as.numeric(sum(x*beta))
    }
    
    offsets <- c(0.5, -2, seq(-100, 100, length.out = 50))
    balanced <- FALSE
    i <- 1
    best_idx <- 1
    min_diff <- 1
    while ((i <= length(offsets)) & !(balanced) ){
      ps <- expit(offsets[i] + f.x)
      treat <- sapply(ps, FUN= function(pr) rbinom(n=1, size=1, prob = pr))
      
      balanced <- (sum(treat==1)/length(treat) > 0.3 & sum(treat==1)/length(treat)<0.7)
      
      diff <- abs(sum(treat==1)/length(treat) - sum(treat==0)/length(treat))
      if (diff < min_diff){
        best_idx <- i
        min_diff <- diff
      }
      i <- i+1
    }
    if (i > length(offsets)){
      ps <- expit(offsets[best_idx] + f.x)
      treat <- sapply(ps, FUN= function(p) rbinom(n=1, size=1, prob = p)) 
    }
  }
  return(c(ps = ps, treat = treat))
}


# corresponding to multilevel svd setting
gen_treat_svd <- function(x, class.interaction = FALSE, link = "nonlinear"){
  if (link == "nonlinear"){
    if (class.interaction){
      class <- as.integer(x[length(x)])
      x <- x[1:length(x)-1]
      f.x <- as.numeric(-((class-1)^2) + 0.03*cos(5*x[1]) - 0.5*sin(x[2]+class) + 
                          x[3]*x[5] - 0.2*x[1]^2 + 0.8*exp(-class) +
                           0.78*(x[4]>0) + 2.5*class*sqrt(abs(x[10])))
    } else {
      f.x <- as.numeric(-2 + 0.03*cos(5*x[1]) - 0.5*sin(x[2]) + x[3]*x[5] - 0.2*x[1]^2 +
                          0.78*(x[4]>0) + 2.5*sqrt(abs(x[10])))
    }
  }
  else if (link == "nonlinear3"){
    f.x <- l0
    for (j in 1:length(x)){
      f.x <- f.x + (mod(j,5)==1)*(x[j]*0.03*cos(5*x[j]) - 0.2*x[j]^2) +
        (mod(j,5)==2)*(-0.5*sin(x[j])) +
        (mod(j,5)==3)*(0.78*(x[j]>0)) +
        (mod(j,5)==4)*(-2.5*sqrt(abs(x[j])))
      if (mod(j,5)==0) f.x <- f.x + (x[j-2]*x[j])
    }
  } else if (link == "linear") {
    if (class.interaction){
      class <- as.integer(x[length(x)])
      x <- x[1:length(x)-1]
      beta <- array(c(class*0.15, -(class-1.5)*0.5, -0.1*class^2, 0.1, -0.53, 2/(class+1)), dim = length(x))
      f.x <- as.numeric(sum(x*beta))
    } else {
      beta <- array(c(0.15, -0.5, -0.1, 1, -0.53, 2), dim = length(x))
      f.x <- as.numeric(sum(x*beta))
    }
  }
  
  offsets <- c(0.5, -2, seq(-100, 100, length.out = 50))
  balanced <- FALSE
  i <- 1
  best_idx <- 1
  min_diff <- 1
  while ((i <= length(offsets)) & !(balanced) ){
    ps <- expit(offsets[i] + f.x)
    treat <- sapply(ps, FUN= function(p) rbinom(n=1, size=1, prob = p))
    
    balanced <- (sum(treat==1)/length(treat) > 0.3 & sum(treat==1)/length(treat)<0.7)
    
    diff <- abs(sum(treat==1)/length(treat) - sum(treat==0)/length(treat))
    if (diff < min_diff){
      best_idx <- i
      min_diff <- diff
    }
    i <- i+1
  }
  if (i > length(offsets)){
    ps <- expit(offsets[best_idx] + f.x)
    treat <- sapply(ps, FUN= function(p) rbinom(n=1, size=1, prob = p)) 
  }
    
  return(c(ps = ps, treat = treat))
}



# corresponding to dlvm setting
gen_treat_deep <- function(x, link = "nonlinear"){
  if (link == "nonlinear"){
    f.x <- as.numeric(-1 + 0.3*cos(5*x[1]) - 0.1*sin(x[2]) - 0.1*x[3]*x[5] +
                        0.78*(x[4]>0) + sqrt(abs(x[10])))
  } 
  else if (link == "nonlinear3"){
    f.x <- 0
    for (j in length(x)){
      f.x <- f.x + (mod(j,5)==1)*(x[j]*0.03*cos(5*x[j]) - 0.2*x[j]^2) +
        (mod(j,5)==2)*(-0.5*sin(x[j])) +
        (mod(j,5)==3)*(0.78*(x[j]>0)) +
        (mod(j,5)==4)*(-2.5*sqrt(abs(x[j])))
      if (mod(j,5)==0) f.x <- f.x + (x[j-2]*x[j])
    }
  }
  else {
    beta <- array(c(0.15, -0.5, -0.1, 1, -0.53, 2), dim = length(x))
    f.x <- as.numeric(sum(x*beta))
  }
  
  offsets <- c(0.5, -1, seq(-100, 100, length.out = 50))
  balanced <- FALSE
  i <- 1
  best_idx <- 1
  min_diff <- 1
  while ((i <= length(offsets)) & !(balanced) ){
    ps <- expit(offsets[i] + f.x)
    treat <- sapply(ps, FUN= function(p) rbinom(n=1, size=1, prob = p))
    
    balanced <- (sum(treat==1)/length(treat) > 0.3 & sum(treat==1)/length(treat)<0.7)
    
    diff <- abs(sum(treat==1)/length(treat) - sum(treat==0)/length(treat))
    if (diff < min_diff){
      best_idx <- i
      min_diff <- diff
    }
    i <- i+1
  }
  if (i > length(offsets)){
    ps <- expit(offsets[best_idx] + f.x)
    treat <- sapply(ps, FUN= function(p) rbinom(n=1, size=1, prob = p)) 
  }
  
  return(c(ps = ps, treat = treat))
}


########################################################################################################################


gen_out_lat <- function(x.w.c, class.interaction = FALSE, sd=0.1, link = "nonlinear"){
  if (link == "nonlinear"){
    if (class.interaction){
      x <- x.w.c[1:(length(x.w.c)-2)]
      w <- x.w.c[length(x.w.c)-1]
      class <- x.w.c[length(x.w.c)]
      eps <- rnorm(n = 1, sd = sd, mean = 0)
      if (length(x) == 3) {
        y <- as.numeric(log(class+1) - (w==0)*(sin(0.4*x[1]*class + 0.154*x[2])+x[3]^2*x[1]*0.5) - (w==1)*((0.25*x[2]^2 - 0.15*x[1]-log(class+1))>0)) + eps
      } else if (length(x) ==10){
        y <- as.numeric(log(class+1) - (w==0)*(sin(0.4*x[1]*class + 0.154*x[2]) + x[3]^2*x[1]*0.5  + 0.5*cos(x[5]^2+1) - 0.152*(x[10]+class)^2) -
                          (w==1)*((0.25*x[2]^2 - 0.15*x[1]-log(class+1) + 0.04*x[10]^2+ 0.12*sqrt(abs(x[5])))>0)) + eps
      }
    } else {
      x <- x.w.c[1:length(x.w.c)-1]
      w <- x.w.c[length(x.w.c)]
      eps <- rnorm(n = 1, sd = sd, mean = 0)
      if (length(x) == 3){
        y <- as.numeric(2.445 - (w==0)*(sin(0.4*x[1] + 0.154*x[2])+x[3]^2*x[1]*0.5) - (w==1)*((0.25*x[2]^2 - 0.15*x[1])>0)) + eps
      } else if (length(x) == 10){
        y <- as.numeric(2.445 - (w==0)*(sin(0.4*x[1] + 0.154*x[2])+x[3]^2*x[1]*0.5 + 0.5*cos(x[5]^2+1) - 0.152*x[10]) - 
                          (w==1)*((0.25*x[2]^2 - 0.15*x[1] + 0.04*x[10]^2 + 0.12*sqrt(abs(x[5])))>0)) + eps
      }
    }
  } else if (link == "nonlinear3"){
    if (class.interaction){
      x <- x.w.c[1:(length(x.w.c)-2)]
      w <- x.w.c[length(x.w.c)-1]
      class <- as.integer(x.w.c[length(x.w.c)])
      eps <- rnorm(n = 1, sd = sd, mean = 0)
      beta <- array(c(-class*0.4, 0.154, 0.5/(class+1), -1.95/(class^2)), dim = length(x))
      y <- as.numeric(exp(log(class+1) + sum(x*beta)) + w) + eps
    } else {
      x <- x.w.c[1:length(x.w.c)-1]
      w <- x.w.c[length(x.w.c)]
      eps <- rnorm(n = 1, sd = sd, mean = 0)
      beta <- array(c(-0.2, 0.154, 0.5, -1.95), dim = length(x))
      y <- as.numeric(exp(0.5 + sum(x*beta)) + w) + eps
    }
    
  } else {
    if (class.interaction){
      x <- x.w.c[1:(length(x.w.c)-2)]
      w <- x.w.c[length(x.w.c)-1]
      class <- as.integer(x.w.c[length(x.w.c)])
      eps <- rnorm(n = 1, sd = sd, mean = 0)
      beta <- array(c(-class*0.4, 0.154, 0.5/(class+1), -1.95/(class^2)), dim = length(x))
      y <- as.numeric(log(class+1) + sum(x*beta) + w) + eps
    } else {
      x <- x.w.c[1:length(x.w.c)-1]
      w <- x.w.c[length(x.w.c)]
      eps <- rnorm(n = 1, sd = sd, mean = 0)
      beta <- array(c(-0.2, 0.154, 0.5, -1.95), dim = length(x))
      y <- as.numeric(0.5 + sum(x*beta) + w) + eps
    }
  }
  return(y)
}


gen_out_svd <- function(x.w.c, class.interaction = FALSE, sd=0.1, link = "nonlinear"){
  if (link == "nonlinear"){
    if (class.interaction){
      x <- x.w.c[1:(length(x.w.c)-2)]
      w <- x.w.c[length(x.w.c)-1]
      class <- x.w.c[length(x.w.c)]
      eps <- rnorm(n = 1, sd = sd, mean = 0)
      y <- as.numeric(log(class+1) - (w==0)*(sin(0.4*x[1]*class + 0.154*x[2]) + x[3]^2*x[1]*0.5  + 0.5*cos(x[5]^2+1) - 0.152*(x[10]+class)^2) -
                        (w==1)*((0.25*x[2]^2 - 0.15*x[1]-log(class+1) + 0.04*x[10]^2+ 0.12*sqrt(abs(x[5])))>0)) + eps
    } else {
      x <- x.w.c[1:length(x.w.c)-1]
      w <- x.w.c[length(x.w.c)]
      eps <- rnorm(n = 1, sd = sd, mean = 0)
      y <- as.numeric(2.445 - (w==0)*(sin(0.4*x[1] + 0.154*x[2])+x[3]^2*x[1]*0.5 + 0.5*cos(x[5]^2+1) - 0.152*x[10]) - 
                        (w==1)*((0.25*x[2]^2 - 0.15*x[1] + 0.04*x[10]^2 + 0.12*sqrt(abs(x[5])))>0)) + eps
    }
  } else if (link == "nonlinear3"){
    if (class.interaction){
      x <- x.w.c[1:(length(x.w.c)-2)]
      w <- x.w.c[length(x.w.c)-1]
      class <- as.integer(x.w.c[length(x.w.c)])
      eps <- rnorm(n = 1, sd = sd, mean = 0)
      beta <- array(c(-class*0.4, 0.154, 0.5/(class+1), -1.95/(class^2), 0.1, 0), dim = length(x))
      y <- as.numeric(exp(log(class+1) + sum(x*beta)) + w) + eps
    } else {
      x <- x.w.c[1:length(x.w.c)-1]
      w <- x.w.c[length(x.w.c)]
      eps <- rnorm(n = 1, sd = sd, mean = 0)
      beta <- array(c(-0.2, 0.154, 0.5, -1.95, 0.1, 0), dim = length(x))
      y <- as.numeric(exp(0.5 + sum(x*beta)) + w) + eps
    }
    
  } else {
    if (class.interaction){
      x <- x.w.c[1:(length(x.w.c)-2)]
      w <- x.w.c[length(x.w.c)-1]
      class <- as.integer(x.w.c[length(x.w.c)])
      eps <- rnorm(n = 1, sd = sd, mean = 0)
      beta <- array(c(-class*0.4, 0.154, 0.5/(class+1), -1.95/(class^2), 0.1, 0), dim = length(x))
      y <- as.numeric(log(class+1) + sum(x*beta) + w) + eps
    } else {
      x <- x.w.c[1:length(x.w.c)-1]
      w <- x.w.c[length(x.w.c)]
      eps <- rnorm(n = 1, sd = sd, mean = 0)
      beta <- array(c(-0.2, 0.154, 0.5, -1.95, 0.1, 0), dim = length(x))
      y <- as.numeric(0.5 + sum(x*beta) + w) + eps
    }
  }
  return(y)
}

gen_out_deep <- function(x.w.c, sd=0.1, link = "nonlinear"){
  x <- x.w.c[1:length(x.w.c)-1]
  w <- x.w.c[length(x.w.c)]
  eps <- rnorm(n = 1, sd = sd, mean = 0)
  
  if (link == "nonlinear"){
    y <- as.numeric(2.445 - (w==0)*((sin(0.4*x[1] + 0.154*x[2]) - 1/(x[3]^2*x[1]) + 0.5*cos(x[5]^2+1) - 0.152*x[10])>0) +
      (w==1)*((0.25*x[2]^2 - 0.15*x[1] + 0.04*x[10]^2 + 0.12*sqrt(abs(x[5])))>0)) + eps 
  } else {
    beta <- array(c(-0.2, 0.154, 0.5, -1.95, 0.1, 0), dim = length(x))
    y <- as.numeric(0.5 + sum(x*beta) + w) + eps
  }
  return(y)
}



########################################################################################################################


gen_latentclass <- function(n, p, nb.class=3, mus = NULL, Sigmas = NULL, class.interaction = FALSE,
                            sd=0.1, 
                            seed = 0,
                            mechanism=FALSE, prop.missing = 0, 
                            cit = FALSE, cio = FALSE,
                            link = "nonlinear"){
  set.seed(seed)
  gamma.mar.x <- (20/p)*(runif(n=floor(p/2))-0.5)

  if (p == 3){
    if (is.null(mus)){
      mus <- list(t(rep(0, p)), t(rep(2.5, p)), t(array(c(-3, -2, 1), dim = p)))
    }
    
    if (is.null(Sigmas)){
      Sigmas <- list(diag(p), 
                     diag(p) + 0.3 *upper.tri(diag(p)) + 0.3*lower.tri(diag(p)),
                     matrix(c(2.87, -1.76, -1, -1.76, 6.56, 0.5, -1, 0.5, 3.54), ncol = 3))
    }
  }
  if (p != 3){
    if (is.null(mus)){
      mus <- list(t(rep(0, p)), 
                  t(rep(2.5, p)),
                  t(array(c(-3, -2, 1), dim = p)))
    }
    if (is.null(Sigmas)){
      Sigmas <- list(diag(p), 
                     diag(p) + 0.3 *upper.tri(diag(p)) + 0.3*lower.tri(diag(p)),
                     spcov::GenerateCliquesCovariance(ncliques=3, cliquesize=p/3, -1)$Sigma)
    }
  }
  
  class = sample(1:nb.class, n, replace = TRUE)
  X <- matrix(0, ncol = p, nrow = n)
  for (c in 1:nb.class){
    X <- X + (class == c) * mvrnorm(n, mu = mus[[c]], Sigma = Sigmas[[c]])
  }
  
  X <- data.frame(X)
  
  # MISSING DATA
  X.incomp <- X
  idx_NA <- NULL
  if (mechanism == "MCAR"){
    idx_NA = matrix(runif(n*dim(X)[2]), nrow=n, ncol=dim(X)[2]) <= prop.missing
    idx_NA[rowSums(idx_NA)==ncol(X), sample(ncol(X),1)] = FALSE # avoid having empty observations
    X.incomp[idx_NA] <- NA
  }
  if (mechanism == "MAR"){
    idx_NA <- c()
    for (j in 1:ceiling(p/2)) {
      m <- sapply(expit(2 + as.matrix(X[,(ceiling(p/2)+1):p])%*%gamma.mar.x), FUN= function(pr) rbinom(n=1, size=1, prob = pr))
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
  
  # TREATMENT ASSIGNMENT
  X.tmp <- X
  if (cit){
    X.tmp[idx_NA] <- 0
  }
  
  if (class.interaction){
    assignment <- apply(cbind(X.tmp, class), MARGIN=1, FUN = function(x) gen_treat_lat(x, class.interaction, link = link))
    ps <- assignment[1,]
    treat <- assignment[2,]
  } else {
    assignment <- apply(X.tmp, MARGIN=1, FUN = function(x) gen_treat_lat(x, FALSE, link = link))
    ps <- assignment[1,]
    treat <- assignment[2,]
  }
  
  X.tmp <- X
  if (cio){
    X.tmp[idx_NA] <- 0
  }
  
  if (class.interaction){
    y <- apply(cbind(X.tmp, treat, class), MARGIN=1, FUN = function(x) gen_out_lat(x, class.interaction, sd=sd, link = link))
  } else {
    y <- apply(cbind(X.tmp, treat), MARGIN=1, FUN = function(x) gen_out_lat(x, FALSE, sd=sd, link = link))
  }
  
  if (class.interaction){
    y.1 <- apply(cbind(X, rep(1, n), class), MARGIN=1, FUN = function(x) gen_out_lat(x, class.interaction, sd=sd, link = link))
    y.0 <- apply(cbind(X, rep(0, n), class), MARGIN=1, FUN = function(x) gen_out_lat(x, class.interaction, sd=sd, link = link))
  } else {
    y.1 <- apply(cbind(X, rep(1, n)), MARGIN=1, FUN = function(x) gen_out_lat(x, FALSE, sd=sd, link = link))
    y.0 <- apply(cbind(X, rep(0, n)), MARGIN=1, FUN = function(x) gen_out_lat(x, FALSE, sd=sd, link = link))
  }
  tau <- y.1 - y.0
  
  return(list("X" = X,
              "ps" = ps,
              "treat" = treat,
              "y" = y,
              "X.incomp" = X.incomp,
              "tau" = tau,
              "class" = class))
}


gen_multisvd <- function(n, p, ngr = 5, ncpW = 2, ncpB = 2,
                       class.interaction = FALSE,
                       sd=0.1, 
                       seed = 0,
                       missing=FALSE, prop.missing = 0, 
                       cit = FALSE, cio = FALSE,
                       link = "nonlinear"){
  
  set.seed(seed)
  gamma.mar.x <- (20/p)*(runif(n=floor(p/2))-0.5)


  effective <- rep(ceiling(n/ngr),ngr)
  group <- rep(1:ngr, effective)
  # between part
  Xb <- matrix(rnorm(ngr*p), nrow=ngr, ncol=p)
  decomp <- svd(Xb)
  d <- diag(decomp$d[1:ncpB])
  Xb <- decomp$u[,1:ncpB]%*%d%*%t(decomp$v[,1:ncpB])
  Xb <- scale(Xb)
  Xb <- lapply(1:ngr, function(i) matrix(rep(Xb[i,],effective[i]), byrow=T, nrow=effective[i]))
  Xb <- do.call(rbind, Xb)
  
  # within part 
  I <- sum(effective)
  Xw = matrix(rnorm(I*p), nrow=I, ncol=p)
  decomp = svd(Xw)
  d = diag(decomp$d[1:ncpW])
  Xw = decomp$u[,1:ncpW]%*%d%*%t(decomp$v[,1:ncpW])
  Xw = scale(Xw)
  X <- Xb + Xw
  X = X  + matrix(rnorm(I*p, 0, sd), I, p)
  X.class = cbind.data.frame(as.factor(group), X)
  
  n <- nrow(X)
  p <- ncol(X)
  
  
  # MISSING DATA
  X.incomp <- X
  idx_NA <- NULL
  if (missing == "MCAR"){
    idx_NA = matrix(runif(n*dim(X)[2]), nrow=n, ncol=dim(X)[2]) <= prop.missing
    idx_NA[rowSums(idx_NA)==ncol(X), sample(ncol(X),1)] = FALSE # avoid having empty observations
    X.incomp[idx_NA] <- NA
  }
  if (mechanism == "MAR"){
    idx_NA <- c()
    for (j in 1:ceiling(p/2)) {
      m <- sapply(expit(2 + X[,(ceiling(p/2)+1):p]%*%gamma.mar.x), FUN= function(pr) rbinom(n=1, size=1, prob = pr))
      X.incomp[which(m==1), j] <- NA
      idx_NA <- cbind(idx_NA, m==1)
    }
    idx_NA <- cbind(idx_NA, matrix(FALSE, nrow = n, ncol = floor(p/2)))
    idx_NA <- matrix(idx_NA, nrow=n, ncol = p)
        
  } 
  if (missing == "MNAR"){
    
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
  
  X.incomp.class = cbind.data.frame(as.factor(group), X.incomp)
  
  # TREATMENT ASSIGNMENT
  X.tmp <- X
  if (cit){
    X.tmp[idx_NA] <- 0
  }
  
  if (class.interaction){
    assignment <- apply(cbind(X.tmp, as.factor(group)), MARGIN=1, FUN = function(x) gen_treat_svd(x, class.interaction, link = link))
    ps <- assignment[1,]
    treat <- assignment[2,]
  } else {
    assignment <- apply(X.tmp, MARGIN=1, FUN = function(x) gen_treat_svd(x, FALSE, link = link))
    ps <- assignment[1,]
    treat <- assignment[2,]
  }
  
  X.tmp <- X
  if (cio){
    X.tmp[idx_NA] <- 0
  }
  
  if (class.interaction){
    y <- apply(cbind(X.tmp, treat, as.factor(group)), MARGIN=1, FUN = function(x) gen_out_svd(x, class.interaction, sd=sd, link = link))
  } else {
    y <- apply(cbind(X.tmp, treat), MARGIN=1, FUN = function(x) gen_out_svd(x, FALSE, sd=sd, link = link))
  }
  
  if (class.interaction){
    y.1 <- apply(cbind(X, rep(1, n), as.factor(group)), MARGIN=1, FUN = function(x) gen_out_svd(x, class.interaction, sd=sd, link = link))
    y.0 <- apply(cbind(X, rep(0, n), as.factor(group)), MARGIN=1, FUN = function(x) gen_out_svd(x, class.interaction, sd=sd, link = link))
  } else {
    y.1 <- apply(cbind(X, rep(1, n)), MARGIN=1, FUN = function(x) gen_out_svd(x, FALSE, sd=sd, link = link))
    y.0 <- apply(cbind(X, rep(0, n)), MARGIN=1, FUN = function(x) gen_out_svd(x, FALSE, sd=sd, link = link))
  }
  tau <- y.1 - y.0
  
  return(list("X" = X,
              "ps" = ps,
              "treat" = treat,
              "y" = y,
              "X.incomp" = X.incomp,
              "tau" = tau,
              "class" = as.factor(group)))
}


gen_dlvm <- function(n, p, d=3, h = 5,
                     sd=0.1, 
                     seed = 0,
                     mechanism=FALSE, prop.missing = 0, 
                     cit = FALSE, cio = FALSE,
                     link = "nonlinear",
                     sigma.structure = "diagonal"){
  set.seed(0)
  gamma.mar.x <- (20/p)*(runif(n=floor(p/2))-0.5)


  V <- mvrnorm(n=p, mu= rep(1,h), Sigma = diag(h))
  W <- matrix(runif(n = h*d), ncol = d, nrow = h)
  a <- runif(n=h)
  b <- rnorm(n=p)
  
  if (sigma.structure == "diagonal"){
    alpha <- rnorm(n = h)
    beta <- runif(n = 1)
  } else {
    alpha <- mvrnorm(n = p, mu = rep(0,h), Sigma = diag(h))
    beta <- runif(n = p)
    U <- randortho(n = p, type="orthonormal")
  }
  set.seed(seed)
  
  get_params <- function(z){
    hu<-W%*%z + a

    if (sigma.structure == "diagonal"){
      mu <- (V%*%tanh(hu) + b)
      sig <- as.numeric(exp(t(alpha)%*%tanh(hu) + beta))
      Sigma <- sig*diag(p)
    } else {
      mu <- U%*%(V%*%tanh(hu) + b)
      Sigma <- U%*%diag(array(exp(alpha%*%tanh(hu) + beta)), nrow = p)%*%t(U)
    }
    return(list(mu=mu, Sigma = Sigma))
  }
  
  
  codes <- mvrnorm(n=n, mu = rep(0, d), Sigma = diag(d))
  
  params <- apply(codes, FUN = get_params, MARGIN = 1)
  
  X.list <- lapply(params, FUN = function(eta) mvrnorm(n=1, mu = eta$mu, Sigma = eta$Sigma))
  X <- do.call(rbind, X.list)
  
  X <- data.frame(X)
  
  # MISSING DATA
  X.incomp <- X
  idx_NA <- NULL
  if (mechanism == "MCAR"){
    idx_NA = matrix(runif(n*dim(X)[2]), nrow=n, ncol=dim(X)[2]) <= prop.missing
    idx_NA[rowSums(idx_NA)==ncol(X), sample(ncol(X),1)] = FALSE # avoid having empty observations
    X.incomp[idx_NA] <- NA
  }
  if (mechanism == "MAR"){
    idx_NA <- c()
    for (j in 1:ceiling(p/2)) {
      m <- sapply(expit(2 + X[,(ceiling(p/2)+1):p]%*%gamma.mar.x), FUN= function(pr) rbinom(n=1, size=1, prob = pr))
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
  
  # TREATMENT ASSIGNMENT
  X.tmp <- X
  if (cit){
    X.tmp[idx_NA] <- 0
  }
  
  
  assignment <- apply(X.tmp, MARGIN=1, FUN = function(x) gen_treat_deep(x, link = link))
  ps <- assignment[1,]
  treat <- assignment[2,]

  
  X.tmp <- X
  if (cio){
    X.tmp[idx_NA] <- 0
  }
  
  y <- apply(cbind(X.tmp, treat), MARGIN=1, FUN = function(x) gen_out_deep(x, sd=sd, link = link))
  
  
  y.1 <- apply(cbind(X, rep(1, n)), MARGIN=1, FUN = function(x) gen_out_deep(x, sd=sd, link = link))
  y.0 <- apply(cbind(X, rep(0, n)), MARGIN=1, FUN = function(x) gen_out_deep(x, sd=sd, link = link))
  
  tau <- y.1 - y.0
  
  return(list("X" = X,
              "ps" = ps,
              "treat" = treat,
              "y" = y,
              "X.incomp" = X.incomp,
              "tau" = tau,
              "class" = NULL))
}


gen_ding <- function(n, set = 1,
                     seed = 0,
                     missing = "MNAR"){
  set.seed(seed)
  
  if (set == 1){
    # DATA
    X <- data.frame(cbind(rnorm(n=n, mean = 1, sd = 1), 
                          rbinom(n = n, size = 1, prob = 0.5)))
    
    # TREATMENT ASSIGNMENT
    ps <- apply(X, MARGIN=1, FUN = function(x) expit(1.25 -0.5*x[1] - 0.5*x[2]))
    treat <- sapply(ps, FUN= function(p) rbinom(n=1, size=1, prob = p))
    
    
    
    y.1 <- apply(X, MARGIN=1, FUN = function(x) 3*x[1]+2*x[2] + rnorm(1,0,1))
    y.0 <- apply(X, MARGIN=1, FUN = function(x) 0.5 + 2*x[1] + x[2] + rnorm(1,0,1))
    
    y <- treat*y.1 + (1-treat)*y.0
    
    tau <- 1
  
    # MISSING DATA
    X.incomp <- X
  
    p1 <- expit(-2+2*X[,1] + treat*(1.5+X[,2]))
    idx_NA <- matrix(cbind(sapply(p1, FUN = function(p) rbinom(n=1, size=1, prob = p))==1, 
                           rep(FALSE, n)),
                     nrow = n, ncol = 2)
    X.incomp[idx_NA] <- NA
  }
    
  return(list("X" = X,
              "ps" = ps,
              "treat" = treat,
              "y" = y,
              "X.incomp" = X.incomp,
              "tau" = rep(1,n),
              "class" = NULL))
}


