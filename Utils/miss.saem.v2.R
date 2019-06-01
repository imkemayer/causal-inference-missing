miss.saem.v2 <- function(X.obs,y,pos_var=1:ncol(X.obs),maxruns=500,tol_em=1e-7,nmcmc=2,tau=1,k1=50, seed=200, print_iter=TRUE, var_cal=FALSE, ll_obs_cal=FALSE) {
  set.seed(seed)
  #judge
  if (class(X.obs) == "data.frame") {
    X.obs <- as.matrix(X.obs)
  }
  if (sum(sapply(X.obs, is.numeric)) < ncol(X.obs)) {
    stop("Error: the variables should be numeric.")
  }
  if (sum(y==1) +  sum(y==0) < nrow(X.obs)) {
    stop("Error: y must be coded by 0 or 1, and there is no missing data in y.")
  }
  
  if (sum(pos_var %in% 1:ncol(X.obs)) < length(pos_var))  {
    stop("Error: index of selected variables must be in the range of covariates.")
  }
  
  if (length(unique(pos_var)) != length(pos_var)){
    stop("Error: index of selected variables must not be repeated.")
  }
  
  p=ncol(X.obs)
  
  #delete rows completely missing
  if(any(apply(is.na(X.obs),1,sum)==p)){
    i_allNA=which(apply(is.na(X.obs),1,sum)==p)
    X.obs = X.obs[-i_allNA,]
    y = y[-i_allNA]
  }
  if(any((is.na(y))==TRUE)){
    i_YNA=which(is.na(y)==TRUE)
    X.obs = X.obs[-i_YNA,]
    y = y[-i_YNA]
  }
  n=length(y)
  
  
  rindic = as.matrix(is.na(X.obs))
  if(sum(rindic)>0){
    whichcolmissing = (1:ncol(rindic))[apply(rindic,2,sum)>0]
    missingcols = length(whichcolmissing)
  }
  if(sum(rindic)==0){missingcols=0}
  
  
  ptm <- Sys.time()
  if(missingcols>0){
    k=0
    cstop=0.1
    seqbeta = matrix(NA,nrow=ncol(X.obs)+1,ncol=(maxruns+1))
    seqbeta_avg = matrix(NA,nrow=ncol(X.obs)+1,ncol=(maxruns+1))
    
    X.mean = X.obs
    for(i in 1:ncol(X.mean)){
      X.mean[is.na(X.mean[,i]), i] <- mean(X.mean[,i], na.rm = TRUE)
    }
    X.sim <- X.mean
    
    mu = apply(X.mean,2,mean)
    Sigma = var(X.mean)*(n-1)/n
    beta= rep(0,p+1)
    beta[c(1,pos_var+1)]= glm(y~ X.mean[,pos_var],family=binomial(link='logit'))$coef
    
    while ((cstop>tol_em)*(k<maxruns)|(k<20)){
      k = k+1
      beta.old = beta
      
      if(k <k1){gamma <- 1}else{gamma <- 1/(k-(k1-1))^tau}
      
      S.inv <- solve(Sigma)
      
      for (i in (1:n)) {
        jna <- which(is.na(X.obs[i,]))
        njna <- length(jna)
        if (njna>0) {
          xi <- X.sim[i,]
          Oi <- solve(S.inv[jna,jna])
          mi <- mu[jna]
          lobs <- beta[1]
          if (njna<p) {
            jobs <- setdiff(1:p,jna)
            mi <- mi - (xi[jobs] - mu[jobs])%*%S.inv[jobs,jna]%*%Oi
            lobs <- lobs + sum(xi[jobs]*beta[jobs+1])
          }
          
          cobs <- exp(lobs)
          if(cobs<1e-16){cobs=1e-16}
          if(cobs>1e+16){cobs=1e+16}
          
          xina <- xi[jna]
          betana <- beta[jna+1]
          for (m in (1:nmcmc)) {
            xina.c <- mi + rnorm(njna)%*%chol(Oi)
            
            if (y[i]==1){
              exp1 <- 1+exp(-sum(xina*betana))
              exp2 <- 1+exp(-sum(xina.c*betana))
              if(exp1<1e-16){exp1=1e-16}
              if(exp2<1e-16){exp2-1e-16}
              if(exp1>1e+16){exp1=1e+16}
              if(exp2>1e+16){exp2=1e+16}
              alpha <- (exp1/cobs)/(exp2/cobs)
            } else {
              exp1 <- 1+exp(sum(xina*betana))
              exp2 <- 1+exp(sum(xina.c*betana))
              if(exp1<1e-16){exp1=1e-16}
              if(exp2<1e-16){exp2-1e-16}
              if(exp1>1e+16){exp1=1e+16}
              if(exp2>1e+16){exp2=1e+16}
              alpha <- (exp1*cobs)/(exp2*cobs)
            }
            if (runif(1) < alpha){
              xina <- xina.c
            }
          }
          X.sim[i,jna] <- xina
        }
      }
      beta_new= rep(0,p+1)
      beta_new[c(1,pos_var+1)]= glm(y~ X.sim[,pos_var],family=binomial(link='logit'))$coef
      
      beta <- (1-gamma)*beta + gamma*beta_new
      cstop = sum((beta-beta.old)^2)
      
      mu <- (1-gamma)*mu + gamma*colMeans(X.sim)
      Sigma <- (1-gamma)*Sigma + gamma*cov(X.sim)
      
      seqbeta[,k]=beta.old
      
      if(k==1){
        seqbeta_avg[,k]=beta.old
      }else{
        seqbeta_avg[,k]= 1/k*rowSums(seqbeta[,1:k])
      }
      
      if(print_iter==TRUE & k %% 10 == 0){
        cat(sprintf('iteration = %i ', k))
        cat(sprintf('beta ='),beta,'\n')
        cat(sprintf('Distance from last iteration ='),cstop,'\n')
      }
    }
    var_obs = ll = std_obs =NULL
    if(var_cal==TRUE){
      var_obs = louis_lr_saem(beta,mu,Sigma,y,X.obs,pos_var,rindic,whichcolmissing,mc.size=1000)
      std_obs <- sqrt(diag(var_obs))
    }
    if(ll_obs_cal==TRUE){
      ll= likelihood_saem(beta,mu,Sigma,y,X.obs,rindic,whichcolmissing,mc.size=1000)
    }
  }
  if(missingcols==0){
    X.obs = matrix(X.obs,nrow=n)
    data.complete <- data.frame(y=y,X.obs)
    model.complete <- glm(y ~.,family=binomial(link='logit'),data=data.complete)
    mu = apply(X.obs,2,mean)
    Sigma = var(X.obs)*(n-1)/n
    beta <- model.complete$coefficients
    var_obs = ll = ll1 =ll2= std_obs =seqbeta_avg= seqbeta=NULL
    if(var_cal==TRUE){
      P <- predict(model.complete, type = "response")
      W <- diag(P*(1-P))
      X <- model.matrix(model.complete)
      
      var_obs <- solve(t(X)%*%W%*%X)
      std_obs <- sqrt(diag(var_obs))
    }
    if(ll_obs_cal==TRUE){
      ll = likelihood_saem(beta,mu,Sigma,y,X.obs,rindic,whichcolmissing,mc.size=1000)
    }
  }
  time_run=Sys.time() - ptm
  return(list(mu=mu, sig2=Sigma, beta=beta,time_run=time_run,seqbeta=seqbeta,seqbeta_avg=seqbeta_avg,ll=ll,var_obs=var_obs,std_obs=std_obs))
}
