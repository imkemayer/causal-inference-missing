###################################################################################
## Define propensity scores estimation
###################################################################################


# Use SAEM for logistic regression treatment ~ X with missing values in X
predict_misaem <- function(X, treat, seed=0, pattern = NULL, use.interaction = FALSE){
  length.covariates <- dim(X)[2]
  col.complete <- which(sapply(X, FUN=function(x) sum(is.na(x))==0))
  length.mask <- 0
  
  if (!is.null(pattern)){
    colname <- colnames(X)
    X <- cbind(X, apply(pattern, FUN=function(x) 2*as.numeric(x) - 1, MARGIN=c(1,2)))
    length.mask <- dim(pattern)[2]
    if (dim(pattern)[2]>0){
      colnames(X) <- c(colname, paste0("R",1:dim(pattern)[2]))
    }
  }
  
  if (use.interaction & !is.null(pattern) & length.mask > 0 & length(col.complete) > 0){
    rnam <- paste("R", 1:length.mask, sep="")
    xnam1 <- paste("X", 1:length.covariates, sep="")
    xnam2 <- paste("X", col.complete, sep="")
    fmla <- as.formula(paste("~ ", paste0(paste(c(xnam1,rnam), collapse= "+"), # "-1+",
                                          "+",
                                          paste(kronecker(xnam2, rnam,function(x,r) paste0(x,"*",r)), collapse= "+"))))
    na.action.default <- getOption("na.action")
    options(na.action = "na.pass")
    X <- model.matrix(fmla, data = X)
    options(na.action = na.action.default)
  } 
  
  test.boundary <- T
  k <- 0 
  while (test.boundary & k < 10){
    list.saem <- miss.saem.v2(as.matrix(X), treat, seed = seed+k, tol_em=1e-06, print_iter=FALSE)
    
    pr.saem <- NULL
    try(pr.saem <- pred_saem(as.matrix(X), list.saem$beta, list.saem$mu, list.saem$sig2, method="map"))
    if (is.null(pr.saem)){
      pr.saem <- pred_saem(as.matrix(X), list.saem$beta, list.saem$mu, list.saem$sig2, method="impute")
    }
    test.boundary <- (min(pr.saem)<1e-6 | max(pr.saem)>1-1e-6)
    k <- k+1
  }
  
  fitted = as.data.frame(pr.saem)
  colnames(fitted) <- c("pscore")
  
  
  return(list(fitted=fitted, mu = list.saem$mu, sigma = list.saem$sig2))
}

# Standard logistic regression on complete X
predict_glm <- function(X, treat, seed, family = binomial, regularize=FALSE){
  set.seed(seed)
  
  if (regularize){
    x <- model.matrix(~., data = X)
    cv.fit <- glmnet::cv.glmnet(x=x, 
                                y=treat, alpha = 0,
                                family = "binomial")
    glm_mod <- glmnet::glmnet(x=x, y=treat, alpha = 0, family = "binomial",
                              lambda = cv.fit$lambda.min)
    pred <- predict(glm_mod, newx = x, type="response")
  } else {
    glm_mod <- train(y = treat,
                     x = X,
                     trControl = trainControl(method = "cv", number = 5),
                     method = "glm", 
                     family="binomial")
    pred <- predict(glm_mod$finalModel, X, type="response")
  }
  fitted = as.data.frame(pred)
  colnames(fitted) <- c("pscore")
  return(fitted)
}

# Random forest (from ranger package) on complete X
predict_random_forests <- function(X, treat, seed){
  set.seed(seed)
  
  rf_model <- ranger(treat ~ . , data = data.frame(cbind(X, treat=treat)), probability = TRUE)
  fitted = as.data.frame(rf_model$predictions[,2])
  colnames(fitted) <- c("pscore")
  return(fitted)
}

# Gradient Boosting on complete X
predict_gbm <- function(X, treat, seed){
  set.seed(seed)
  
  gbm.r <- train(y= treat,
                 x=X, method="gbm",
                 trControl = trainControl(method="cv",number=5), verbose = FALSE)
  
  pred <- predict(gbm.r, data = data.frame(X), type = "prob")[,2]
  fitted = as.data.frame(pred)
  colnames(fitted) <- c("pscore")
  return(fitted)
}

# Regression forest (of grf package) on complete X
predict_grf <- function(X, treat, seed){
  set.seed(seed)
  na.action <- options()$na.action
  options(na.action='na.pass')
  if (is.data.frame(X)){
    X.m = model.matrix(~., data=X)
  } else {
    X.m = model.matrix(~., data=data.frame(X))
  }
  if (max(as.integer(treat))==2){
    treat=as.integer(treat)-1
  } else {
    treat = as.integer(treat)
  }
  
  forest.W = regression_forest(X.m, treat, tune.parameters = "all")
  pred = predict(forest.W)$predictions
  fitted = as.data.frame(pred)
  colnames(fitted) <- c("pscore")
  options(na.action=na.action)
  return(fitted)
}

# Regression forest (of ranger package) on complete X
predict_ranger <- function(X, treat, seed){
  set.seed(seed)
  df <- data.frame(cbind(X, "w" = treat))
  pred <- ranger::ranger(w ~., data=df,probability = TRUE)$predictions[,2]
  fitted = as.data.frame(pred)
  colnames(fitted) <- c("pscore")
  return(fitted)
}

###################################################################################
## Define IPW method
###################################################################################

ipw <- function(X, outcome, treat, 
                ps.method="glm", 
                target= "all", 
                seed = 0,
                trimming_weight = 1,
                fitted=NULL,
                mask= NULL,
                use.interaction = FALSE,
                regularize = FALSE){
  # @param X [data.frame] confounders \eqn{n \times p}{n * p}
  # @param outcome Response vector \eqn{n \times 1}{n * 1}
  # @param treat Binary treatment assignment vector \eqn{n \times 1}{n * 1}
  # @param ps.method either one of "glm", "rf", "gbm", "grf"
  # @param target either one of "all", "treated", "control", "overlap", "matching"
  # @param seed random values generator seed (optional)
  # @param trimming_weight quantile used for trimming weights (if set to 1, then no trimming, if set to 0.9, then trimming above 0.9 quantile and below 0.1 quantile) (optinal)  
  # @param fitted [data.frame] Data.frame with one variable "pscore" containing (estimated) propensity scores \eqn{n \times 1}{n * 1} (optional)
  # @param mask [data.frame] response pattern if it is to be added to the confounders for the propensity model (optional)
  # @param use.interaction [boolean] whether to add interactions between fully observed confounders and the mask (optional)

  # @return The average treatment effect estimate and its sample standard deviation

  ##################
  # Propensity model
  ##################
  length.covariates <- dim(X)[2]
  col.complete <- which(sapply(X, FUN=function(x) sum(is.na(x))==0))
  length.mask <- 0
  
  # Add response mask
  if (!is.null(mask)){
    colname <- colnames(X)
    X <- cbind(X, mask)
    length.mask <- dim(mask)[2]
  } 
  if (length.mask>0){
    colnames(X) <- c(colname, paste0("R",1:length.mask))
  }
  
  # Add interaction terms between fully observed X and mask
  if (use.interaction & !is.null(mask) & length.mask > 0 & length(col.complete) > 0){
    xnam1 <- paste("X", 1:length.covariates, sep="")
    xnam2 <- paste("X", col.complete, sep="")
    rnam <- paste("R", 1:length.mask, sep="")
    fmla <- as.formula(paste("~ ", paste0(paste(c(xnam1,rnam), collapse= "+"),
                                          "+",
                                          paste(kronecker(xnam2, rnam,function(x,r) paste0(x,"*",r)), collapse= "+"))))
    na.action.default <- getOption("na.action")
    options(na.action = "na.pass")
    X <- data.frame(model.matrix(fmla, data = X))
    options(na.action = na.action.default)
  }
  
  # Estimate propensity scores
  if (is.null(fitted)){
    if (ps.method == "glm") {
      fitted <- predict_glm(X, as.factor(treat), seed, regularize=regularize)
    } else if (ps.method == "rf") {
      fitted <- predict_random_forests(X, as.factor(treat), seed)
    } else if (ps.method == "gbm") {
      fitted <- predict_gbm(X, as.factor(treat), seed)
    } else if (ps.method %in% c("grf","grf.ate")) {
      fitted <- predict_grf(X, as.factor(treat), seed)
    } else if (ps.method == "ranger"){
      fitted <- predict_ranger(X, as.factor(treat), seed)
    }
  }
  
  # Compute weights depending on the estimand
  if (is.numeric(treat)){
    W <- as.logical(treat)
  } else {
    W <- as.logical(as.character(treat))
  }
  if (target == "all"){
    fitted$weight <- (W)/fitted$pscore + (1 - W)/ (1 - fitted$pscore)
  } else if (target == "treated") {
    fitted$weight <- W + (1-W)*fitted$pscore / (1 - fitted$pscore)
  } else if (target == "control") {
    fitted$weight <- W * (1 - fitted$pscore) / fitted$pscore + (1-W)
  } else if (target == "overlap") {
    fitted$weight <- W * (1 - fitted$pscore) + (1-W)*fitted$pscore
  } else if (target == "matching") {
    fitted$weight <- W * apply(cbind(fitted$pscore, 1-fitted$pscore), 1, min) / fitted$pscore + 
                      (1-fitted$treated)* apply(cbind(fitted$pscore, 1-fitted$pscore), 1, min) / (1-fitted$pscore)
  }
  
  
  
  # Trim the weights (default: no trimming)
  if (any(is.na(fitted$weight)) | any(is.infinite(fitted$weight))){
    weightMax <- quantile(fitted$weight, c(trimming_weight), na.rm=TRUE) 
    fitted[is.na(fitted$weight) ,"weight"] <- weightMax # If the weight is NA this means that pscore was 0
    fitted[fitted$weight > weightMax,"weight"] <- weightMax
    weightMin <- quantile(fitted$weight, c(1-trimming_weight), na.rm=TRUE) 
    fitted[fitted$weight < weightMin,"weight"] <- weightMin
  }
  
  
  # Compute HT verion of IPW
  if (!is.numeric(outcome)){
    if (length(unique(outcome))==2){
      Y <- as.logical(as.character(outcome))
    } else
      Y <- as.numeric(as.character(outcome))
  } else{
    Y <- outcome
  }
  ipw1 <- 1/length(Y)*(sum(Y[which(W==1)]*fitted$weight[which(W==1)]) - sum(Y[which(!(W==1))]*fitted$weight[which(!(W==1))]))
  
  # Compute normalized version of IPW
  ipw2 <- 1/sum(fitted$weight[which(W==1)]) * sum(Y[which(W==1)] * fitted$weight[which(W==1)]) - 1/sum(fitted$weight[which(!(W==1))]) * sum(Y[which(!(W==1))] * fitted$weight[which(!(W==1))])
  
  delta_i <- 1/sum(fitted$weight[which(W==1)]) * W*(Y * fitted$weight) - 1/sum(fitted$weight[which(W==0)]) * (1-W)*(Y * fitted$weight)
  se.ipw2 <- sqrt(var(length(Y)*delta_i) / (length(Y) - 1))
  
  return(cbind(ipw1 = ipw1,
               ipw2 = ipw2,
               se.ipw2 = se.ipw2))
}




###################################################################################
## Define DR method
###################################################################################


dr <- function(X, 
               X.for.ps=NULL,
               X.for.outcome=NULL, 
               outcome, treat, 
               ps.method="glm", 
               target= "all", 
               seed = 0,
               trimming_weight = 1,
               fitted=NULL,
               out.method = "glm",
               mask = NULL,
               mask.for.ps = NULL,
               mask.for.outcome = NULL,
               use.interaction = FALSE,
               subset = NULL,
               regularize = FALSE,
               clusters=NULL){

  #' @param X [data.frame] confounders \eqn{n \times p}{n * p}
  #' @param X.for.ps [data.frame] confounders and other predictors of treatment assignment \eqn{n \times (p + p_p)}{n * (p+pp)}
  #' @param X.for.outcome [data.frame] confounders and other predictors of outcome \eqn{n \times (p + p_o)}{n * (p+po)}
  #' @param outcome Response vector \eqn{n \times 1}{n * 1}
  #' @param treat Binary treatment assignment vector \eqn{n \times 1}{n * 1}
  #' @param ps.method either one of "glm", "rf", "gbm", "grf", "saem"
  #' @param target either one of "all", "treated", "control", "overlap", "matching"
  #' @param seed random values generator seed (optional)
  #' @param trimming_weight quantile used for trimming weights (if set to 1, then no trimming, if set to 0.9, then trimming above 0.9 quantile and below 0.1 quantile) (optinal)  
  #' @param fitted [data.frame] Data.frame with one variable "pscore" containing (estimated) propensity scores \eqn{n \times 1}{n * 1} (optional)
  #' @param out.method either one of "glm", "rf", "gbm", "grf", "grf.ate"
  #' @param mask [data.frame] response pattern if it is to be added to the confounders for the propensity model (optional)
  #' @param mask.for.prop [data.frame] response pattern if it is to be added to the propensity model (optional)
  #' @param mask.for.outcome [data.frame] response pattern if it is to be added to the outcome model (optional)
  #' @param use.interaction [boolean] whether to add interactions between fully observed confounders and the mask (optional)
  #' @param regularize [boolean] whether to regularize logistic or linear regressions
  #' 
  #' @return The average treatment effect estimate and an estimate of its standard deviation
  
  ##################
  # Propensity model
  ##################
  X.orig <- X
  mask.orig <- mask
  
  if (!(is.null(X.for.ps))){
    X <- X.for.ps
    mask <- mask.for.ps
  }
  
  length.covariates <- dim(X)[2]
  if (!is.null(mask)){
    col.complete <- which(sapply(mask, FUN=function(x) length(unique(x))==1))
  } else {
    col.complete <- c()
  }
  length.mask <- 0
  
  # Add response mask
  if (is.null(mask)){
    X1 <- X
  } else {
    colname <- colnames(X)
    length.mask <- dim(mask)[2]
    X1 <- cbind(X, mask)
    if (length.mask>0){
      colnames(X1) <- c(colname, paste0("R",1:length.mask))
    }
  }
  
  # Add interaction terms between fully observed confounders and mask
  if (use.interaction & !is.null(mask) & length.mask > 0 & length(col.complete) > 0){
    xnam1 <- paste("X", 1:length.covariates, sep="")
    xnam2 <- paste("X", col.complete, sep="")
    rnam <- paste("R", 1:length.mask, sep="")
    fmla <- as.formula(paste("~ ", paste0(paste(c(xnam1,rnam), collapse= "+"),
                                          "+",
                                          paste(kronecker(xnam2, rnam,function(x,r) paste0(x,"*",r)), collapse= "+"))))
    na.action.default <- getOption("na.action")
    options(na.action = "na.pass")
    X1 <- data.frame(model.matrix(fmla, data = X1))
    options(na.action = na.action.default)
  }
  
  # Estimate propensity scores
  if (is.null(fitted)) {
    if (ps.method == "glm") {
      fitted <- predict_glm(X1, as.factor(treat), seed, regularize=regularize)
    } else if (ps.method == "rf") {
      fitted <- predict_random_forests(X1, as.factor(treat), seed)
    } else if (ps.method == "gbm") {
      fitted <- predict_gbm(X1, as.factor(treat), seed)
    } else if (ps.method %in% c("grf")) {
      fitted <- predict_grf(X1, as.factor(treat), seed)
    }
  }
  # Compute weights depending on the estimand
  if (is.numeric(treat)){
    W <- as.logical(treat)
  } else {
    W <- as.logical(as.character(treat))
  }
  if (!is.null(fitted)){
    if (target == "all"){
      fitted$weight <- (W)/fitted$pscore + (1 - W)/ (1 - fitted$pscore)
    } else if (target == "treated") {
      fitted$weight <- W + (1-W)*fitted$pscore / (1 - fitted$pscore)
    } else if (target == "control") {
      fitted$weight <- W* (1 - fitted$pscore) / fitted$pscore + (1-W)
    } else if (target == "overlap") {
      fitted$weight <- W * (1 - fitted$pscore) + (1-W)*fitted$pscore
    } else if (target == "matching") {
      fitted$weight <- W * apply(cbind(fitted$pscore, 1-fitted$pscore), 1, min) / fitted$pscore + (1-fitted$treated)* apply(cbind(fitted$pscore, 1-fitted$pscore), 1, min) / (1-fitted$pscore)
    }
  }
  
  if (!is.null(fitted) & (any(is.nan(fitted$weight)) | any(is.infinite(fitted$weight)))) {
    # Trim the weights (default: no trimming)
    weightMax <- quantile(fitted$weight, c(trimming_weight), na.rm=TRUE) 
    fitted[is.na(fitted$weight) ,"weight"] <- weightMax # If the weight is NA this means that pscore was 0
    fitted[fitted$weight > weightMax,"weight"] <- weightMax
    weightMin <- quantile(fitted$weight, c(1-trimming_weight), na.rm=TRUE) 
    fitted[fitted$weight < weightMin,"weight"] <- weightMin
  }

  ###############
  # Outcome model
  ###############
  if (!is.null(fitted)){
    W.hat <- fitted$weight
  }

  X <- X.orig
  mask <- mask.orig
  
  if (!is.null(X.for.outcome)){
    X <- X.for.outcome
    mask <- mask.for.outcome
  }
  length.covariates <- dim(X)[2]
  if (!is.null(mask)){
    col.complete <- which(sapply(mask, FUN=function(x) length(unique(x))==1))
  } else{
    col.complete <- c()
  }
  length.mask <- 0

  colnames(X) <- paste("X", 1:dim(X)[2], sep="")
  X2 <- X
  if (!is.null(mask)) {
    X2 <- cbind(X, mask)
    length.mask <- dim(mask)[2]
  }
  if (length.mask>0){
    colnames(X2) <- c(colnames(X), paste0("R",1:length.mask))
  }

  # Use grf regression forests to fit propensity and outcome models and then use AIPW formula for treatment effect estimation
  if (ps.method == "grf" & out.method == "grf"){
    na.action <- options()$na.action
    options(na.action='na.pass')
    var.factor <- colnames(X2)[sapply(X2, is.factor)]
    X.1.m = model.matrix(~. -1, data=data.frame(X2[which(treat==1),]), 
                         contrasts = as.list(setNames(rep("contr.sum", length(var.factor)), var.factor)))
    X.0.m = model.matrix(~. -1, data=data.frame(X2[which(treat==0),]),
                         contrasts = as.list(setNames(rep("contr.sum", length(var.factor)), var.factor)))
    
    X.m = model.matrix(~. -1, data=data.frame(X2), 
                       contrasts = as.list(setNames(rep("contr.sum", length(var.factor)), var.factor)))
    
    forest.1.Y = regression_forest(X.1.m, outcome[which(treat==1)], tune.parameters = "all")
    y_1.hat = predict(forest.1.Y, X.m)$predictions
    
    forest.0.Y = regression_forest(X.0.m, outcome[which(treat==0)], tune.parameters = "all")
    y_0.hat = predict(forest.0.Y, X.m)$predictions
    
    if (!is.null(subset)){
      X2 <- X2[subset,]
      y_1.hat <- y_1.hat[subset]
      y_0.hat <- y_0.hat[subset]
      treat <- treat[subset]
      outcome <- outcome[subset]
      W.hat <- W.hat[subset]
    } 
    n.sample <- dim(X2)[1]
    delta_i <- y_1.hat -  y_0.hat + treat*(outcome-y_1.hat)*W.hat - (1-treat)*(outcome-y_0.hat)*W.hat
    treatment_effect_dr <- mean(delta_i)
    treatment_effect_dr <- c(treatment_effect_dr, sqrt(sum((delta_i-treatment_effect_dr)^2))/n.sample)
    options(na.action=na.action)
  } 

  # Use grf's average_treatment_effect function
  if (ps.method == "grf.ate" & out.method == "grf.ate") {  
    na.action <- options()$na.action
    options(na.action='na.pass')
    y.hat <- NULL
    if (!is.null(X.for.outcome)) {
      X.m = model.matrix(~. , data=data.frame(X2))
      forest.Y = regression_forest(X.m, outcome, tune.parameters = "all")
      y.hat = predict(forest.Y, X.m)$predictions
    }
    w.hat <- NULL
    if (!is.null(X.for.ps)) {
      X.m = model.matrix(~. , data=data.frame(X1))
      forest.W = regression_forest(X.m, treat, tune.parameters = "all")
      w.hat = predict(forest.W, X.m)$predictions
    }
    
    X <- X.orig
    mask <- mask.orig
    
    length.covariates <- dim(X)[2]
    if (!is.null(mask)){
      col.complete <- which(sapply(mask, FUN=function(x) length(unique(x))==1))
    } else{
      col.complete <- c()
    }
    length.mask <- 0
    
    colnames(X) <- paste("X", 1:dim(X)[2], sep="")
    Xconf <- X
    if (!is.null(mask)) {
      Xconf <- cbind(X, mask)
      length.mask <- dim(mask)[2]
    }
    if (length.mask>0){
      colnames(Xconf) <- c(colnames(X), paste0("R",1:length.mask))
    }
    X.m <- model.matrix(~. , data=data.frame(Xconf)) # Xconf corresponds to confounders
    tau.forest = causal_forest(X.m, outcome, treat, Y.hat = y.hat, W.hat =w.hat, clusters=clusters)
  
    treatment_effect_dr <- average_treatment_effect(tau.forest, target.sample = target, subset = subset)
    options(na.action=na.action)
  } 
  if (ps.method == "glm.grf" & out.method == "glm.grf") {  
    if (length(unique(outcome))==2){ y = as.factor(outcome) } else { y = outcome }
    if (regularize){
      x <- model.matrix(~., data = data.frame(X2))
      if (is.factor(y)){
        cv.fit <- glmnet::cv.glmnet(x=x, 
                                    y=y, alpha = 0,
                                    family = "binomial")
        lm.fit<- glmnet::glmnet(x=x, y=y, alpha = 0, family = "binomial",
                                     lambda = cv.fit$lambda.min)
        y.hat <- predict(lm.fit, newx=x, type="response")
      } else {
        cv.fit <- glmnet::glmnet(x=x, 
                                 y=y, alpha = 0, nlambda = 50)
        lm.fit <- glmnet::glmnet(x=x, y=y, alpha = 0,
                                     lambda = cv.fit$lambda.min)
        y.hat <- predict(lm.fit, newx=x)
      }
      
    } else { 
      if (is.factor(y)){
        lm.fit <- train(y=y,
                            x=X2, method="glm",
                            trControl = trainControl(method="cv",number=5),
                            family="binomial")
        y.hat <- predict(lm.fit, X2, type = "prob")[,2]
      } else {
        lm.fit <- train(y=y,
                        x=X2, method="glm",
                        trControl = trainControl(method="cv",number=5))
        y.hat <- predict(lm.fit, X2)
      }
    }
    
    if (regularize){
      x <- model.matrix(~., data = data.frame(X1))
      cv.fit <- glmnet::cv.glmnet(x=x, 
                                  y=as.factor(treat), alpha = 0,
                                  family = "binomial")
      lm.fit<- glmnet::glmnet(x=x, y=as.factor(treat), alpha = 0, family = "binomial",
                                lambda = cv.fit$lambda.min)
      w.hat <- predict(lm.fit, newx=x, type="response")
    } else { 
      lm.fit <- train(y=as.factor(treat),
                      x=X1, method="glm",
                      trControl = trainControl(method="cv",number=5),
                      family="binomial")
      w.hat <- predict(lm.fit, X1, type="prob")[,2]
    }
    
    X <- X.orig
    mask <- mask.orig
    length.covariates <- dim(X)[2]
    if (!is.null(mask)){
      col.complete <- which(sapply(mask, FUN=function(x) length(unique(x))==1))
    } else{
      col.complete <- c()
    }
    length.mask <- 0
    colnames(X) <- paste("X", 1:dim(X)[2], sep="")
    Xconf <- X
    if (!is.null(mask)) {
      Xconf <- cbind(X, mask)
      length.mask <- dim(mask)[2]
    }
    if (length.mask>0){
      colnames(Xconf) <- c(colnames(X), paste0("R",1:length.mask))
    }
    X.m <- model.matrix(~. , data=data.frame(Xconf)) # Xconf corresponds to confounders
    tau.forest = causal_forest(X.m, outcome, treat, Y.hat = y.hat, W.hat =w.hat, clusters=clusters)
    
    treatment_effect_dr <- average_treatment_effect(tau.forest, target.sample = target, subset = subset)
  }
  
  # Use EM for linear outcome model (either on complete(d) or incomplete X)
  if (out.method %in% c("glm", "gbm")){ # "gbm" is only available for complete data

    # Regular EM for linear regression Y ~ X to predict potential outcomes
    if (sum(is.na(X))==0){
      if (use.interaction & length.mask > 0 & length(col.complete) > 0){
        xnam1 <- paste("X", 1:length.covariates, sep="")
        xnam2 <- paste("X", col.complete, sep="")
        rnam <- paste("R", 1:length.mask, sep="")
        fmla <- as.formula(paste("~ ", paste0(paste(c(xnam1,rnam), collapse= "+"),
                                              "+",
                                              paste(kronecker(xnam2, rnam,function(x,r) paste0(x,"*",r)), collapse= "+"))))
        na.action.default <- getOption("na.action")
        options(na.action = "na.pass")
        X2 <- data.frame(model.matrix(fmla, data = X2))
        options(na.action = na.action.default)
      }

      df.treated <- data.frame(y = outcome[which(treat==1)], X = X2[which(treat==1),])
      colnames(df.treated) <- c("y", paste("X", 1:dim(X2)[2], sep=""))
      df.control <- data.frame(y = outcome[which(treat==0)], X = X2[which(treat==0),])
      colnames(df.control) <- c("y", paste("X", 1:dim(X2)[2], sep=""))
      
      df.treated <- df.treated[!duplicated(df.treated[,2:ncol(df.treated)]),]
      one.level.treated <- sapply(df.treated, FUN = function(x) length(unique(x))==1)
      df.treated <- df.treated[,!one.level.treated]
      fmla <- as.formula(paste0("y ~ ", paste( colnames(df.treated[,2:ncol(df.treated)]), collapse= "+")))
      if (out.method == "glm"){
        if (length(unique(df.treated[,1]))==2){ y = as.factor(df.treated[,1]) } else { y = df.treated[,1] }
        if (regularize){
          x <- model.matrix(~., data = data.frame(df.treated[,-1]))
          if (is.factor(y)){
            cv.fit <- glmnet::cv.glmnet(x=x, 
                                        y=y, alpha = 0,
                                        family = "binomial")
            lm.treated <- glmnet::glmnet(x=x, y=y, alpha = 0, family = "binomial",
                                         lambda = cv.fit$lambda.min)
          } else {
            cv.fit <- glmnet::glmnet(x=x, 
                                  y=y, alpha = 0, nlambda = 50)
            lm.treated <- glmnet::glmnet(x=x, y=y, alpha = 0,
                                         lambda = cv.fit$lambda.min)
          }
        } else {
          if (is.factor(y)){
            lm.treated <- train(y=y,
                                x=data.frame(df.treated[,-1]), method="glm", family="binomial",
                                trControl = trainControl(method="cv",number=5))
          } else {
            lm.treated <- train(y=y,
                                x=data.frame(df.treated[,-1]), method="glm",
                                trControl = trainControl(method="cv",number=5))
          }
        }
      } else {
        if (length(unique(df.treated[,1]))==2){ y = as.factor(df.treated[,1]) } else { y = df.treated[,1] }
        lm.treated <- train(y= y,
                       x=data.frame(df.treated[,-1]), method="gbm",
                       trControl = trainControl(method="cv",number=5))
      }

      df.control <- df.control[!duplicated(df.control[,2:ncol(df.control)]),]
      one.level.control <- sapply(df.control, FUN = function(x) length(unique(x))==1)
      df.control <- df.control[,!one.level.control]
      fmla <- as.formula(paste0("y ~ ", paste( colnames(df.control[,2:ncol(df.control)]), collapse= "+")))
      if (out.method == "glm"){
        if (length(unique(df.control[,1]))==2) { y = as.factor(df.control[,1]) } else { y = df.control[,1] }
        if (regularize){
          x <- model.matrix(~., data = data.frame(df.control[,-1]))
          if (is.factor(y)){
            cv.fit <- glmnet::cv.glmnet(x=x, 
                                        y=y, alpha = 0,
                                        family = "binomial")
            lm.control <- glmnet::glmnet(x=x, y=y, alpha = 0, family = "binomial",
                                         lambda = cv.fit$lambda.min)
          } else {
            cv.fit <- glmnet::glmnet(x=x, 
                                     y=y, alpha = 0, nlambda = 50)
            lm.control <- glmnet::glmnet(x=x, y=y, alpha = 0,
                                         lambda = cv.fit$lambda.min)
          }
        } else {
          if (is.factor(y)){
            lm.control <- train(y=y,
                              x=data.frame(df.control[,-1]), method="glm", family="binomial",
                              trControl = trainControl(method="cv",number=5))
          } else {
            lm.control <- train(y=y,
                                x=data.frame(df.control[,-1]), method="glm", 
                                trControl = trainControl(method="cv",number=5))
          }
        }
      } else {
        if (length(unique(df.control[,1]))==2) { y = as.factor(df.control[,1]) } else { y = df.control[,1] }
        lm.control <- train(y= y,
                       x=data.frame(df.control[,-1]), method="gbm",
                       trControl = trainControl(method="cv",number=5))
      }

      colnames(X2) <- paste("X", 1:dim(X2)[2], sep="")
      if (regularize){
        if (length(unique(outcome))==2){
          x <- model.matrix(~., data = X2[,!one.level.treated[-1]])
          y_1.hat <- predict(lm.treated, newx=x, type="response")
          x <- model.matrix(~., data = X2[,!one.level.control[-1]])
          y_0.hat <- predict(lm.control, newx=x, type="response")
        } else {
          x <- model.matrix(~., data = X2[,!one.level.treated[-1]])
          y_1.hat <- predict(lm.treated, newx=x)
          x <- model.matrix(~., data = X2[,!one.level.control[-1]])
          y_0.hat <- predict(lm.control, newx=x)
        }
      } else {
        if (length(unique(outcome))==2){
          y_1.hat <- predict(lm.treated, X2[,!one.level.treated[-1]], type="prob")[,2]
          y_0.hat <- predict(lm.control, X2[,!one.level.control[-1]], type="prob")[,2]
        } else {
          y_1.hat <- as.numeric(predict(lm.treated, X2[,!one.level.treated[-1]]))
          y_0.hat <- as.numeric(predict(lm.control, X2[,!one.level.control[-1]]))
        }
      }
    }
    # EM for linear regression Y ~ X with missing values to predict potential outcomes 
    else {
      # Parameter estimation of potential outcome model under treatment
      yX <- as.matrix(cbind(outcome[which(treat==max(as.numeric(treat)))], X[which(treat==max(as.numeric(treat))),]))
      
      if (is.null(mask)){
        X2 <- yX
      } else {
        X2 <- as.matrix(X[which(treat==max(as.numeric(treat))),])
        colname <- colnames(data.frame(X2))
        length.mask <- dim(mask)[2]
        X2 <- data.frame(cbind(X2, as.matrix(mask[which(treat==max(as.numeric(treat))),])))
        if (length.mask>0){
          colnames(X2) <- c(colname, paste0("R",1:length.mask))
        }
        if (use.interaction & length.mask > 0 & length(col.complete) > 0){
          xnam1 <- paste("X", 1:length.covariates, sep="")
          xnam2 <- paste("X", col.complete, sep="")
          rnam <- paste("R", 1:length.mask, sep="")
          fmla <- as.formula(paste("~ ", paste0(paste(c(xnam1,rnam), collapse= "+"),
                                                "+",
                                                paste(kronecker(xnam2, rnam,function(x,r) paste0(x,"*",r)), collapse= "+"))))
          na.action.default <- getOption("na.action")
          options(na.action = "na.pass")
          X2 <- data.frame(model.matrix(fmla, data = data.frame(X2)))
          X2 <- cbind(outcome[which(treat==max(as.numeric(treat)))], X2)
          options(na.action = na.action.default)
        } else {
          colname <- colnames(X2)
          X2 <- cbind(outcome[which(treat==max(as.numeric(treat)))], X2)
          colnames(X2) <- c("y", colname)
        }
      }
      dim2.X <- dim(X2)[2]
      tmp <- X2
      tmp[is.na(tmp)] <- -1e5
      rm.cols.treated <- c()
      rm.cols.treated <- caret::findLinearCombos(tmp)$remove
      dim2.X <- dim2.X - length(rm.cols.treated)
      if (length(rm.cols.treated)>0){
        mod.t <- prelim.norm(as.matrix(X2[,-rm.cols.treated]))
      } else{
        mod.t <- prelim.norm(as.matrix(X2))
      }
      thetahat <- em.norm(mod.t, showits = F)
      pars <- getparam.norm(mod.t, thetahat)
      sig.inv <- solve(pars$sigma[2:dim2.X,2:dim2.X])
      beta_treated <- c(pars$mu[1] - 
                           pars$sigma[1,2:dim2.X] %*% sig.inv %*% pars$mu[2:dim2.X], 
                           pars$sigma[1,2:dim2.X] %*% sig.inv)
      
      
      # Parameter estimation of potential outcome model under control
      yX <- as.matrix(cbind(outcome[which(treat==min(as.numeric(treat)))], X[which(treat==min(as.numeric(treat))),]))
      
      if (is.null(mask)){
        X2 <- yX
      } else {
        X2 <- as.matrix(X[which(treat==min(as.numeric(treat))),])
        colname <- colnames(data.frame(X2))
        length.mask <- dim(mask)[2]
        X2 <- data.frame(cbind(X2, as.matrix(mask[which(treat==min(as.numeric(treat))),])))
        if (length.mask>0){
          colnames(X2) <- c(colname, paste0("R",1:length.mask))
        }
        if (use.interaction & length.mask > 0 & length(col.complete) > 0){
          xnam1 <- paste("X", 1:length.covariates, sep="")
          xnam2 <- paste("X", col.complete, sep="")
          rnam <- paste("R", 1:length.mask, sep="")
          fmla <- as.formula(paste("~ ", paste0(paste(c(xnam1,rnam), collapse= "+"),
                                                "+",
                                                paste(kronecker(xnam2, rnam,function(x,r) paste0(x,"*",r)), collapse= "+"))))
          na.action.default <- getOption("na.action")
          options(na.action = "na.pass")
          X2 <- data.frame(model.matrix(fmla, data = data.frame(X2)))
          X2 <- cbind(outcome[which(treat==min(as.numeric(treat)))], X2)
          options(na.action = na.action.default)
        } else {
          colname <- colnames(X2)
          X2 <- cbind(outcome[which(treat==min(as.numeric(treat)))], X2)
          colnames(X2) <- c("y", colname)
        }
      }
      dim2.X <- dim(X2)[2]
      tmp <- X2
      tmp[is.na(tmp)] <- -1e5
      rm.cols.control <- c()
      rm.cols.control <- caret::findLinearCombos(tmp)$remove
      dim2.X <- dim2.X - length(rm.cols.control)
      if (length(rm.cols.control)>0){
        mod.t <- prelim.norm(as.matrix(X2[,-rm.cols.control]))
      } else{
        mod.t <- prelim.norm(as.matrix(X2))
      }
      thetahat <- em.norm(mod.t, showits = F)
      pars <- getparam.norm(mod.t, thetahat)
      sig.inv <- solve(pars$sigma[2:dim2.X,2:dim2.X])
      beta_control <- c(pars$mu[1] - 
                          pars$sigma[1,2:dim2.X] %*% sig.inv %*% pars$mu[2:dim2.X], 
                        pars$sigma[1,2:dim2.X] %*% sig.inv)
      
      # Prediction of potential outcomes
      X.imp.ml <- get_imputeEM(as.matrix(X))$Ximp
      
      if (is.null(mask)){
        X2 <- X.imp.ml
      } else {
        colname <- colnames(X.imp.ml)
        length.mask <- dim(mask)[2]
        X2 <- data.frame(cbind(X.imp.ml, as.matrix(mask)))
        if (length.mask>0){
          colnames(X2) <- c(colname, paste0("R",1:length.mask))
        }
        if (use.interaction & length.mask > 0 & length(col.complete) > 0){
          xnam1 <- paste("X", 1:length.covariates, sep="")
          xnam2 <- paste("X", col.complete, sep="")
          rnam <- paste("R", 1:length.mask, sep="")
          fmla <- as.formula(paste("~ ", paste0(paste(c(xnam1,rnam), collapse= "+"),
                                                "+",
                                                paste(kronecker(xnam2, rnam,function(x,r) paste0(x,"*",r)), collapse= "+"))))
          na.action.default <- getOption("na.action")
          options(na.action = "na.pass")
          X2 <- data.frame(model.matrix(fmla, data = data.frame(X2)))
          options(na.action = na.action.default)
        }
      }
      
      X2.treated <- X2
      X2.control <- X2
      if (length(rm.cols.treated)>0){
        X2.treated <- X2[,-rm.cols.treated]
      }
      if (length(rm.cols.control)>0){
        X2.control <- X2[,-rm.cols.control]
      }
      y_1.hat <- as.matrix(X2.treated, ncol=dim(X2.treated)[2])%*%beta_treated[2:(dim(X2.treated)[2]+1)] + beta_treated[1]
      y_0.hat <- as.matrix(X2.control, ncol=dim(X2.control)[2])%*%beta_control[2:(dim(X2.control)[2]+1)] + beta_control[1]
    }
    
    if (!is.null(subset)){
      X <- X[subset,]
      y_1.hat <- y_1.hat[subset]
      y_0.hat <- y_0.hat[subset]
      treat <- treat[subset]
      outcome <- outcome[subset]
      W.hat <- W.hat[subset]
    }
    n.sample <- dim(X)[1]
    delta_i <- y_1.hat -  y_0.hat + treat*(outcome-y_1.hat)*W.hat - (1-treat)*(outcome-y_0.hat)*W.hat
    treatment_effect_dr <- mean(delta_i)
    treatment_effect_dr <- c(treatment_effect_dr, sqrt(sum((delta_i-treatment_effect_dr)^2))/n.sample)
  
  }
  return(data.frame("dr"=treatment_effect_dr[1], "se"=treatment_effect_dr[2]))
}
