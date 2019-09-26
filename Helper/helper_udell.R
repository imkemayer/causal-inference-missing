

make_V <- function(r, p){
  matrix(rnorm(r*p), nrow = p, ncol = r)
}

design_matrix <- function(V, n, r, p){
  # this function constructs the design list with U, V, X components 
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

compute_folds <- function(Xm, nfolds = 3){
  # create cross-validation folds for gaussian matrix factorization 
  
  n = nrow(Xm); p = ncol(Xm)
  nfold_train = createFolds(1:n, nfolds, returnTrain = T)
  pfold_train = createFolds(1:p, nfolds, returnTrain = T)
  nfold_test = lapply(nfold_train, function(x) setdiff(1:n, x))
  pfold_test = lapply(pfold_train, function(x) setdiff(1:p, x))
  list(nfold_train = nfold_train, pfold_train = pfold_train, 
       nfold_test = nfold_test, pfold_test = pfold_test)
}

cross_valid <- function(X, r, warm, folds, nfolds = 3){
  # compute the cross validation error for gaussian matrix factorization with different ranks 
  
  assert_that(length(folds$nfold_train) == nfolds)
  
  nfold_train = folds$nfold_train
  pfold_train = folds$pfold_train
  nfold_test = folds$nfold_test
  pfold_test = folds$pfold_test
  
  error_folds = numeric(nfolds)
  fit_folds = list()
  
  for (f in 1:nfolds){
    temp_data = X
    temp_data[nfold_test[[f]], pfold_test[[f]]] = NA
    fit = softImpute(temp_data, rank.max = r, type = "als", maxit = 1000, warm.start = warm)
    pred = impute(fit, i = rep(nfold_test[[f]], length(pfold_test[[f]])), 
                  j = rep(pfold_test[[f]], each = length(nfold_test[[f]])))
    assert_that(length(c(X[nfold_test[[f]], pfold_test[[f]]])) == length(pred))
    error = mean((c(X[nfold_test[[f]], pfold_test[[f]]]) - pred)^2, na.rm = T)
    error_folds[f] = error
    fit_folds[[f]] = fit
  }
  list(error = mean(error_folds), fit = fit_folds[[which.min(error_folds)]])
}

recover_pca_gaussian_cv <- function(X, r_seq, nfolds = 3){
  # Gaussian matrix factorization on the noisy proxy matrix design$Xm 
  #     with rank chosen by cross validation from r_seq 
  #     the matrix factorization is carried out by the softImpute package  
  
  cv_error = numeric(length(r_seq))
  warm_list = list()
  folds = compute_folds(X)
  
  for (r in 1:length(r_seq)){
    if (r == 1){
      temp = cross_valid(X, r_seq[r], warm = NULL, folds = folds, nfolds = nfolds)
      cv_error[r] = temp$error
      warm_list[[r]] = temp$fit
    } else{
      temp = cross_valid(X, r_seq[r], warm = warm_list[[r-1]], folds = folds, nfolds = nfolds)
      cv_error[r] = temp$error
      warm_list[[r]] = temp$fit
    }
  }
  
  best_r = r_seq[which.min(cv_error)]
  warm = warm_list[[which.min(cv_error)]]
  result = softImpute(X, rank.max = best_r, type = "als", maxit = 1000, warm.start = warm)
  list(Uhat = result$u, best_r = best_r, Xhat = result$u%*%diag(result$d)%*%t(result$v))
}

recover_pca_gaussian <- function(X, r){
  # Gaussian matrix factorization on the noisy proxy matrix design$Xm 
  #     with given rank
  #     the matrix factorization is carried out by the softImpute package  
  
  result = softImpute(X, rank.max = r, type = "als", maxit = 1000)
  list(Uhat = result$u, Xhat = result$u%*%diag(result$d)%*%t(result$v))
}