#suppressMessages(library("adaptMT"))
#suppressMessages(library("SNPknock"))
#suppressMessages(library("gam"))
#suppressMessages(library("randomForest"))
#suppressMessages(library("expm"))
#suppressMessages(library("rdist"))



#' @noRd
obj_fun <- function(lam, lambda_0, beta, R, n, d = 20){
  p <- length(beta)/2
  r <- ncol(R)
  
  eta <- as.vector(R %*% lam) / d
  max_eta = max(eta)
  w = exp(eta-max_eta) * p / (sum(exp(eta-max_eta))) 
  w <- w * lambda_0
  combined_beta <- abs(beta[1:p]) + abs(beta[(p+1):(2*p)])
  
  # part1 = - 2 * sum(eta)
  part2<- sum( w * combined_beta )
  part3 <- 0.5 * sum(lam^2) / n
  
  return(part2 + part3)
}


#' calculate power and FDR
#' @import knockoff
#' @noRd
calc_power_fdr = function(beta, p, alpha, nonzero) {
  T = abs(beta[1:p])
  T_tilde = abs(beta[(p+1):(2*p)])
  T_max = pmax(T,T_tilde)
  W = T - T_tilde
  
  tau = knockoff.threshold(W,fdr = alpha,offset = 1)
  rej2 = as.numeric(which(W >= tau))
  power_temp = power_cal(rej2, nonzero)
  fdr_tmp <- fdr_cal(rej2, nonzero)
  c(power_temp, fdr_tmp)
}


#' Calculate power
#' @param index_result index for selected covariates
#' @param index_true index for true relevant covariates
#' @return statistical power
#' @export
power_cal <- function(index_result, index_true){
  if(length(index_result) == 0){
    return(0)
  }
  index_common <- intersect(index_result, index_true)
  return(length(index_common)/max(length(index_true), 1))
}

#' Calculate FDR
#' @param index_result index for selected covariates
#' @param index_true index for true relevant covariates
#' @return realized FDR
#' @export
fdr_cal <- function(index_result, index_true){
  if(length(index_result) == 0){
    # return(NA)
    return(0)
  }
  index_common <- intersect(index_result, index_true)
  return(1-length(index_common)/max(length(index_result), 1))
}


################# knockoff_anno

#' Annotation-informed Knockoffs
#' @param X Orignal design matrix
#' @param Xk Knockoff copies of X
#' @param y response vector
#' @param R matrix of standardized annotation matrix
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom stats var sd
#' @return A list contains the following results: Prior_final, final weight phi_j; beta, final beta in adaptive lasso; lambda_s, estimated lambda_l parameter;
#' @export
AnnoKn <- function(X, Xk, y, lams = NULL, attempts = NULL, R){
  # knockoff_anno want to improve the efficiency.
  
  ##### adjust the situation when R is all NAs
  var_para <- apply(R, 2, var)
  
  if(all(is.na(R)) | sum(var_para!=0) == 0){
    print("no informative annotations")
    p = ncol(X)
    mdl = cv.glmnet(cbind(X, Xk), y, alpha=1)
    cvlambda = mdl$lambda.min
    beta = mdl$glmnet.fit$beta[, mdl$lambda == cvlambda]
    
    return(list(lams = NA, # A (set of lambda_0)
                lambda_BEST = cvlambda, # optimal lambda_0
                Prior_final = rep(0, 2*p), # final weight phi_j
                beta = beta, # final beta in adaptive lasso
                cv_error = NA, # cv error for each lamdba_0 in A
                obj_fun_all = NA,
                att_which = NA,
                lambda_s = rep(0, ncol(R))))
  }
  
  R <- R[, var_para!=0]
  R = scale(R)
  #####
  
  X = scale(X)
  Xk = scale(Xk)
  y = (y-mean(y))/sd(y)
  
  
  lambda_s_final = NA
  if(is.null(attempts)){
    attempts = c(0)
  }
  
  n = nrow(X)
  p = ncol(X)
  # R = matrix(scale(R))
  power <- numeric()
  fdr <- numeric()
  
  obj_fun_all <- numeric()
  att_which <- numeric()
  cv_list <- cv_list_mse <- numeric()
  
  nattempts = length(attempts)
  Prior_final = numeric(2 * p) + 1
  CV_BEST = Inf
  lambda_BEST = 0
  
  if(is.null(lams)){
    mdl = cv.glmnet(cbind(X,Xk),y,alpha=1, penalty.factor = Prior_final)
    lambda <- mdl$lambda.min
    lams <- exp(seq(log(lambda)-1, log(lambda)+2, 0.1))
  }
  
  init_lam <- numeric(ncol(R)) + 5
  
  
  for (lam in lams) {
    
    print(paste0("starting lambda0 = : ", lam))
    
    # lambda_0 = lam
    BEST_OBJ = Inf
    
    for (i in 1:nattempts) {
      sol = sol_beta_and_lam(X = X, Xk = Xk, y = y, R = R, lambda_0 = lam, 
                             init_scale = attempts[i], verbose = TRUE, init_lam = init_lam)  
      obj_curr = sol$objs[length(sol$objs)]
      if (obj_curr < BEST_OBJ) {
        att_tmp <- i
        BEST_OBJ = obj_curr
        Prior = sol$Prior
        lambda_s = sol$lambda_s_all[length(sol$lambda_s_all)]
      }
      #print the result: cat("======= attempt ", i, " obj = ", obj_curr, "best_obj = ", BEST_OBJ, "=========\n")
    }
    
    obj_fun_all <- append(obj_fun_all, BEST_OBJ)
    att_which <- append(att_which, att_tmp)
    
    
    mse_cv <- cv_mse(X = cbind(X, Xk), y = y, lambda = lam, nfolds = 10, Prior = Prior)
    cv_list_mse <- append(cv_list_mse, mse_cv)
    
  
    if(mse_cv < CV_BEST){
      lambda_BEST = lam
      CV_BEST = mse_cv
      Prior_final = Prior
      lambda_s_final = lambda_s
    }
    
  }
  
  beta <- as.vector(glmnet(cbind(X,Xk), y, lambda = lambda_BEST * sum(Prior_final) / (2 * p), 
                           penalty.factor = Prior_final, intercept = FALSE, standardize = FALSE)$beta)

  
  list(lams = lams, # A (set of lambda_0)
       lambda_BEST = lambda_BEST, # optimal lambda_0
       Prior_final = Prior_final, # final weight phi_j
       beta = beta, # final beta in adaptive lasso
       cv_error = cv_list_mse, # cv error for each lamdba_0 in A
       obj_fun_all = obj_fun_all,
       att_which = att_which,
       lambda_s =lambda_s_final,
       R = R) # target parameter
  
}




#' @importFrom glmnet glmnet
#' @importFrom stats optim
#' @noRd
sol_beta_and_lam = function(X, Xk, y, R, lambda_0, init_scale = NULL, maxiter = 100, 
                            verbose = TRUE, init_lam = NULL) {
  
  if(is.null(init_scale)){
    init_scale = 0
  }
  
  p = ncol(X)
  n = nrow(X)
  
  Prior = numeric(2 * p)+1
  beta <- as.vector(glmnet(cbind(X,Xk), y, lambda = lambda_0 * init_scale, penalty.factor = Prior)$beta)
  
  beta_old = beta
  lambda_s_old = numeric(ncol(R)) + 100
  
  lambda_s_all <- numeric()
  objs = c()
  
  for(t in 1:maxiter){
    
    res <- optim(
      par   = init_lam,
      fn    = obj_fun,   
      method = "BFGS", 
      lambda_0 = lambda_0,
      beta  = beta,
      R     = R,
      n     = n
    )
    
    lambda_s = res$par
    lambda_s_all <- append(lambda_s_all, list(lambda_s))
    
    w = weight_standardized(lambda_s, R)
    
    Prior = rep(w, 2)
    
    beta <- as.vector(glmnet(cbind(X,Xk), y, lambda = lambda_0 * sum(Prior) / (2 * p), penalty.factor = Prior, intercept = FALSE, standardize = FALSE)$beta)
    
    obj_curr = obj_fun(lambda_s, lambda_0, beta, R, n)
    obj_old = obj_fun(lambda_s_old, lambda_0, beta_old, R, n)
    obj_value = obj_curr - obj_old
    objs = c(objs, obj_curr)
    
    if (verbose) {
      
      cat("iter = ", t, " lambda_s = ", lambda_s, " diff lambda_s = ", abs(lambda_s - lambda_s_old), " obj diff = ", obj_value, " sum(Prior) = ", sum(Prior), "\n")
      if(mean(abs(lambda_s - lambda_s_old)) < 0.001 || abs(obj_value) < 1e-5){
        print(paste0("convergence achieved at ", t, "-th iteration"))
        break
      }
    }
    
    lambda_s_old = lambda_s
    
    beta_old = beta
  }
  list(beta = beta, objs = objs, lambda_s_all = lambda_s_all, Prior = Prior)
}


############################## objective function used in the AnnoKn method: ##############################

#' @importFrom glmnet glmnet
#' @importFrom stats predict
#' @noRd
cv_mse <- function(X, y, lambda, nfolds = 5, family = "gaussian", Prior = Prior) {
  # set.seed(seed)
  n <- nrow(X)
  folds <- sample(rep(1:nfolds, length.out = n))
  mse_list <- numeric(nfolds)
  
  for (k in 1:nfolds) {
    test_idx <- which(folds == k)
    train_idx <- setdiff(1:n, test_idx)
    
    X_train <- X[train_idx, , drop = FALSE]
    y_train <- y[train_idx]
    X_test  <- X[test_idx, , drop = FALSE]
    y_test  <- y[test_idx]
    
    # Fit model using glmnet with fixed lambda
    fit <- glmnet(X_train, y_train, lambda = lambda, family = family, penalty.factor = Prior)
    
    # Predict on the test fold
    y_pred <- predict(fit, newx = X_test)
    
    # Compute MSE for this fold
    mse_list[k] <- mean((y_test - y_pred)^2)
  }
  
  # Return the average CV MSE
  return(mean(mse_list))
}




################################# knockoff_simple #################################

#' Annotation-informed Knockoffs
#' @param X Orignal design matrix
#' @param Xk Knockoff copies of X
#' @param y response vector
#' @param R matrix of standardized annotation matrix
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom stats var sd
#' @return A list contains the following results: Prior, weight phi_j; beta, beta in adaptive lasso; lambda_s, estimated lambda_l parameter;
#' @export
AnnoKn_lite = function(X, Xk, y, R) {
  
  var_para <- apply(R, 2, var)
  
  if(all(is.na(R)) | sum(var_para!=0) == 0){
    print("no informative annotations")
    p = ncol(X)
    mdl = cv.glmnet(cbind(X, Xk), y, alpha=1)
    cvlambda = mdl$lambda.min
    beta = mdl$glmnet.fit$beta[, mdl$lambda == cvlambda]
    
    return(list(beta = beta, 
                objs = NA, 
                lambda_s = rep(0, ncol(R)),
                Prior = rep(0, 2*p)))
  }
  
  R <- R[, var_para!=0]
  R = scale(R)
  #####
  
  
  
  X = scale(X)
  Xk = scale(Xk)
  y = (y-mean(y))/sd(y)
  
  n = nrow(X)
  p = ncol(X)
  Prior = numeric(2 * p)+1
  mdl = cv.glmnet(cbind(X, Xk),y,alpha=1, penalty.factor = Prior)
  lambda_0 <- mdl$lambda.min
  
  init_lam <- numeric(ncol(R)) + 5
  result = sol_beta_and_lam_simple(X, Xk, y, lambda_0 = lambda_0, 
                                   R = R, init_lam = init_lam)
  return(result)
}





#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom stats optim
#' @noRd
sol_beta_and_lam_simple = function(X, Xk, y, R, lambda_0, init_scale = NULL, maxiter = 100, 
                            verbose = TRUE, init_lam = NULL) {
  # calculate beta if lambda_0 is given, and then update lambda_0 
  if(is.null(init_scale)){
    init_scale = 0
  }
  
  p = ncol(X)
  n = nrow(X)
  
  Prior = numeric(2 * p)+1
  beta <- as.vector(glmnet(cbind(X,Xk), y, lambda = lambda_0 * init_scale, penalty.factor = Prior)$beta)
  
  beta_old = beta
  lambda_s_old = numeric(ncol(R)) + 100
  
  lambda_s_all <- numeric()
  objs = c()
  
  for(t in 1:maxiter){
    
    res <- optim(
      par   = init_lam,
      fn    = obj_fun,   
      method = "BFGS", 
      lambda_0 = lambda_0,
      beta  = beta,
      R     = R,
      n     = n
    )
    
    lambda_s = res$par
    lambda_s_all <- append(lambda_s_all, lambda_s)
    
    w = weight_standardized(lambda_s, R)
    
    Prior = rep(w, 2)
    
    beta <- as.vector(glmnet(cbind(X,Xk), y, lambda = lambda_0 * sum(Prior) / (2 * p), penalty.factor = Prior, intercept = FALSE, standardize = FALSE)$beta)
    
    obj_curr = obj_fun(lambda_s, lambda_0, beta, R, n)
    obj_old = obj_fun(lambda_s_old, lambda_0, beta_old, R, n)
    obj_value = obj_curr - obj_old
    objs = c(objs, obj_curr)
    
    if (verbose) {
      
      cat("iter = ", t, " lambda_s = ", lambda_s, " diff lambda_s = ", abs(lambda_s - lambda_s_old), " obj diff = ", obj_value, " sum(Prior) = ", sum(Prior), "\n")
    }
    
    if(mean(abs(lambda_s - lambda_s_old)) < 0.001 || abs(obj_value) < 1e-4){
      print(paste0("convergence achieved at ", t, "-th iteration"))
      break
    }
    
    lambda_s_old = lambda_s
    
    beta_old = beta
  }
  
  mdl = cv.glmnet(cbind(X,Xk),y,alpha=1, penalty.factor = Prior)
  lambda <- mdl$lambda.min
  beta <- as.vector(glmnet(cbind(X,Xk), y, lambda = lambda, penalty.factor = Prior)$beta)
  
  list(beta = beta, objs = objs, lambda_s = lambda_s, Prior = Prior)
}



# sourced from ghostknockoff
#' @noRd
create_gaussian_s <- function(X,Sigma,s,chol_V){
  n = nrow(X)
  p = ncol(X)
  D = diag(s)
  E = matrix(rnorm(p*n), n, p)
  P = diag(p) - solve(Sigma, D)
  V = 2*D - D%*%solve(Sigma, D)
  return(X%*%P + E%*%chol_V)
}




