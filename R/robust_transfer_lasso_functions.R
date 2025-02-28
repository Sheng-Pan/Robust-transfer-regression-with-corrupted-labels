

#   Calculate the **Signal-to-Noise Ratio (SNR)** between the true signal and its estimate.
# Formula explanation:
#   
#   SNR=10log 10(Power of true signal/Power of error)
# Handles division by zero by returning infinity (Inf) if power or error is zero.
ser <- function(signal_true, signal_estimated) {
  error <- sum((signal_true - signal_estimated)^2)
  power <- sum(signal_true^2)
  if (power != 0 && error != 0) {
    return(10 * log10(power / error))
  } else {
    return(Inf)
  }
}

###############################
#generate target dataset
###############################
#generate dataset
generate_target_data <- function(i,n,p,s,sigma,k){
  set.seed(i)
  # Create Gaussian measurement matrix Phi normalized by sqrt(n)
  Phi <- matrix(rnorm(n * p), nrow = n, ncol = p) / sqrt(n)
  # Generate sparse signal with s non-zero entries
  signal_true <- rep(0, p)
  nonzero_indices0 <- sample(p, s)
  signal_true[nonzero_indices0] <- sign(rnorm(s))
  # Clean measurement without corruption
  measurements_clean <- Phi %*% signal_true
  # Add Cauchy noise (alternative to α-stable noise)
  noise <- rnorm(n,0,sigma)#
  # Create corruption vector with m/2 zeros and uniform values in others
  vector <- rep(NA, n)
  zero_indices <- sample(n, size = floor((1-k)*n), replace = FALSE)
  vector[zero_indices] <- 0
  non_zero_indices <- setdiff(seq_len(n), zero_indices)
  vector[non_zero_indices] <-  runif(length(non_zero_indices), min = 0.5, max = 1)
  # Return dataset components
  return(list(X = Phi, 
              y =  measurements_clean + noise + vector,
              signal_true = signal_true,
              non_zero_indices = nonzero_indices0,
              corruption =  vector))
}

###############################
#generate source data
###############################

create_auxiliary_data <- function(j,e,n,p) {
  set.seed(j)  
  # Generate auxiliary dataset with modified signal_true
  signal_true1 <- c(signal_true, rep(0,p-length(signal_true)))
  a = sample(e,1)
  signal_true1[sample(nonzero_indices, a)] <- 0
  signal_true1[sample(setdiff(1:n,nonzero_indices), a)] <- sign(rnorm(a))
  tau = rbinom(1, size = 1, prob = 1/L)
  signal_true1[sample(which(signal_true1==0),8)] <- tau*sign(rnorm(8))
  # Create Gaussian measurement matrix Phi_au normalized by sqrt(n)
  Phi_au <- matrix(rnorm(n * p), nrow = n, ncol = p) / sqrt(n)
  # Clean measurement without corruption
  measurements_clean <- Phi_au %*% signal_true1
  # Add small Gaussian noise
  noise <- rnorm(n,0,0.01)#
  return(list(X = Phi_au, 
              y =  measurements_clean + noise,
              signal = signal_true1,
              error = a))
}

###############################
#Rlasso
###############################
#Implements the ​soft thresholding operator​ used in ADMM for L1 regularization.
soft_threshold <- function(a, kappa) {
  return(pmax(0, a - kappa) - pmax(0, -a - kappa))
}
Rlasso <- function(y0, Phi, tau_x, tau_r,method, rho=1, max_iter = 20, tol = 1e-6) {
  if(method == "ADMM"){
    
    n <- length(y0)
    p <- ncol(Phi)
    
    x <- rep(0, p)
    r <- rep(0, n)
    z <- rep(0, p)
    v <- rep(0, n)
    u <- rep(0, p)
    w <- rep(0, n)
    
    
    Phi_T <- t(Phi)
    I_p <- diag(p)
    I_n <- diag(n)
    
    # ADMM 
    for (k in 1:max_iter) {
      # 
      A <- rbind(cbind(Phi_T %*% Phi + rho * I_p, Phi_T),
                 cbind(Phi, I_n + rho * I_n))
      b <- c(Phi_T %*% y0 + rho * (z - u),
             y0 + rho * (v - w))
      xr_update <- solve(A, b)
      x_next <- xr_update[1:p]
      r_next <- xr_update[(p+1):(p+n)]
      
      # 
      z_next <- soft_threshold(x_next + u, tau_x / rho)
      
      #
      v_next <- soft_threshold(r_next + w, tau_r / rho)
      
      # 
      u_next <- u + x_next - z_next
      w_next <- w + r_next - v_next
      
      # 
      primal_residual <- sqrt(sum((x_next - z_next)^2) + sum((r_next - v_next)^2))
      dual_residual <- rho * sqrt(sum((z_next - z)^2) + sum((v_next - v)^2))
      
      if (primal_residual < tol && dual_residual < tol) {
        cat("Converged at iteration", k, "\n")
        break
      }
      
      # 
      x <- x_next
      r <- r_next
      z <- z_next
      v <- v_next
      u <- u_next
      w <- w_next
    }
    
    return(list(x=x,r=r))
  }else if(method == 'cd'){ n <- length(y0)
  p <- ncol(Phi)
  
  
  x <- rep(0, p)
  r <- rep(0, n)
  
  for (iter in 1:max_iter) {
    x_old <- x
    r_old <- r
    
    
    for (j in 1:p) {
      
      residuals <- y0 - Phi %*% x - r
      aj <- crossprod(Phi[, j], Phi[, j])
      cj <- crossprod(Phi[, j], residuals + Phi[, j] * x[j])
      
      
      x[j] <- soft_threshold(cj / aj, tau_x / aj)
      
    }
    
    
    for (i in 1:n) {
      
      residuals <- y0 - Phi %*% x - r
      ai <- 1
      ci <- residuals[i] + r[i]
      
      
      r[i] <- soft_threshold(ci / ai, tau_r / ai)
    }
    print(sum(abs(result0$x-x)))
    print(iter )
    # 
    primal_residual <- sqrt(sum((x - x_old)^2) + sum((r - r_old)^2))
    if (primal_residual < tol) {
      cat("Converged at iteration", iter, "\n")
      break
    }
  }
  
  return(list(x=x, r=r))}else if(method == 'PDIPM'){
    n = dim(x0)[1]
    p = dim(x0)[2]
    x <- Variable(p)
    r <- Variable(n)
    
    #
    
    objective <- (1/2) * sum((y0 - x0 %*% x - r)^2) + tau_x * sum(abs(x))+ tau_r* sum(abs(r))
    
    # 
    problem <- Problem(Minimize(objective))
    result <- solve(problem, solver = "ECOS")
    x = result$getValue(x)
    r = result$getValue(r)
    return(list(x=x, r=r))
  }else {
    stop("Unsupported method: ", method)
  }
  
}
###############################
#distributed_lasso
###############################
distributed_lasso <- function(Xs, ys, Ahat, max_iter = 40, rho = 0.1) {
  
  if(length(Ahat)==1){
    cv.lasso <- cv.glmnet(Xs[[1]], ys[[1]], alpha=1)  # alpha=1 for Lasso
    fit_w <- glmnet(Xs[[1]], ys[[1]], alpha = 1, lambda = cv.lasso$lambda.min, standardize = FALSE)
    lasso <- as.numeric(coef(fit_w)[-1])
    return(lasso)
  }else{
    beta_hat <- rep(0, p)  
    z <- rep(0, p)         
    u <- rep(0, p)         
    X_test = Xs[[test_data]]
    y_test = ys[[test_data]]
    for(t0 in 1:1){
      for (tt in 1:max_iter) {
        
        XtX <- crossprod(X_test)
        Xty <- crossprod(X_test, y_test)
        
        
        gradient_term <- 0
        for (j in Ahat[-which(Ahat == test_data)]) {
          Xj <- Xs[[j]]
          yj <- ys[[j]]
          gradient_term <- gradient_term + t(Xj %*% beta_hat - yj) %*% Xj
        }
        g0 <- t(X_test %*% beta_hat - y_test) %*% X_test
        a0 <- (g0 + gradient_term) / length(Ahat) - g0
        a0 <- a0/n
        
        XtX <- crossprod(X_test)
        Xty <- crossprod(X_test, y_test - (a0 %*% beta_hat)[1])
        inv_part <- solve(XtX + rho * diag(p))
        beta_hat <- as.vector(inv_part %*% (Xty + rho * (z - u)))
        
        # 
        z_old <- z
        beta_hat_rho <- beta_hat + u
        lambda_t <- 0.01 * sqrt(log(p) / (length(Ahat) * n)) + sqrt(log(p) / n) * (0.01 * 12 * sqrt(log(p) / n))^t0
        z <- soft_threshold(beta_hat_rho, lambda_t / rho)
        
        # 
        u <- u + beta_hat - z

        # 
        if (sqrt(sum((z - z_old)^2)) < 1e-4) break
      }
    }
    
    
    return(z)
  }
  
}###############################
#adaptive cv
###############################
###cv
cv_RTL2 <- function(Xs,ys,X0,y0,folds,ctilde,k0,r,r2,test_data,w_hat_A0,lambda_beta_grid, lambda_e_grid){
  if(r2 > r){
    X_test=Xs[[test_data]] 
    y_test = ys[[test_data]]
    opf2 <- function(z){
      tau_x = z[1]
      tau_r = z[2]
      result1 =  RTL_admm_solver_cpp(y0, X0, w_hat_A0, tau_x, tau_r, rho=0.1, max_iter = 40, tol = 1e-6)
      ab0 = mean((y_test  - X_test %*% result1$beta_hat)^2)
      return(ab0 + 10*(sum(abs(result1$beta_hat-w_hat_A0))<ctilde) )
    }
  }else{
    opf2 <- function(z){
      tau_x = z[1]
      tau_r = z[2]
      
      scores=NULL
      for (i in 1:k0) {
        set = which(folds == i)
        X_test=X0[ set, ]
        y_test = y0[ set ]
        X_train=X0[ -set, ]
        y_train = y0[-set ]
        result = RTL_admm_solver_cpp(y_train, X_train, w_hat_A0, tau_x, tau_r, rho=0.1, max_iter = 40, tol = 1e-6)
        score = mean((y_test  - c(X_test %*% result$beta_hat))^2)
        scores <- c(scores, score)
      }
      return( mean(scores))
    }
  }
  
  
  best_score <- Inf
  best_params <- list(lambda_beta = NA, lambda_e = NA)
  
  for (lambda_beta in lambda_beta_grid) {
    for (lambda_e in lambda_e_grid) {
      score <- tryCatch({opf2(c(lambda_beta ,lambda_e))},
                        error = function(e) {
                          1000
                        })
      #print(c(lambda_beta ,lambda_e,score))
      if (score < best_score) {
        best_score <- score
        best_params <- list(lambda_beta = lambda_beta, lambda_e = lambda_e)
      }
    }
  }
  return(best_params)
}

apply2 = function(x){
  if(length(dim(x))==0){x}else{apply(x,1,mean)}
}

robust_trans_lasso <- function(X0, y0,Xs,ys,method,Ahat,w_hat_A0,test_data, h,lambda_beta,  lambda_e) {
  if(method == 'ADMM'){
    # ADMM-based robust Trans-Lasso
    n = dim(X0)[1]
    p = dim(X0)[2]
    
    if(length(Ahat)==0){
      return(rep(0,p))
    }else{
      result = RTL_admm_solver_cpp(y0, X0, w_hat_A0, lambda_beta, lambda_e, rho=0.1, max_iter = 100, tol = 1e-6)
      return(result)
    }
  } else if(method == 'PDIPM'){
    # PDIPM-based robust Trans-Lasso
    x <- Variable(p)
    r <- Variable(n)
    objective <- (1/2) * sum((y0 -  X0 %*% w_hat_A0 -  X0 %*% x - r)^2) + lambda_beta * sum(abs(x))+ lambda_e* sum(abs(r))
    problem <- Problem(Minimize(objective))
    result <- solve(problem, solver = "ECOS")
    return(list(beta_hat=result$getValue(x),r=result$getValue(r)))
  } else {
    stop("Unsupported method")
  }
}


