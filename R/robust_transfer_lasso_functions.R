
msefun = function(x,y){
  sum(x-y)^2
}
l1fun = function(x,y){
  sum(abs((x-y)))
}
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
library(MASS)
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

generate_target_data2 <- function(i,scenario ,n,p,s,rr,sigma,k,distr){
  set.seed(i)
  if(scenario == 3){
    sigma0 <- diag(p) # 
    mu <- rep(0, p)
    # 
    for (i in 1:(p-1)) {
      sigma0[i, i+1] <- rr # 
      sigma0[i+1, i] <- rr # 
    }
    Phi <- mvrnorm(n , mu = mu, Sigma = sigma0)#/ sqrt(n)
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
    
    if(distr == "unif"){
      er = runif(length(non_zero_indices), min = 5, max = 8)
      # er = rlnorm(length(non_zero_indices), meanlog = 0, sdlog = 1)
      vector[non_zero_indices] <-  er
      noise1 =  vector + noise
    }else if(distr == "alpha-stable"){
      er = abs(rstable(length(non_zero_indices), alpha = 1, beta = 0, gamma = 1, delta = 0))
      #er[which(er>10)] = 10
      vector[non_zero_indices] <-  er
      noise1 =  vector + noise
      #noise1 =  rstable(n, alpha = 1.5-k, beta = 0, gamma = 1, delta = 0)
    }else if(distr == "xb"){
  
      signal_er <- rep(0, p)
      id1 <- sample(p, 2*s)
      signal_er[id1] <- sign(rnorm(s))
      ma = matrix(rnorm(n * p), nrow = n, ncol = p)
      vector[non_zero_indices] <-  4*abs((ma%*%signal_er)[non_zero_indices] )
      noise1 =  vector + noise
    }
    
    # Return dataset components
    return(list(X = Phi, 
                y =  measurements_clean + noise1,
                signal_true = signal_true,
                non_zero_indices = nonzero_indices0,
                corruption =  vector))
  }else if(scenario == 1){
    Phi <- matrix(rnorm(n * p), nrow = n, ncol = p)# / sqrt(n)

    signal_true <- rep(0, p)
    nonzero_indices0 <- sample(p, s)
    signal_true[nonzero_indices0] <- sign(rnorm(s))
    # 
    measurements_clean <- Phi %*% signal_true
    
    # 
    noise <- rnorm(n,0,sigma)#
    
    #
    vector <- rep(NA, n)
    
    # 
    zero_indices <- sample(n, size = floor((1-k)*n), replace = FALSE)
    vector[zero_indices] <- 0
    
    # 
    non_zero_indices <- setdiff(seq_len(n), zero_indices)
    vector[non_zero_indices] <-  runif(length(non_zero_indices), min = 0.5, max = 1)#rcauchy(length(non_zero_indices),location = 0.2, scale = sigma)
    #
    error <- vector
    
    return(list(X = Phi, 
                y =  measurements_clean + noise + error,
                signal_true = signal_true,
                non_zero_indices = nonzero_indices0,
                corruption =  error)
    )
  }else if(scenario == 4){
 
    # 
    Sigma <- matrix(rr, nrow = p, ncol = p)
    diag(Sigma) <- 1  # 对角线为1
    
    # 
    gauss_cop <- normalCopula(P2p(Sigma), dim = p, dispstr = "un")
    #
    Phi  <- qnorm(rCopula(n, gauss_cop) )
    
    # 
    signal_true <- rep(0, p)
    nonzero_indices0 <- sample(p, s)
    signal_true[nonzero_indices0] <- sign(rnorm(s))
    # 
    measurements_clean <- Phi %*% signal_true
    
    # 
    noise <- rnorm(n,0,sigma)#
    
    # 
    vector <- rep(NA, n)
    
    # 
    zero_indices <- sample(n, size = floor((1-k)*n), replace = FALSE)
    vector[zero_indices] <- 0
    
    # 
    non_zero_indices <- setdiff(seq_len(n), zero_indices)
    vector[non_zero_indices] <-  runif(length(non_zero_indices), min = 5, max = 15)#rcauchy(length(non_zero_indices),location = 0.2, scale = sigma)
    #
    error <- vector
    
    return(list(X = Phi, 
                y =  measurements_clean + noise + error,
                signal_true = signal_true,
                non_zero_indices = nonzero_indices0,
                corruption =  error)
    )
  }

  # Create Gaussian measurement matrix Phi normalized by sqrt(n)
  
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
create_auxiliary_data2 <- function(j,scenario ,rr,sigma,e,n,p) {
  set.seed(j)  
  if(scenario ==3){
    # Generate auxiliary dataset with modified signal_true
    signal_true1 <- c(signal_true, rep(0,p-length(signal_true)))
    a = sample(e,1)
    signal_true1[sample(nonzero_indices, a)] <- 0
    signal_true1[sample(setdiff(1:n,nonzero_indices), a)] <- sign(rnorm(a))
    tau = rbinom(1, size = 1, prob = 1/L)
    signal_true1[sample(which(signal_true1==0),8)] <- tau*sign(rnorm(8))
    # Create Gaussian measurement matrix Phi_au normalized by sqrt(n)
    sigma0 <-  diag(p) # 
    mu <- rep(0, p)
    #
    for (i in 1:(p-1)) {
      sigma0[i, i+1] <- rr #
      sigma0[i+1, i] <- rr # 
    }
    Phi_au <- mvrnorm(n , mu = mu, Sigma = sigma0)#/ sqrt(n)
    # Clean measurement without corruption
    measurements_clean <- Phi_au %*% signal_true1
    # Add small Gaussian noise
    noise <- rnorm(n,0,sigma)#
    return(list(X = Phi_au, 
                y =  measurements_clean + noise,
                signal = signal_true1,
                error = a))
  }else if(scenario == 1){
    signal_true1 <- c(signal_true,rep(0,p-length(signal_true)))
    a = sample(e,1)
    signal_true1[sample(nonzero_indices, a)] <- 0
    signal_true1[sample(setdiff(1:n,nonzero_indices), a)] <- sign(rnorm(a))
    tau = rbinom(1, size = 1, prob = 1/L)
    signal_true1[sample(which(signal_true1==0),8)] <- tau*sign(rnorm(8))
    # 
    #Phi_au <- matrix(rnorm(m * n), nrow = m, ncol = n) / sqrt(m)
    Phi_au <- matrix(rnorm(n * p), nrow = n, ncol = p)# / sqrt(n)
    #
    measurements_clean <- Phi_au %*% signal_true1
    # 
    noise <- rnorm(n,0,sigma)#
    return(list(X = Phi_au, 
                y =  measurements_clean + noise,
                signal = signal_true1,
                error = a))
  }else if(scenario == 4){
    signal_true1 <- c(signal_true,rep(0,p-length(signal_true)))
    a = sample(e,1)
    signal_true1[sample(nonzero_indices, a)] <- 0
    signal_true1[sample(setdiff(1:n,nonzero_indices), a)] <- sign(rnorm(a))
    tau = rbinom(1, size = 1, prob = 1/L)
    signal_true1[sample(which(signal_true1==0),8)] <- tau*sign(rnorm(8))
    # 
    #Phi_au <- matrix(rnorm(m * n), nrow = m, ncol = n) / sqrt(m)
    
    Sigma <- matrix(rr, nrow = p, ncol = p)
    diag(Sigma) <- 1  # 
    
    # 2. 
    gauss_cop <- normalCopula(P2p(Sigma), dim = p, dispstr = "un")
    # 4. 
   

    Phi_au <-  qnorm(rCopula(n, gauss_cop) )
    

    # 
    measurements_clean <- Phi_au %*% signal_true1
    # 
    noise <- rnorm(n,0,sigma)#
    return(list(X = Phi_au, 
                y =  measurements_clean + noise,
                signal = signal_true1,
                error = a))
  }
  
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
#
proximal_gradient_descent <- function(X, y, a0, xint, lambda, max_iter = 40, step_size = 0.01, tol = 1e-6) {
  n <- nrow(X)
  p <- ncol(X)
  beta <- xint
  
  for (iter in 1:max_iter) {
    beta_old <- beta
    
    # 
    residual <- y - X %*% beta
    gradient <- (-2/n) * t(X) %*% residual + a0
    
    # 
    step_ok <- FALSE
    temp_step <- step_size
    for (backtrack in 1:100) {
      beta_temp <- beta - temp_step * gradient
      beta_new <- sign(beta_temp) * pmax(abs(beta_temp) - lambda * temp_step, 0)
      
      # 
      old_loss <- mean((y - X %*% beta)^2) + sum(a0 * beta) + lambda * sum(abs(beta))
      new_loss <- mean((y - X %*% beta_new)^2) + sum(a0 * beta_new) + lambda * sum(abs(beta_new))
      
      if (new_loss <= old_loss) {
        step_ok <- TRUE
        break
      }
      #adaptive step size
      temp_step <- temp_step * 0.1  # 
    }
    
    if (!step_ok) {
      break
    }
    
    beta <- beta_new

    if (max(abs(beta - beta_old)) < tol) break
  }
  return(beta)
}
distributed_lasso2 <- function(Xs, ys,test_data, Ahat,nv, max_iter = 40, rho = 0.1) {
  
  if(length(Ahat)==1){
    cv.lasso <- cv.glmnet(Xs[[test_data]], ys[[test_data]], alpha=1)  # alpha=1 for Lasso
    fit_w <- glmnet(Xs[[test_data]], ys[[test_data]], alpha = 1, lambda = cv.lasso$lambda.min, standardize = FALSE)
    lasso <- as.numeric(coef(fit_w)[-1])
    return(lasso)
  }else{
    X_test = Xs[[test_data]]
    y_test = ys[[test_data]]
    X=X_test
    y=y_test
    cv.lasso <- cv.glmnet(X, y, alpha=1)  # alpha=1 for Lasso
    fit_w <- glmnet(X, y, alpha = 1,lambda = cv.lasso$lambda.min)
    beta_hat <- as.numeric(coef(fit_w)[-1])
    z <- rep(0, p)         
    u <- rep(0, p)         

    z <- beta_hat
    #print(ser(s0,beta_hat))
    for(t0 in 1:1){
      for (tt in 1:max_iter) {
        
        gradient_term <- 0
        for (j in Ahat[-which(Ahat == test_data)]) {
          Xj <- Xs[[j]]
          yj <- ys[[j]]
          gradient_term <- gradient_term + t(Xj %*% z - yj) %*% Xj
        }
        g0 <- t(X_test %*% z - y_test) %*% X_test
        a0 <- (g0 + gradient_term) / length(Ahat) - g0
        a0 <- as.vector(a0/n)
        
        #     #adaptive step size
        if(max(abs(a0)) > 2) {
          warning(paste("break：max(abs(a0)) =", max(abs(a0)), "> 2"))
          break
        }
        
        lambda_t <- 1* sqrt(log(p) / (length(Ahat) * n)) + sqrt(log(p) / n) * (0.01 * 12 * sqrt(log(p) / n))^t0
        
        z = proximal_gradient_descent(X_test, y_test, a0, z, lambda_t, max_iter = 40, step_size = 0.01, tol = 1e-6)

      }
    }
    if(max(abs(a0)) > 2){
      z <- beta_hat
      for(t0 in 1:1){
        for (tt in 1:max_iter) {
          
          gradient_term <- 0
          for (j in Ahat[-which(Ahat == test_data)]) {
            Xj <- Xs[[j]]
            yj <- ys[[j]]
            gradient_term <- gradient_term + t(Xj %*% z - yj) %*% Xj
          }
          g0 <- t(X_test %*% z - y_test) %*% X_test
          a0 <- (g0 + gradient_term) / length(Ahat) - g0
          a0 <- as.vector(a0/n)
          
          # 
          if(max(abs(a0)) > 2) {
            warning(paste("break：max(abs(a0)) =", max(abs(a0)), "> 2"))
            break
          }
          
          lambda_t <- 1* sqrt(log(p) / (length(Ahat) * n)) + sqrt(log(p) / n) * (0.01 * 12 * sqrt(log(p) / n))^t0
          
          z = proximal_gradient_descent(X_test, y_test, a0, z, lambda_t, max_iter = 400, step_size = 0.001, tol = 1e-6)
          
          print(ser(s0,z))
          
          # if (sqrt(sum((z - z_old)^2)) < 1e-4) break
        }
      }
    }
    
    
    return(z)
  }
  
}
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
  
}
distributed_lasso_fun <- function(Xs, ys,method,test_data, Ahat,nv, max_iter = 40, rho = 0.1)
{
  if(method == "ADMM"){
    distributed_lasso(Xs, ys, Ahat, max_iter = 40, rho = 0.1) 
  }else if(method == "proximal"){
    distributed_lasso2(Xs, ys,test_data, Ahat,nv, max_iter = 40, rho = 0.1) 
  }
}
###############################
#adaptive cv
###############################
###cv
cv_RTL2 <- function(Xs,ys,X0,y0,Ahat,folds,method,ctilde,k0,r,r2,test_data,w_hat_A0,lambda_beta_grid, lambda_e_grid){
  if(r2 > r){
    X_test=Xs[[test_data]] 
    y_test = ys[[test_data]]
    opf2 <- function(z){
      tau_x = z[1]
      tau_r = z[2]

      result1 =  robust_trans_lasso(X0, y0,Ahat,Xs,ys,method,w_hat_A0,test_data, 30,tau_x,  tau_r)
      ab0 = mean((y_test  - X_test %*% result1$beta_hat)^2)
      print(ser(signal_true , result1$beta_hat))
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
      print(c(lambda_beta ,lambda_e,score))
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

robust_trans_lasso <- function(X0, y0,Ahat,Xs,ys,method,w_hat_A0,test_data, h,lambda_beta,  lambda_e) {
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
    y1 = y0  - result$getValue(r)
    cv.lasso <- cv.glmnet(X0, y1, alpha=1)  # alpha=1 for Lasso
    fit_w <- glmnet(X0, y1, alpha = 1,lambda = cv.lasso$lambda.min)
    #beta_hat=result$getValue(x)+ w_hat_A0
    beta_hat = as.numeric(coef(fit_w)[-1])
   # print(ser(signal_true , beta_hat))
    return(list(beta_hat= beta_hat,r=result$getValue(r)))
  } else {
    stop("Unsupported method")
  }
}
delta_estimation_fun <- function(Xs, ys,methed, tau_x, tau_r, w_hat_A,  X0,   y0,  L,  p,  n){
  if(method == "ADMM" ){
    delta_estimation(Xs, ys, tau_r, tau_x,w_hat_A,  X0,   y0,  L,  p,  n, 1, 40)
  }else{
    delta_estimationR(Xs, ys, tau_x, tau_r,w_hat_A,  X0,   y0,  L,  p,  n)
  }
}
delta_estimationR <- function(Xs, ys, tau_x, tau_r,w_hat_A,  X0,   y0,  L,  p,  n) {
  hhat = shat = h0 = rep(0,L)
  for(j in 1:L){
    w_hat <- w_hat_A[,j]
    X1 = rbind(Xs[[j]],X0)
    Y1 = c(ys[[j]],y0)
    p = dim(X1)[2]
    n = dim(X1)[1]
    q <- Variable(p)
    r <- Variable(n)

    objective <- Minimize((1/2) * sum((Y1  - X1 %*%q - r)^2) +  tau_x * sum(abs(q))+ tau_r* sum(abs(r)))
    # 
    problem <- Problem( objective)
    result <- solve(problem)
    x0 = result$getValue(q)
    hhat[j]  <-  sum(abs(x0 - w_hat ))

  }

  Ahat = order(hhat)[1:10]#intersect(order(hhat)[1:10],intersect(which(hhat<h),which(shat<15)))
  return(list(hhat=hhat, 
              test_data = which.min(hhat)))
}


