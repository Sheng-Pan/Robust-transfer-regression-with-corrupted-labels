
setwd('path_to_your_project') 
source('R/robust_transfer_lasso_functions.R')
library(glmnet)  
library(CVXR)
library(Rcpp)
library(stabledist)
library(copula)
setwd('path_to_your_project/scr')

sourceCpp('RTL_admm_solver_cpp.cpp')
sourceCpp('SDL_cpp.cpp')
setwd('C:/Users/Administrator/Downloads/github/cv_data')
####scenario 1
load('cv3.Rdata')
cv_lambda0 = get("cv_lambda0")
load('oracal_cv.Rdata')
cv_lambda = get("cv_lambda")
N = 500
r_vector = seq(0,0.9,0.1)
ser_lasso = ser_rlasso = ser_transfer_lasso = ser_trimmed=ser_robust_transfer_lasso  = matrix(0,length(r_vector),N)
oracal_robust_transfer_lasso_signal = array(0, dim = c(400, length(r_vector),N))
ser_oracal_robust_transfer_lasso  = matrix(0,length(r_vector),N)
lasso_signal = rlasso_signal = transfer_lasso_signal = array(0, dim = c(400, length(r_vector),N))
robust_transfer_lasso_signal  = array(0, dim = c(400, length(r_vector),N))

for(u in 2:10){
  k = r_vector[u]
  for(m in 1:N){
    ###############################
    #generate target dataset
    ###############################
    # ptm <- Sys.time()
    ptm <- Sys.time()
    s <- 12           
    p <- 400          
    n <- 100          
    alpha <- 1       
    sigma <- 0.01    
    gdata = generate_target_data(m,n,p,s,sigma,k)
    X0 = Phi = gdata$X
    y0 = gdata$y
    signal_true = gdata$signal_true
    nonzero_indices = gdata$non_zero_indices
    
    ###############################
    #generate source data
    ###############################
    
    e=8
    L = 20
    auxiliary_data <- lapply(seq_len(L), function(i) create_auxiliary_data(m+i,8,n,p))
    
    get_X_safe <- function(data_list, component_name) {
      if (!is.null(data_list) && is.list(data_list) && component_name %in% names(data_list)) {
        return(data_list[[component_name]])
      } else {
        stop(paste("Component", component_name, "not found in data list."))
      }
    }
    
    
    Xs <- lapply(auxiliary_data, get_X_safe, component_name = "X")
    ys <- lapply(auxiliary_data, get_X_safe, component_name = "y")
    signal <- lapply(auxiliary_data, get_X_safe, component_name = "signal")
    error <- unlist(lapply(auxiliary_data, get_X_safe, component_name = "error"))
    # ###############################
    # #distributed lasso on source data
    # ###############################
    w_hat_A = w_hat_A_debias = NULL
    nv = rep(0,L)
    for(j in 1:L){
      cv.lasso <- cv.glmnet(Xs[[j]], ys[[j]], alpha=1)  # alpha=1 for Lasso
      fit_w <- glmnet(Xs[[j]], ys[[j]], alpha = 1,lambda = cv.lasso$lambda.min)
      w <- as.numeric(coef(fit_w)[-1])
      w_hat_A <- cbind(w_hat_A,w)
      nv[j] = mean((ys[[j]] - Xs[[j]]%*%w)^2)
    }
    nv = sqrt(mean(nv))
    ###############################
    #lasso
    ###############################
    cv.lasso <- cv.glmnet(X0, y0, alpha=1)  # alpha=1 for Lasso
    fit_w <- glmnet(X0, y0, alpha = 1,lambda = cv.lasso$lambda.min)
    lasso <- as.numeric(coef(fit_w)[-1])
    ser_lasso[u,m] = ser(signal_true,  lasso)
    lasso_signal[,u,m] = lasso
    ser(signal_true,  lasso)
    # ###############################
    # #tran-lasso(Li)
    # ###############################
    # ###generate the data
    ## The following code of tran-lasso(Li) is adapted from  https://github.com/saili0103/TransLasso
    X <- NULL
    y <- NULL
    X = rbind(X, Phi)
    y = c(y,y0)
    n.vec <- c(n, rep(100, L))
    for (j in 1:(L )) {
      X <- rbind(X,Xs[[j]])
      ind.k <- ind.set(n.vec, j)
      y <- c(y, ys[[j]])
    }
    ###
    l1=T
    size.A0 = 12
    h=25
    A0 = 1:size.A0
    
    ser_transfer_lasso[u,m] <- tryCatch({
      prop.re1 <- Trans.lasso(X, y, n.vec, I.til = 1:50, l1 = l1)
      prop.re2 <- Trans.lasso(X, y, n.vec, I.til = 101:n.vec[1], l1 = l1)
      
      if(size.A0 > 0 & size.A0 < L) { #Rank.re characterizes the performance of the sparsity index Rk
        Rank.re <- (sum(prop.re1$rank.pi[1:size.A0] <= size.A0) +
                      sum(prop.re2$rank.pi[1:size.A0] <= size.A0)) / 2 / size.A0
      } else {
        Rank.re <- 1
      }
      tran_lasso <- (prop.re1$beta.hat + prop.re2$beta.hat) / 2
      ser_transfer_lasso[u,m] <- ser(signal_true, tran_lasso)
      
      
      ser_transfer_lasso[u,m]
      
    }, error = function(e) {
      
      cat("Error occurred:", e$message, "\n")
      0
    })
    ser(signal_true,  tran_lasso )
    ser_transfer_lasso[u,m] = ser(signal_true,  tran_lasso)
    transfer_lasso_signal[,u,m] = tran_lasso
    # ###############################
    #robust lasso
    ###############################
    tau_x = sqrt((log(p)/n))
    tau_r = sqrt((log(n)/n))
    #method for rlasso: 'admm','cd','PDIPM'
    method = 'ADMM'
    result0 <- Rlasso(y0, X0, tau_x, tau_r, method,rho=0.1, max_iter = 40, tol = 1e-6)
    rlasso = result0$x
    r2 = sum(abs(result0$r)>0.2*sqrt(log(n)/n))
    ser_rlasso[u,m] = ser(signal_true,  rlasso)
    rlasso_signal[,u,m] = rlasso
    ser(signal_true,  rlasso)
    # ###############################
    # #roust_trans_lasso
    # ###############################
    method = 'ADMM'
    # result = delta_estimation(Xs, ys,w_hat_A,  X0,   y0,  L,  p,  n, 1, 40)
    result = delta_estimation_fun(Xs, ys,method, tau_x, tau_r, w_hat_A,  X0,   y0,  L,  p,  n)
    delta = result$delta
    hhat = result$hhat
    Ahat = intersect(order(hhat)[1:10],which(hhat<30))#which(s0<15)#intersect(J,which(s0<15))
    if(length(Ahat)==0){
      RTL  =rep(0,p)
    }else{
      test_data = which.min(hhat)
      X_test = Xs[[test_data]]
      y_test = ys[[test_data]]
      #method = "ADMM" or "proximal"
      method = "ADMM"
      w_hat_A0 = distributed_lasso_fun(Xs, ys,method,test_data, Ahat,nv, max_iter = 40, rho = 0.1)
      lambda_beta_grid <- seq(cv_lambda0[u,1]-0.01, cv_lambda0[u,1], by = 0.005)
      lambda_e_grid <- seq(cv_lambda0[u,2]-0.005, cv_lambda0[u,2], by =  0.003)
      k0 = 2
      r0 = 15
      ctilde = 2
      set.seed(m)
      folds <- sample(rep(1:k0, length.out = nrow(X0))) 
      # lambda_beta_grid <- seq(0.02, 0.05, by = 0.002)
      # lambda_e_grid <- seq(0.007, 0.03, by =  0.002)
      #cross-validation using R
      method = 'ADMM'
      cv_para = cv_RTL2(Xs,ys,X0,y0,Ahat,folds,method,ctilde,k0,r0,r2,test_data,w_hat_A0,lambda_beta_grid, lambda_e_grid)
      #cross-validation in parallel using C++
      # cv_para2 = cv_RTL_cpp( folds, X_test, y_test , X0, y0, ctilde, r, k0, r2, w_hat_A0, lambda_beta_grid, lambda_e_grid)
      # load(cv_para2,'RTL_cv_k.Rdata')
      lambda_beta = cv_para$lambda_beta
      lambda_e = cv_para$lambda_e
      #method for RTL: 'admm','PDIPM'
      method = 'ADMM'
      result <- robust_trans_lasso(X0, y0,Ahat,Xs,ys,method,w_hat_A0,test_data, 25,lambda_beta,  lambda_e)
      RTL = result$beta_hat
    }
    ser_robust_transfer_lasso[u,m] = ser(signal_true,  RTL)
    robust_transfer_lasso_signal[,u,m] = RTL
    # 
    #  
    #  # # ###############################
    #  # # #oracal roust_trans_lasso
    #  # # ###############################
    h0 = rep(0,L)
    for(j in 1:L){
      h0[j] = sum(abs(signal_true - signal[[j]]))
    }
    Ah = which(h0<10)
    test_data = which.min(h0)
    #method = "ADMM" or "proximal"
    method = "ADMM"
    w_hat_A0 = distributed_lasso_fun(Xs, ys,method,test_data, Ah,nv, max_iter = 40, rho = 0.1)
    lambda_beta_grid <- seq(cv_lambda[u,1]-0.01, cv_lambda[u,1], by = 0.005)
    lambda_e_grid <- seq(cv_lambda[u,2]-0.005, cv_lambda[u,2], by =  0.003)
    k0 = 2
    r0 = 15
    ctilde = 2
    set.seed(m)
    folds <- sample(rep(1:k0, length.out = nrow(X0))) 
    #cross-validation using R
    method = 'ADMM'
    cv_para = cv_RTL2(Xs,ys,X0,y0,Ah,folds,method,ctilde,k0,r0,r2,test_data,w_hat_A0,lambda_beta_grid, lambda_e_grid)
    #cross-validation in parallel using C++
    # cv_para2 = cv_RTL_cpp( folds, X_test, y_test , X0, y0, ctilde, r, k0, r2, w_hat_A0, lambda_beta_grid, lambda_e_grid)
    lambda_beta = cv_para$lambda_beta
    lambda_e = cv_para$lambda_e
    result <- robust_trans_lasso(X0, y0,A,Xs,ys,'ADMM',w_hat_A0,test_data, 12,lambda_beta,  lambda_e)
    oracalRTL = result$beta_hat
    ser_oracal_robust_transfer_lasso[u,m] = ser(signal_true,  oracalRTL)
    oracal_robust_transfer_lasso_signal[,u,m] = oracalRTL
    print(c(u,m,ser(signal_true,  lasso),
            ser(signal_true,   tran_lasso),
            ser(signal_true,  RTL),
            ser(signal_true,  oracalRTL)))
  }
  a = data.frame(k = r_vector,
                 ser_lasso = apply(ser_lasso ,1,mean),
                 ser_rlasso = apply(ser_rlasso ,1,mean),
                 ser_transfer_lasso = apply(ser_transfer_lasso ,1,mean),
                 ser_robust_transfer_lasso = apply(ser_robust_transfer_lasso ,1,mean) ,
                 ser_oracal_robust_transfer_lasso = apply(ser_oracal_robust_transfer_lasso ,1,mean))
  # filename <- paste0("result_0.1_ser1",u,".rds")
  # saveRDS(a, filename)
  b=  data.frame(k = r_vector,
                 mse_lasso = apply(mse_lasso ,1,mean),
                 mse_rlasso = apply(mse_rlasso ,1,mean),
                 mse_transfer_lasso = apply(mse_transfer_lasso ,1,mean),
                 mse_robust_transfer_lasso = apply(mse_robust_transfer_lasso ,1,mean) ,
                 mse_oracal_robust_transfer_lasso = apply(mse_oracal_robust_transfer_lasso ,1,mean))
  # filename <-paste0("result_0.1_mse1",u,".rds")
  # saveRDS(b, filename)
  
}
####scenario 2
load('cv3.Rdata')
cv_lambda0 = get("cv_lambda0")
load('oracal_cv.Rdata')
cv_lambda = get("cv_lambda")
N = 500
r_vector = seq(0,0.9,0.1)
ser_lasso = ser_rlasso = ser_transfer_lasso = ser_trimmed=ser_robust_transfer_lasso  = matrix(0,length(r_vector),N)
oracal_robust_transfer_lasso_signal = array(0, dim = c(400, length(r_vector),N))
ser_oracal_robust_transfer_lasso  = matrix(0,length(r_vector),N)
lasso_signal = rlasso_signal = transfer_lasso_signal = array(0, dim = c(400, length(r_vector),N))
robust_transfer_lasso_signal  = array(0, dim = c(400, length(r_vector),N))

for(u in 2:10){
  k = r_vector[u]
  for(m in 1:N){
    ###############################
    #generate target dataset
    ###############################
    # ptm <- Sys.time()
    ptm <- Sys.time()
    s <- 12           
    p <- 400          
    n <- 100          
    alpha <- 1       
    sigma <- 0.1    
    gdata = generate_target_data(m,n,p,s,sigma,k)
    X0 = Phi = gdata$X
    y0 = gdata$y
    signal_true = gdata$signal_true
    nonzero_indices = gdata$non_zero_indices
    
    ###############################
    #generate source data
    ###############################
    
    e=8
    L = 20
    auxiliary_data <- lapply(seq_len(L), function(i) create_auxiliary_data(m+i,8,n,p))
    
    get_X_safe <- function(data_list, component_name) {
      if (!is.null(data_list) && is.list(data_list) && component_name %in% names(data_list)) {
        return(data_list[[component_name]])
      } else {
        stop(paste("Component", component_name, "not found in data list."))
      }
    }
    
    
    Xs <- lapply(auxiliary_data, get_X_safe, component_name = "X")
    ys <- lapply(auxiliary_data, get_X_safe, component_name = "y")
    signal <- lapply(auxiliary_data, get_X_safe, component_name = "signal")
    error <- unlist(lapply(auxiliary_data, get_X_safe, component_name = "error"))
    # ###############################
    # #distributed lasso on source data
    # ###############################
    w_hat_A = w_hat_A_debias = NULL
    nv = rep(0,L)
    for(j in 1:L){
      cv.lasso <- cv.glmnet(Xs[[j]], ys[[j]], alpha=1)  # alpha=1 for Lasso
      fit_w <- glmnet(Xs[[j]], ys[[j]], alpha = 1,lambda = cv.lasso$lambda.min)
      w <- as.numeric(coef(fit_w)[-1])
      w_hat_A <- cbind(w_hat_A,w)
      nv[j] = mean((ys[[j]] - Xs[[j]]%*%w)^2)
    }
    nv = sqrt(mean(nv))
    ###############################
    #lasso
    ###############################
    cv.lasso <- cv.glmnet(X0, y0, alpha=1)  # alpha=1 for Lasso
    fit_w <- glmnet(X0, y0, alpha = 1,lambda = cv.lasso$lambda.min)
    lasso <- as.numeric(coef(fit_w)[-1])
    ser_lasso[u,m] = ser(signal_true,  lasso)
    lasso_signal[,u,m] = lasso
    ser(signal_true,  lasso)
    # ###############################
    # #tran-lasso(Li)
    # ###############################
    # ###generate the data
    ## The following code of tran-lasso(Li) is adapted from  https://github.com/saili0103/TransLasso
    X <- NULL
    y <- NULL
    X = rbind(X, Phi)
    y = c(y,y0)
    n.vec <- c(n, rep(100, L))
    for (j in 1:(L )) {
      X <- rbind(X,Xs[[j]])
      ind.k <- ind.set(n.vec, j)
      y <- c(y, ys[[j]])
    }
    ###
    l1=T
    size.A0 = 12
    h=25
    A0 = 1:size.A0
    
    ser_transfer_lasso[u,m] <- tryCatch({
      prop.re1 <- Trans.lasso(X, y, n.vec, I.til = 1:50, l1 = l1)
      prop.re2 <- Trans.lasso(X, y, n.vec, I.til = 101:n.vec[1], l1 = l1)
      
      if(size.A0 > 0 & size.A0 < L) { #Rank.re characterizes the performance of the sparsity index Rk
        Rank.re <- (sum(prop.re1$rank.pi[1:size.A0] <= size.A0) +
                      sum(prop.re2$rank.pi[1:size.A0] <= size.A0)) / 2 / size.A0
      } else {
        Rank.re <- 1
      }
      tran_lasso <- (prop.re1$beta.hat + prop.re2$beta.hat) / 2
      ser_transfer_lasso[u,m] <- ser(signal_true, tran_lasso)
      
      
      ser_transfer_lasso[u,m]
      
    }, error = function(e) {
      
      cat("Error occurred:", e$message, "\n")
      0
    })
    ser(signal_true,  tran_lasso )
    ser_transfer_lasso[u,m] = ser(signal_true,  tran_lasso)
    transfer_lasso_signal[,u,m] = tran_lasso
    # ###############################
    #robust lasso
    ###############################
    tau_x = sqrt((log(p)/n))
    tau_r = sqrt((log(n)/n))
    #method for rlasso: 'admm','cd','PDIPM'
    method = 'ADMM'
    result0 <- Rlasso(y0, X0, tau_x, tau_r, method,rho=0.1, max_iter = 40, tol = 1e-6)
    rlasso = result0$x
    r2 = sum(abs(result0$r)>0.2*sqrt(log(n)/n))
    ser_rlasso[u,m] = ser(signal_true,  rlasso)
    rlasso_signal[,u,m] = rlasso
    ser(signal_true,  rlasso)
    # ###############################
    # #roust_trans_lasso
    # ###############################
    method = 'ADMM'
    # result = delta_estimation(Xs, ys,w_hat_A,  X0,   y0,  L,  p,  n, 1, 40)
    result = delta_estimation_fun(Xs, ys,method, tau_x, tau_r, w_hat_A,  X0,   y0,  L,  p,  n)
    delta = result$delta
    hhat = result$hhat
    Ahat = intersect(order(hhat)[1:10],which(hhat<30))#which(s0<15)#intersect(J,which(s0<15))
    if(length(Ahat)==0){
      RTL  =rep(0,p)
    }else{
      test_data = which.min(hhat)
      X_test = Xs[[test_data]]
      y_test = ys[[test_data]]
      #method = "ADMM" or "proximal"
      method = "ADMM"
      w_hat_A0 = distributed_lasso_fun(Xs, ys,method,test_data, Ahat,nv, max_iter = 40, rho = 0.1)
      lambda_beta_grid <- seq(cv_lambda0[u,1]-0.01, cv_lambda0[u,1], by = 0.005)
      lambda_e_grid <- seq(cv_lambda0[u,2]-0.005, cv_lambda0[u,2], by =  0.003)
      k0 = 2
      r0 = 15
      ctilde = 2
      set.seed(m)
      folds <- sample(rep(1:k0, length.out = nrow(X0))) 
      # lambda_beta_grid <- seq(0.02, 0.05, by = 0.002)
      # lambda_e_grid <- seq(0.007, 0.03, by =  0.002)
      #cross-validation using R
      method = 'ADMM'
      cv_para = cv_RTL2(Xs,ys,X0,y0,Ahat,folds,method,ctilde,k0,r0,r2,test_data,w_hat_A0,lambda_beta_grid, lambda_e_grid)
      #cross-validation in parallel using C++
      # cv_para2 = cv_RTL_cpp( folds, X_test, y_test , X0, y0, ctilde, r, k0, r2, w_hat_A0, lambda_beta_grid, lambda_e_grid)
      # load(cv_para2,'RTL_cv_k.Rdata')
      lambda_beta = cv_para$lambda_beta
      lambda_e = cv_para$lambda_e
      #method for RTL: 'admm','PDIPM'
      method = 'ADMM'
      result <- robust_trans_lasso(X0, y0,Ahat,Xs,ys,method,w_hat_A0,test_data, 25,lambda_beta,  lambda_e)
      RTL = result$beta_hat
    }
    ser_robust_transfer_lasso[u,m] = ser(signal_true,  RTL)
    robust_transfer_lasso_signal[,u,m] = RTL
    # 
    #  
    #  # # ###############################
    #  # # #oracal roust_trans_lasso
    #  # # ###############################
    h0 = rep(0,L)
    for(j in 1:L){
      h0[j] = sum(abs(signal_true - signal[[j]]))
    }
    Ah = which(h0<10)
    test_data = which.min(h0)
    #method = "ADMM" or "proximal"
    method = "ADMM"
    w_hat_A0 = distributed_lasso_fun(Xs, ys,method,test_data, Ah,nv, max_iter = 40, rho = 0.1)
    lambda_beta_grid <- seq(cv_lambda[u,1]-0.01, cv_lambda[u,1], by = 0.005)
    lambda_e_grid <- seq(cv_lambda[u,2]-0.005, cv_lambda[u,2], by =  0.003)
    k0 = 2
    r0 = 15
    ctilde = 2
    set.seed(m)
    folds <- sample(rep(1:k0, length.out = nrow(X0))) 
    #cross-validation using R
    method = 'ADMM'
    cv_para = cv_RTL2(Xs,ys,X0,y0,Ah,folds,method,ctilde,k0,r0,r2,test_data,w_hat_A0,lambda_beta_grid, lambda_e_grid)
    #cross-validation in parallel using C++
    # cv_para2 = cv_RTL_cpp( folds, X_test, y_test , X0, y0, ctilde, r, k0, r2, w_hat_A0, lambda_beta_grid, lambda_e_grid)
    lambda_beta = cv_para$lambda_beta
    lambda_e = cv_para$lambda_e
    result <- robust_trans_lasso(X0, y0,A,Xs,ys,'ADMM',w_hat_A0,test_data, 12,lambda_beta,  lambda_e)
    oracalRTL = result$beta_hat
    ser_oracal_robust_transfer_lasso[u,m] = ser(signal_true,  oracalRTL)
    oracal_robust_transfer_lasso_signal[,u,m] = oracalRTL
    print(c(u,m,ser(signal_true,  lasso),
            ser(signal_true,   tran_lasso),
            ser(signal_true,  RTL),
            ser(signal_true,  oracalRTL)))
  }
  a = data.frame(k = r_vector,
                 ser_lasso = apply(ser_lasso ,1,mean),
                 ser_rlasso = apply(ser_rlasso ,1,mean),
                 ser_transfer_lasso = apply(ser_transfer_lasso ,1,mean),
                 ser_robust_transfer_lasso = apply(ser_robust_transfer_lasso ,1,mean) ,
                 ser_oracal_robust_transfer_lasso = apply(ser_oracal_robust_transfer_lasso ,1,mean))
  # filename <- paste0("result_0.1_ser1",u,".rds")
  # saveRDS(a, filename)
  b=  data.frame(k = r_vector,
                 mse_lasso = apply(mse_lasso ,1,mean),
                 mse_rlasso = apply(mse_rlasso ,1,mean),
                 mse_transfer_lasso = apply(mse_transfer_lasso ,1,mean),
                 mse_robust_transfer_lasso = apply(mse_robust_transfer_lasso ,1,mean) ,
                 mse_oracal_robust_transfer_lasso = apply(mse_oracal_robust_transfer_lasso ,1,mean))
  # filename <-paste0("result_0.1_mse1",u,".rds")
  # saveRDS(b, filename)
  
}
###scenario 3
load('cv3_xb.Rdata')
cv_lambda0 = get("cv_lambda0")
load('oracal_cv_xb.Rdata')
cv_lambda = get("cv_lambda")
N = 500
r_vector = seq(0.1,0.9,0.1)
ser_lasso = ser_rlasso = ser_transfer_lasso = ser_trimmed=ser_robust_transfer_lasso  = matrix(0,length(r_vector),N)
oracal_robust_transfer_lasso_signal = array(0, dim = c(400, length(r_vector),N))
ser_oracal_robust_transfer_lasso  = matrix(0,length(r_vector),N)
lasso_signal = rlasso_signal = transfer_lasso_signal = array(0, dim = c(400, length(r_vector),N))
robust_transfer_lasso_signal  = array(0, dim = c(400, length(r_vector),N))
mse_lasso = mse_rlasso = mse_transfer_lasso = mse_trimmed=mse_robust_transfer_lasso  = matrix(0,length(r_vector),N)
L1_lasso = L1_rlasso = L1_transfer_lasso =L1_robust_transfer_lasso  = matrix(0,length(r_vector),N)
mse_oracal_robust_transfer_lasso = L1_oracal_robust_transfer_lasso = matrix(0,length(r_vector),N)
  
for(u in 1:length(r_vector)){
  k = r_vector[u]
  for(m in 1:N){
    ###############################
    #generate target dataset
    ###############################
    # ptm <- Sys.time()
    ptm <- Sys.time()
    s <- 12           # 稀疏度
    p <- 400          # 信号长度
    n <- 100          # 测量数量
    sigma <- 1   # 噪声尺度参数
    distr = 'xb'
    rr = 0.3
    sce = 2
    gdata = generate_target_data2(m,sce,n,p,s,rr,sigma,k,distr)
    X0 = Phi = gdata$X
    X0 = scale(X0)
    y0 = gdata$y
    signal_true = gdata$signal_true
    nonzero_indices = gdata$non_zero_indices
    
    ###############################
    #generate source data
    ###############################
    # 创建辅助数据集的函数
    e=8
    L = 20
    auxiliary_data <- lapply(seq_len(L), function(i) create_auxiliary_data2(m+i,sce,rr,sigma,8,n,p))
    # 安全地提取 X 和 y
    get_X_safe <- function(data_list, component_name) {
      if (!is.null(data_list) && is.list(data_list) && component_name %in% names(data_list)) {
        return(data_list[[component_name]])
      } else {
        stop(paste("Component", component_name, "not found in data list."))
      }
    }
    
    # 提取所有辅助数据集中的 X 和 y
    Xs <- lapply(auxiliary_data, get_X_safe, component_name = "X")
    ys <- lapply(auxiliary_data, get_X_safe, component_name = "y")
    signal <- lapply(auxiliary_data, get_X_safe, component_name = "signal")
    error <- unlist(lapply(auxiliary_data, get_X_safe, component_name = "error"))
    # ###############################
    # #distributed lasso on source data
    # ###############################
    w_hat_A = w_hat_A_debias = NULL
    nv = rep(0,L)
    for(j in 1:L){
      cv.lasso <- cv.glmnet(Xs[[j]], ys[[j]], alpha=1)  # alpha=1 for Lasso
      fit_w <- glmnet(Xs[[j]], ys[[j]], alpha = 1,lambda = cv.lasso$lambda.min)
      w <- as.numeric(coef(fit_w)[-1])
      w_hat_A <- cbind(w_hat_A,w)
      nv[j] = mean((ys[[j]] - Xs[[j]]%*%w)^2)
    }
    nv = sqrt(mean(nv))
    ###############################
    #lasso
    ###############################
    cv.lasso <- cv.glmnet(X0, y0, alpha=1)  # alpha=1 for Lasso
    fit_w <- glmnet(X0, y0, alpha = 1,lambda = cv.lasso$lambda.min)
    lasso <- as.numeric(coef(fit_w)[-1])
    ser(signal_true,  lasso)
    ser_lasso[u,m] = ser(signal_true,  lasso)
    mse_lasso[u,m] = msefun(signal_true,  lasso)
    L1_lasso[u,m] = l1fun(signal_true,  lasso)
    lasso_signal[,u,m] = lasso
    ser(signal_true,  lasso)
    # ###############################
    # #tran-lasso(Li)
    # ###############################
    ###generate the data
    X <- NULL
    y <- NULL
    X = rbind(X, Phi)
    y = c(y,y0)
    n.vec <- c(n, rep(100, L))
    for (j in 1:(L )) {
      X <- rbind(X,Xs[[j]])
      ind.k <- ind.set(n.vec, j)
      y <- c(y, ys[[j]])
    }
    ###
    l1=T
    size.A0 = 12
    h=25
    A0 = 1:size.A0

    tran_lasso <-tryCatch({
      prop.re1 <- Trans.lasso(X, y, n.vec, I.til = 1:50, l1 = l1)
      prop.re2 <- Trans.lasso(X, y, n.vec, I.til = 101:n.vec[1], l1 = l1)

      if(size.A0 > 0 & size.A0 < L) { #Rank.re characterizes the performance of the sparsity index Rk
        Rank.re <- (sum(prop.re1$rank.pi[1:size.A0] <= size.A0) +
                      sum(prop.re2$rank.pi[1:size.A0] <= size.A0)) / 2 / size.A0
      } else {
        Rank.re <- 1
      }
      (prop.re1$beta.hat + prop.re2$beta.hat) / 2
      }, error = function(e) {

      cat("Error occurred:", e$message, "\n")
      return(rep(0,p))
    })
    ser(signal_true,  tran_lasso )
    ser_transfer_lasso[u,m] = ser(signal_true,  tran_lasso)
    transfer_lasso_signal[,u,m] = tran_lasso
    mse_transfer_lasso[u,m] = msefun(signal_true,  tran_lasso)
    L1_transfer_lasso[u,m] = l1fun(signal_true,  tran_lasso)
    # ###############################
    #robust lasso
    ###############################
     tau_x = 10*n*nv*sqrt((log(p)/n))
     tau_r = 2*n*nv*sqrt((log(n)/n))
     #method for rlasso: 'admm','cd','PDIPM'
     method = 'ADMM'
     result0 <- Rlasso(y0, X0, tau_x, tau_r, method,rho=0.1, max_iter = 40, tol = 1e-6)
     rlasso = result0$x
     r2 = sum(abs(result0$r)>0.2*sqrt(log(n)/n))
     ser_rlasso[u,m] = ser(signal_true,  rlasso)
     mse_rlasso[u,m] = msefun(signal_true,  rlasso)
     L1_rlasso[u,m] = l1fun(signal_true,  rlasso)
     rlasso_signal[,u,m] = rlasso
     ser(signal_true,  rlasso)
    # ###############################
    # #roust_trans_lasso
    # ###############################
     h0 = rep(0,L)
     for(j in 1:L){
       h0[j] = sum(abs(signal_true - signal[[j]]))
     }
     tau_x = 10*2*n*nv*sqrt((log(p)/n))
     tau_r = 2*n*nv*sqrt((log(n)/n))
     result = tryCatch({
      
       method = 'PDIPM'
       # result = delta_estimation(Xs, ys,w_hat_A,  X0,   y0,  L,  p,  n, 1, 40)
       delta_estimation_fun(Xs, ys,method, tau_x, tau_r, w_hat_A,  X0,   y0,  L,  p,  n)
     }, error = function(e) {
       
       cat("Error occurred:", e$message, "\n")
       return( delta_estimation_fun(Xs, ys,'ADMM', tau_x, tau_r, w_hat_A,  X0,   y0,  L,  p,  n))
     })
     #delta = result$delta
     hhat = result$hhat

     Ahat = order(hhat)[1:10]#intersect(order(hhat)[1:10],which(hhat<30))#which(s0<15)#intersect(J,which(s0<15))
     h0[Ahat]
     if(length(Ahat)==0){
       RTL  =rep(0,p)
     }else{
       test_data = which.min(hhat)
       X_test = Xs[[test_data]]
       y_test = ys[[test_data]]
      # w_hat_A0 = distributed_lasso(Xs, ys,test_data,nv, Ahat, max_iter = 40, rho = 0.1)
       if(length(Ahat)==1){
         cv.lasso <- cv.glmnet(Xs[[1]], ys[[1]], alpha=1)  # alpha=1 for Lasso
         fit_w <- glmnet(Xs[[1]], ys[[1]], alpha = 1, lambda = cv.lasso$lambda.min, standardize = FALSE)
         lasso <- as.numeric(coef(fit_w)[-1])
         return(lasso)
       }else{
        
         w_hat_A0 =  distributed_lasso(Xs, ys,test_data, Ahat,nv, max_iter = 40, rho = 0.1) 
       }
       ser(signal_true,  w_hat_A0)
       lambda_beta_grid <- seq(cv_lambda0[u,1]-0.01, cv_lambda0[u,1], by = 0.05)
       lambda_e_grid <- seq(cv_lambda0[u,2]-0.005, cv_lambda0[u,2], by =  0.003)
       k0 = 2
       r = 5
       ctilde = 2
       set.seed(m)
       folds <- sample(rep(1:k0, length.out = nrow(X0)))
       # lambda_beta_grid <- seq(0.02, 0.05, by = 0.002)
       # lambda_e_grid <- seq(0.007, 0.03, by =  0.002)
       #cross-validation using R
       method = 'PDIPM'
       cv_para = cv_RTL2(Xs,ys,X0,y0,Ahat,folds,method,ctilde,k0,r,r2,test_data,w_hat_A0,lambda_beta_grid, lambda_e_grid)
       #cross-validation in parallel using C++
      # cv_para2 = cv_RTL_cpp( folds, X_test, y_test , X0, y0, ctilde, r, k0, r2, w_hat_A0, lambda_beta_grid, lambda_e_grid)
       # load(cv_para2,'RTL_cv_k.Rdata')
       lambda_beta = cv_para$lambda_beta
       lambda_e = cv_para$lambda_e
       #method for RTL: 'admm','PDIPM'

       result <- robust_trans_lasso(X0, y0,Ahat,Xs,ys,method,w_hat_A0,test_data, 25,lambda_beta,  lambda_e)
       RTL = result$beta_hat
     }
     ser(signal_true,  RTL)
    ser_robust_transfer_lasso[u,m] = ser(signal_true,  RTL)
    robust_transfer_lasso_signal[,u,m] = RTL
    mse_robust_transfer_lasso[u,m] = msefun(signal_true,  RTL)
    L1_robust_transfer_lasso[u,m] = l1fun(signal_true,  RTL)
    # 
    #  
    #  # # ###############################
    #  # # #oracal roust_trans_lasso
    #  # # ###############################
    h0 = rep(0,L)
    for(j in 1:L){
      h0[j] = sum(abs(signal_true - signal[[j]]))
    }
    Ah = which(h0<10)
    test_data = which.min(h0)
    #w_hat_A0 = distributed_lasso(Xs, ys, Ah,nv, max_iter = 40, rho = 0.1)
    w_hat_A0 =  distributed_lasso(Xs, ys,test_data, Ahat,nv, max_iter = 40, rho = 0.1) 
  ser(signal_true,  w_hat_A0)
    lambda_beta_grid <- seq(cv_lambda[u,1]-0.01, cv_lambda[u,1], by = 0.005)
    lambda_e_grid <- seq(cv_lambda[u,2]-0.005, cv_lambda[u,2], by =  0.003)
    k0 = 2
    r = 15
    ctilde = 2
    set.seed(m)
    folds <- sample(rep(1:k0, length.out = nrow(X0)))
    #cross-validation using R
    cv_para = cv_RTL2(Xs,ys,X0,y0,Ah,folds,method,ctilde,k0,r,r2,test_data,w_hat_A0,lambda_beta_grid, lambda_e_grid)
    #cross-validation in parallel using C++
    #cv_para2 = cv_RTL_cpp( folds, X_test, y_test , X0, y0, ctilde, r, k0, r2, w_hat_A0, lambda_beta_grid, lambda_e_grid)
    lambda_beta = cv_para$lambda_beta
    lambda_e = cv_para$lambda_e
    result <- robust_trans_lasso(X0, y0,Ah,Xs,ys,'PDIPM',w_hat_A0,test_data, 12,lambda_beta,  lambda_e)
    oracalRTL = result$beta_hat
    ser_oracal_robust_transfer_lasso[u,m] = ser(signal_true,  oracalRTL)
    mse_oracal_robust_transfer_lasso[u,m] = msefun(signal_true,  oracalRTL)
    L1_oracal_robust_transfer_lasso[u,m] = l1fun(signal_true,  oracalRTL)
    oracal_robust_transfer_lasso_signal[,u,m] = oracalRTL
    print(c(u,m,ser(signal_true,  lasso),
            ser(signal_true,   tran_lasso),
            ser(signal_true,  RTL),
            ser(signal_true,  oracalRTL)))
    ptm1 <- Sys.time()
    print(ptm1 -ptm )
  }
  a = data.frame(k = r_vector,
                 ser_lasso = apply(ser_lasso ,1,mean),
                 ser_rlasso = apply(ser_rlasso ,1,mean),
                 ser_transfer_lasso = apply(ser_transfer_lasso ,1,mean),
                 ser_robust_transfer_lasso = apply(ser_robust_transfer_lasso ,1,mean) ,
                 ser_oracal_robust_transfer_lasso = apply(ser_oracal_robust_transfer_lasso ,1,mean))
  # filename <- paste0("result_0.1_ser1",u,".rds")
  # saveRDS(a, filename)
  b=  data.frame(k = r_vector,
                 mse_lasso = apply(mse_lasso ,1,mean),
                 mse_rlasso = apply(mse_rlasso ,1,mean),
                 mse_transfer_lasso = apply(mse_transfer_lasso ,1,mean),
                 mse_robust_transfer_lasso = apply(mse_robust_transfer_lasso ,1,mean) ,
                 mse_oracal_robust_transfer_lasso = apply(mse_oracal_robust_transfer_lasso ,1,mean))
  # filename <-paste0("result_0.1_mse1",u,".rds")
  # saveRDS(b, filename)

  
}
###scenario 4
load('cv3_0.1.Rdata')
cv_lambda0 = get("cv_lambda0")
load('oracal_cv_0.1.Rdata')
cv_lambda = get("cv_lambda")
N = 500
r_vector = seq(0.1,0.9,0.1)
ser_lasso = ser_rlasso = ser_transfer_lasso = ser_trimmed=ser_robust_transfer_lasso  = matrix(0,length(r_vector),N)
oracal_robust_transfer_lasso_signal = array(0, dim = c(400, length(r_vector),N))
ser_oracal_robust_transfer_lasso  = matrix(0,length(r_vector),N)
lasso_signal = rlasso_signal = transfer_lasso_signal = array(0, dim = c(400, length(r_vector),N))
robust_transfer_lasso_signal  = array(0, dim = c(400, length(r_vector),N))
mse_lasso = mse_rlasso = mse_transfer_lasso = mse_trimmed=mse_robust_transfer_lasso  = matrix(0,length(r_vector),N)
L1_lasso = L1_rlasso = L1_transfer_lasso =L1_robust_transfer_lasso  = matrix(0,length(r_vector),N)
mse_oracal_robust_transfer_lasso = L1_oracal_robust_transfer_lasso = matrix(0,length(r_vector),N)

for(u in 1:length(r_vector)){
  k = r_vector[u]
  for(m in 1:N){
    ###############################
    #generate target dataset
    ###############################
    # ptm <- Sys.time()
    ptm <- Sys.time()
    s <- 12           # 稀疏度
    p <- 400          # 信号长度
    n <- 100          # 测量数量
    sigma <- 0.1  # 噪声尺度参数
    distr = 'xb'
    rr = 0.3
    sce = 3
    gdata = generate_target_data2(m,sce,n,p,s,rr,sigma,k,distr)
    X0 = Phi = gdata$X
    X0 = scale(X0)
    y0 = gdata$y
    signal_true = gdata$signal_true
    nonzero_indices = gdata$non_zero_indices
    
    ###############################
    #generate source data
    ###############################
    # 创建辅助数据集的函数
    e=8
    L = 20
    auxiliary_data <- lapply(seq_len(L), function(i) create_auxiliary_data2(m+i,sce,rr,sigma,8,n,p))
    # 安全地提取 X 和 y
    get_X_safe <- function(data_list, component_name) {
      if (!is.null(data_list) && is.list(data_list) && component_name %in% names(data_list)) {
        return(data_list[[component_name]])
      } else {
        stop(paste("Component", component_name, "not found in data list."))
      }
    }
    
    # 提取所有辅助数据集中的 X 和 y
    Xs <- lapply(auxiliary_data, get_X_safe, component_name = "X")
    ys <- lapply(auxiliary_data, get_X_safe, component_name = "y")
    signal <- lapply(auxiliary_data, get_X_safe, component_name = "signal")
    error <- unlist(lapply(auxiliary_data, get_X_safe, component_name = "error"))
    # ###############################
    # #distributed lasso on source data
    # ###############################
    w_hat_A = w_hat_A_debias = NULL
    nv = rep(0,L)
    for(j in 1:L){
      cv.lasso <- cv.glmnet(Xs[[j]], ys[[j]], alpha=1)  # alpha=1 for Lasso
      fit_w <- glmnet(Xs[[j]], ys[[j]], alpha = 1,lambda = cv.lasso$lambda.min)
      w <- as.numeric(coef(fit_w)[-1])
      w_hat_A <- cbind(w_hat_A,w)
      nv[j] = mean((ys[[j]] - Xs[[j]]%*%w)^2)
    }
    nv = sqrt(mean(nv))
    ###############################
    #lasso
    ###############################
    cv.lasso <- cv.glmnet(X0, y0, alpha=1)  # alpha=1 for Lasso
    fit_w <- glmnet(X0, y0, alpha = 1,lambda = cv.lasso$lambda.min)
    lasso <- as.numeric(coef(fit_w)[-1])
    ser(signal_true,  lasso)
    ser_lasso[u,m] = ser(signal_true,  lasso)
    mse_lasso[u,m] = msefun(signal_true,  lasso)
    L1_lasso[u,m] = l1fun(signal_true,  lasso)
    lasso_signal[,u,m] = lasso
    ser(signal_true,  lasso)
    # ###############################
    # #tran-lasso(Li)
    # ###############################
    ###generate the data
    X <- NULL
    y <- NULL
    X = rbind(X, Phi)
    y = c(y,y0)
    n.vec <- c(n, rep(100, L))
    for (j in 1:(L )) {
      X <- rbind(X,Xs[[j]])
      ind.k <- ind.set(n.vec, j)
      y <- c(y, ys[[j]])
    }
    ###
    l1=T
    size.A0 = 12
    h=25
    A0 = 1:size.A0
    
    tran_lasso <-tryCatch({
      prop.re1 <- Trans.lasso(X, y, n.vec, I.til = 1:50, l1 = l1)
      prop.re2 <- Trans.lasso(X, y, n.vec, I.til = 101:n.vec[1], l1 = l1)
      
      if(size.A0 > 0 & size.A0 < L) { #Rank.re characterizes the performance of the sparsity index Rk
        Rank.re <- (sum(prop.re1$rank.pi[1:size.A0] <= size.A0) +
                      sum(prop.re2$rank.pi[1:size.A0] <= size.A0)) / 2 / size.A0
      } else {
        Rank.re <- 1
      }
      (prop.re1$beta.hat + prop.re2$beta.hat) / 2
    }, error = function(e) {
      
      cat("Error occurred:", e$message, "\n")
      return(rep(0,p))
    })
    ser(signal_true,  tran_lasso )
    ser_transfer_lasso[u,m] = ser(signal_true,  tran_lasso)
    transfer_lasso_signal[,u,m] = tran_lasso
    mse_transfer_lasso[u,m] = msefun(signal_true,  tran_lasso)
    L1_transfer_lasso[u,m] = l1fun(signal_true,  tran_lasso)
    # ###############################
    #robust lasso
    ###############################
    tau_x = 10*n*nv*sqrt((log(p)/n))
    tau_r = 2*n*nv*sqrt((log(n)/n))
    #method for rlasso: 'admm','cd','PDIPM'
    method = 'ADMM'
    result0 <- Rlasso(y0, X0, tau_x, tau_r, method,rho=0.1, max_iter = 40, tol = 1e-6)
    rlasso = result0$x
    r2 = sum(abs(result0$r)>0.2*sqrt(log(n)/n))
    ser_rlasso[u,m] = ser(signal_true,  rlasso)
    mse_rlasso[u,m] = msefun(signal_true,  rlasso)
    L1_rlasso[u,m] = l1fun(signal_true,  rlasso)
    rlasso_signal[,u,m] = rlasso
    ser(signal_true,  rlasso)
    # ###############################
    # #roust_trans_lasso
    # ###############################
    h0 = rep(0,L)
    for(j in 1:L){
      h0[j] = sum(abs(signal_true - signal[[j]]))
    }
    tau_x = 10*2*n*nv*sqrt((log(p)/n))
    tau_r = 2*n*nv*sqrt((log(n)/n))
    result = tryCatch({
      
      method = 'PDIPM'
      # result = delta_estimation(Xs, ys,w_hat_A,  X0,   y0,  L,  p,  n, 1, 40)
      delta_estimation_fun(Xs, ys,method, tau_x, tau_r, w_hat_A,  X0,   y0,  L,  p,  n)
    }, error = function(e) {
      
      cat("Error occurred:", e$message, "\n")
      return( delta_estimation_fun(Xs, ys,'ADMM', tau_x, tau_r, w_hat_A,  X0,   y0,  L,  p,  n))
    })
    #delta = result$delta
    hhat = result$hhat
    
    Ahat = order(hhat)[1:10]#intersect(order(hhat)[1:10],which(hhat<30))#which(s0<15)#intersect(J,which(s0<15))
    h0[Ahat]
    if(length(Ahat)==0){
      RTL  =rep(0,p)
    }else{
      test_data = which.min(hhat)
      X_test = Xs[[test_data]]
      y_test = ys[[test_data]]
      # w_hat_A0 = distributed_lasso(Xs, ys,test_data,nv, Ahat, max_iter = 40, rho = 0.1)
      if(length(Ahat)==1){
        cv.lasso <- cv.glmnet(Xs[[1]], ys[[1]], alpha=1)  # alpha=1 for Lasso
        fit_w <- glmnet(Xs[[1]], ys[[1]], alpha = 1, lambda = cv.lasso$lambda.min, standardize = FALSE)
        lasso <- as.numeric(coef(fit_w)[-1])
        return(lasso)
      }else{
        
        w_hat_A0 =  distributed_lasso(Xs, ys,test_data, Ahat,nv, max_iter = 40, rho = 0.1) 
      }
      ser(signal_true,  w_hat_A0)
      lambda_beta_grid <- seq(cv_lambda0[u,1]-0.01, cv_lambda0[u,1], by = 0.05)
      lambda_e_grid <- seq(cv_lambda0[u,2]-0.005, cv_lambda0[u,2], by =  0.003)
      k0 = 2
      r = 5
      ctilde = 2
      set.seed(m)
      folds <- sample(rep(1:k0, length.out = nrow(X0)))
      # lambda_beta_grid <- seq(0.02, 0.05, by = 0.002)
      # lambda_e_grid <- seq(0.007, 0.03, by =  0.002)
      #cross-validation using R
      method = 'PDIPM'
      cv_para = cv_RTL2(Xs,ys,X0,y0,Ahat,folds,method,ctilde,k0,r,r2,test_data,w_hat_A0,lambda_beta_grid, lambda_e_grid)
      #cross-validation in parallel using C++
      # cv_para2 = cv_RTL_cpp( folds, X_test, y_test , X0, y0, ctilde, r, k0, r2, w_hat_A0, lambda_beta_grid, lambda_e_grid)
      # load(cv_para2,'RTL_cv_k.Rdata')
      lambda_beta = cv_para$lambda_beta
      lambda_e = cv_para$lambda_e
      #method for RTL: 'admm','PDIPM'
      
      result <- robust_trans_lasso(X0, y0,Ahat,Xs,ys,method,w_hat_A0,test_data, 25,lambda_beta,  lambda_e)
      RTL = result$beta_hat
    }
    ser(signal_true,  RTL)
    ser_robust_transfer_lasso[u,m] = ser(signal_true,  RTL)
    robust_transfer_lasso_signal[,u,m] = RTL
    mse_robust_transfer_lasso[u,m] = msefun(signal_true,  RTL)
    L1_robust_transfer_lasso[u,m] = l1fun(signal_true,  RTL)
    # 
    #  
    #  # # ###############################
    #  # # #oracal roust_trans_lasso
    #  # # ###############################
    h0 = rep(0,L)
    for(j in 1:L){
      h0[j] = sum(abs(signal_true - signal[[j]]))
    }
    Ah = which(h0<10)
    test_data = which.min(h0)
    #w_hat_A0 = distributed_lasso(Xs, ys, Ah,nv, max_iter = 40, rho = 0.1)
    w_hat_A0 =  distributed_lasso(Xs, ys,test_data, Ahat,nv, max_iter = 40, rho = 0.1) 
    ser(signal_true,  w_hat_A0)
    lambda_beta_grid <- seq(cv_lambda[u,1]-0.01, cv_lambda[u,1], by = 0.005)
    lambda_e_grid <- seq(cv_lambda[u,2]-0.005, cv_lambda[u,2], by =  0.003)
    k0 = 2
    r = 15
    ctilde = 2
    set.seed(m)
    folds <- sample(rep(1:k0, length.out = nrow(X0)))
    #cross-validation using R
    cv_para = cv_RTL2(Xs,ys,X0,y0,Ah,folds,method,ctilde,k0,r,r2,test_data,w_hat_A0,lambda_beta_grid, lambda_e_grid)
    #cross-validation in parallel using C++
    #cv_para2 = cv_RTL_cpp( folds, X_test, y_test , X0, y0, ctilde, r, k0, r2, w_hat_A0, lambda_beta_grid, lambda_e_grid)
    lambda_beta = cv_para$lambda_beta
    lambda_e = cv_para$lambda_e
    result <- robust_trans_lasso(X0, y0,Ah,Xs,ys,'PDIPM',w_hat_A0,test_data, 12,lambda_beta,  lambda_e)
    oracalRTL = result$beta_hat
    ser_oracal_robust_transfer_lasso[u,m] = ser(signal_true,  oracalRTL)
    mse_oracal_robust_transfer_lasso[u,m] = msefun(signal_true,  oracalRTL)
    L1_oracal_robust_transfer_lasso[u,m] = l1fun(signal_true,  oracalRTL)
    oracal_robust_transfer_lasso_signal[,u,m] = oracalRTL
    print(c(u,m,ser(signal_true,  lasso),
            ser(signal_true,   tran_lasso),
            ser(signal_true,  RTL),
            ser(signal_true,  oracalRTL)))
    ptm1 <- Sys.time()
    print(ptm1 -ptm )
  }
  a = data.frame(k = r_vector,
                 ser_lasso = apply(ser_lasso ,1,mean),
                 ser_rlasso = apply(ser_rlasso ,1,mean),
                 ser_transfer_lasso = apply(ser_transfer_lasso ,1,mean),
                 ser_robust_transfer_lasso = apply(ser_robust_transfer_lasso ,1,mean) ,
                 ser_oracal_robust_transfer_lasso = apply(ser_oracal_robust_transfer_lasso ,1,mean))
  # filename <- paste0("result_0.1_ser1",u,".rds")
  # saveRDS(a, filename)
  b=  data.frame(k = r_vector,
                 mse_lasso = apply(mse_lasso ,1,mean),
                 mse_rlasso = apply(mse_rlasso ,1,mean),
                 mse_transfer_lasso = apply(mse_transfer_lasso ,1,mean),
                 mse_robust_transfer_lasso = apply(mse_robust_transfer_lasso ,1,mean) ,
                 mse_oracal_robust_transfer_lasso = apply(mse_oracal_robust_transfer_lasso ,1,mean))
  # filename <-paste0("result_0.1_mse1",u,".rds")
  # saveRDS(b, filename)
  
  
}

