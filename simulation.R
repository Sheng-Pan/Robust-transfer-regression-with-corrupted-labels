
setwd('path_to_your_project') 

source('R/robust_transfer_lasso-functions.R')
library(glmnet)  
library(CVXR)
library(Rcpp)
sourceCpp('src/RTL_admm_solver_cpp.cpp')
sourceCpp('src/SDL_cpp.cpp')
load('data/cv3.Rdata')
cv_lambda0 = get("cv_lambda0")
load('data/oracal_cv.Rdata')
cv_lambda = get("cv_lambda")
N = 500
r_vector = seq(0,0.9,0.1)
ser_lasso = ser_rlasso = ser_transfer_lasso = ser_trimmed=ser_robust_transfer_lasso  = matrix(0,length(r_vector),N)
oracal_robust_transfer_lasso_signal = array(0, dim = c(400, length(r_vector),N))
ser_oracal_robust_transfer_lasso  = matrix(0,length(r_vector),N)
lasso_signal = rlasso_signal = transfer_lasso_signal = array(0, dim = c(400, length(r_vector),N))
robust_transfer_lasso_signal  = array(0, dim = c(400, length(r_vector),N))

for(u in 1:10){
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
    alpha <- 1        # α-stable噪声的尾部参数
    sigma <- 0.01     # 噪声尺度参数
    gdata = generate_target_data(m,n,p,s,sigma,k)
    X0 = Phi = gdata$X
    y0 = gdata$y
    signal_true = gdata$signal_true
    nonzero_indices = gdata$non_zero_indices
    
    ###############################
    #generate source data
    ###############################
    # 创建辅助数据集的函数
    e=8
    L = 20
    auxiliary_data <- lapply(seq_len(L), function(i) create_auxiliary_data(m+i,8,n,p))
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

    for(j in 1:L){
      cv.lasso <- cv.glmnet(Xs[[j]], ys[[j]], alpha=1)  # alpha=1 for Lasso
      fit_w <- glmnet(Xs[[j]], ys[[j]], alpha = 1,lambda = cv.lasso$lambda.min)
      w <- as.numeric(coef(fit_w)[-1])
      w_hat_A <- cbind(w_hat_A,w)
    }
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

     result = delta_estimation(Xs, ys,w_hat_A,  X0,   y0,  L,  p,  n, 1, 40)
     delta = result$delta
     hhat = result$hhat
     Ahat = intersect(order(hhat)[1:10],which(hhat<30))#which(s0<15)#intersect(J,which(s0<15))
     if(length(Ahat)==0){
       RTL  =rep(0,p)
     }else{
       test_data = which.min(hhat)
       X_test = Xs[[test_data]]
       y_test = ys[[test_data]]
       w_hat_A0 = distributed_lasso(Xs, ys, Ah, max_iter = 40, rho = 0.1)
       lambda_beta_grid <- seq(cv_lambda0[u,1]-0.01, cv_lambda0[u,1], by = 0.005)
       lambda_e_grid <- seq(cv_lambda0[u,2]-0.005, cv_lambda0[u,2], by =  0.003)
       k0 = 2
       r = 15
       ctilde = 2
       set.seed(m)
       folds <- sample(rep(1:k0, length.out = nrow(X0))) 
       # lambda_beta_grid <- seq(0.02, 0.05, by = 0.002)
       # lambda_e_grid <- seq(0.007, 0.03, by =  0.002)
       #cross-validation using R
       #cv_para = cv_RTL2(Xs,ys,X0,y0,folds,ctilde,k0,r,r2,test_data,w_hat_A0,lambda_beta_grid, lambda_e_grid)
       #cross-validation in parallel using C++
       cv_para2 = cv_RTL_cpp( folds, X_test, y_test , X0, y0, ctilde, r, k0, r2, w_hat_A0, lambda_beta_grid, lambda_e_grid)
       # load(cv_para2,'RTL_cv_k.Rdata')
       lambda_beta = cv_para$lambda_beta
       lambda_e = cv_para$lambda_e
       #method for RTL: 'admm','PDIPM'
       method = 'ADMM'
       result <- robust_trans_lasso(X0, y0,Xs,ys,method,Ahat,w_hat_A0,test_data, 25,lambda_beta,  lambda_e)
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
    Ah = intersect(order(h0)[1:10],intersect(which(h0<12),which(shat<15)))
    test_data = which.min(h0)
    w_hat_A0 = distributed_lasso(Xs, ys, Ah, max_iter = 40, rho = 0.1)
    lambda_beta_grid <- seq(cv_lambda[u,1]-0.01, cv_lambda[u,1], by = 0.005)
    lambda_e_grid <- seq(cv_lambda[u,2]-0.005, cv_lambda[u,2], by =  0.003)
    k0 = 2
    r = 15
    ctilde = 2
    set.seed(m)
    folds <- sample(rep(1:k0, length.out = nrow(X0))) 
    #cross-validation using R
    #cv_para = cv_RTL2(Xs,ys,X0,y0,folds,ctilde,k0,r,r2,test_data,w_hat_A0,lambda_beta_grid, lambda_e_grid)
    #cross-validation in parallel using C++
    cv_para2 = cv_RTL_cpp( folds, X_test, y_test , X0, y0, ctilde, r, k0, r2, w_hat_A0, lambda_beta_grid, lambda_e_grid)
    lambda_beta = cv_para$lambda_beta
    lambda_e = cv_para$lambda_e
    result <- robust_trans_lasso(X0, y0,Xs,ys,'ADMM',Ah,w_hat_A0,test_data, 12,lambda_beta,  lambda_e)
    oracalRTL = result$beta_hat
    ser_oracal_robust_transfer_lasso[u,m] = ser(signal_true,  oracalRTL)
    oracal_robust_transfer_lasso_signal[,u,m] = oracalRTL
  }
  a = data.frame(k = r_vector,
                 ser_lasso = ser_lasso[u,],
                 ser_rlasso = ser_rlasso[u,],
                 ser_transfer_lasso = ser_transfer_lasso[u,],
                 ser_robust_transfer_lasso = ser_robust_transfer_lasso[u,] ,
                 ser_oracal_robust_transfer_lasso = ser_oracal_robust_transfer_lasso[u,])
  filename <- sprintf("C:/Users/Administrator/OneDrive/ynu/papers/corrupted Data/result_%d.rds", u)  
  saveRDS(a, filename)  
  
}
