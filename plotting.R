
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
for(u in 1:length(r_vector)){
  k = r_vector[u]
  ###############################
  #generate target dataset
  ###############################
  s <- 12           
  p <- 400        
  n <- 100         
  alpha <- 1       
  sigma <- 0.01    
  
  gdata = generate_target_data(123,n,p,s,sigma,k)
  X0 = Phi = gdata$X
  y0 = gdata$y
  signal_true = gdata$signal_true
  nonzero_indices = gdata$non_zero_indices
  
  ###############################
  #generate source data
  ###############################
  e=8
  L = 20
  auxiliary_data <- lapply(seq_len(L), function(i) create_auxiliary_data(i+u,8,n,p))
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
  ###############################
  #lasso
  ###############################
  cv.lasso <- cv.glmnet(X0, y0, alpha=1)  # alpha=1 for Lasso
  fit_w <- glmnet(X0, y0, alpha = 1,lambda = cv.lasso$lambda.min)
  lasso <- as.numeric(coef(fit_w)[-1])
  ser(signal_true,  lasso)
  ###############################
  #tran-lasso(Li)
  ###############################
  ###generate the data
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
  h=6
  A0 = 1:size.A0
  
  prop.re1 <- Trans.lasso(X, y, n.vec, I.til = 1:50, l1 = l1)
  prop.re2 <- Trans.lasso(X, y, n.vec, I.til = 101:n.vec[1], l1=l1)
  if(size.A0 > 0 & size.A0< L){ #Rank.re characterizes the performance of the sparsity index Rk
    Rank.re<- (sum(prop.re1$rank.pi[1:size.A0]<=size.A0) +
                 sum(prop.re2$rank.pi[1:size.A0]<=size.A0))/2/size.A0
  }else{ Rank.re <- 1 }
  tran_lasso <- (prop.re1$beta.hat + prop.re2$beta.hat) / 2
  ser(signal_true,  tran_lasso )
  ###############################
  #robust lasso
  ###############################
  ##
  tau_x = sqrt((log(p)/n))
  tau_r = sqrt((log(n)/n))
  #method for rlasso: 'admm','cd','PDIPM'
  method = 'ADMM'
  result0 <- Rlasso(y0, X0, tau_x, tau_r, method,rho=0.1, max_iter = 40, tol = 1e-6)
  rlasso = result0$x
  r2 = sum(abs(result0$r)>0.2*sqrt(log(n)/n))
  
  ###############################
  #distributed lasso on source data
  ###############################
  w_hat_A = w_hat_A_debias = NULL
  
  for(j in 1:L){
    cv.lasso <- cv.glmnet(Xs[[j]], ys[[j]], alpha=1)  # alpha=1 for Lasso
    fit_w <- glmnet(Xs[[j]], ys[[j]], alpha = 1,lambda = cv.lasso$lambda.min)
    w <- as.numeric(coef(fit_w)[-1])
    w_hat_A <- cbind(w_hat_A,w)
  }
  ###############################
  #roust_trans_lasso
  ###############################
  
  result = delta_estimation(Xs, ys,w_hat_A,  X0,   y0,  L,  p,  n, 1, 40)
  delta = result$delta
  hhat = result$hhat
  Ahat = intersect(order(hhat)[1:10],which(hhat<30))
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
  ###############################
  #plotting
  ###############################
  par(mfrow = c(5, 1), mar = c(4, 4, 2, 1), oma = c(0, 0, 2, 0))
  plot(signal_true, type = "l", col = "blue", main = "True Signal", xlab = "", ylab = "",
       cex.main = 2.5)
  # plot(y0, type = "l", col = "green", main = "True Signal", xlab = "", ylab = "",
  #        cex.main = 2.5)
  plot(lasso, type = "l", col = "magenta", main = "lasso", xlab = "", ylab = "",
       cex.main = 2.5)
  plot(rlasso, type = "l", col = "orange", main = "rlasso", xlab = "", ylab = "",
       cex.main = 2.5)
  plot(tran_lasso, type = "l", col = "purple", main = "tran_lasso", xlab = "", ylab = "",
       cex.main = 2.5)
  plot(RTL, type = "l", col =  "red", main = "robust_tran_lasso", xlab = "", ylab = "",
       cex.main = 2.5)
  
  ####
  #save
  ####
  jpeg(filename = paste0("y0_", u-1, ".jpg"), width = 800, height = 600, quality = 100)
  par(mfrow = c(1, 1), mar = c(4, 4, 2, 1), oma = c(0, 0, 2, 0))
  plot(y0, type = "l", col = "green", main = "Corrupted measurements", xlab = "", ylab = "",
       cex.main = 2.5, cex.axis = 1.7, ylim = c(-1, 1))
  dev.off()
  jpeg(filename = paste0("signal_true", ".jpg"), width = 800, height = 600, quality = 100)
  par(mfrow = c(1, 1), mar = c(4, 4, 2, 1), oma = c(0, 0, 2, 0))
  plot(signal_true, type = "l", col = "blue", main = "True Signal", xlab = "", ylab = "",
       cex.main = 2.5, cex.axis = 1.7, ylim = c(-1, 1))
  dev.off()
  jpeg(filename = paste0("lasso_", u-1, ".jpg"), width = 800, height = 600, quality = 100)
  par(mfrow = c(1, 1), mar = c(4, 4, 2, 1), oma = c(0, 0, 2, 0))
  
  plot(lasso, type = "l", col = "magenta", main = "Lasso", xlab = "", ylab = "",
       cex.main = 2.5, cex.axis = 1.7, ylim = c(-1, 1))
  dev.off()
  jpeg(filename = paste0("rlasso_", u-1, ".jpg"), width = 800, height = 600, quality = 100)
  par(mfrow = c(1, 1), mar = c(4, 4, 2, 1), oma = c(0, 0, 2, 0))
  
  plot(rlasso, type = "l", col = "orange", main = "Rlasso", xlab = "", ylab = "",
       cex.main = 2.5, cex.axis = 1.7, ylim = c(-1, 1))
  dev.off()
  jpeg(filename = paste0("tran_lasso_", u-1, ".jpg"), width = 800, height = 600, quality = 100)
  par(mfrow = c(1, 1), mar = c(4, 4, 2, 1), oma = c(0, 0, 2, 0))
  plot(tran_lasso, type = "l", col = "purple", main = "Transfer lasso", xlab = "", ylab = "",
       cex.main = 2.5, cex.axis = 1.7, ylim = c(-1, 1))
  dev.off()
  jpeg(filename = paste0("robust_tran_lasso_", u-1, ".jpg"), width = 800, height = 600, quality = 100)
  par(mfrow = c(1, 1), mar = c(4, 4, 2, 1), oma = c(0, 0, 2, 0))
  plot(RTL, type = "l", col =  "red", main = "Robust transfer lasso", xlab = "", ylab = "",
       cex.main = 2.5, cex.axis = 1.7, ylim = c(-1, 1))
  dev.off()
  #plot(signal_l1_debias, type = "l", col =  "red", main = "signal_l1_debias", xlab = "", ylab = "")
}
h0 = rep(0,L)
for(j in 1:L){
  h0[j] = sum(abs(signal_true - signal[[j]]))
}

signal_test = signal[[1]]
jpeg(filename = paste0("y_test", u-1, ".jpg"), width = 800, height = 600, quality = 100)
par(mfrow = c(1, 1), mar = c(4, 4, 2, 1), oma = c(0, 0, 2, 0))
plot(signal_test, type = "l", col = "green", main = "Signal(source data)", xlab = "", ylab = "",
     cex.main = 2.5, cex.axis = 1.7, ylim = c(-1, 1))
dev.off()