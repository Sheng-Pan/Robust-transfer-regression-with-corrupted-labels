#include <RcppArmadillo.h>
#include <omp.h>

using namespace Rcpp;
using namespace arma;


// Soft thresholding function
arma::vec soft_threshold(const arma::vec& x, double lambda_over_rho) {
  return sign(x) % max(abs(x) - lambda_over_rho, zeros(x.size()));
}

// ADMM solver function
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List RTL_admm_solver_cpp(const arma::vec& y0, const arma::mat& X0, const arma::vec& w_hat_A0, 
                         double lambda_beta, double lambda_e, double rho = 0.1, 
                         int max_iter = 40, double tol = 1e-6) {
  int n = y0.n_elem;
  int p = X0.n_cols;
  
  // Initialize variables
  arma::vec x = zeros<arma::vec>(p);
  arma::vec r = zeros<arma::vec>(n);
  arma::vec z = zeros<arma::vec>(p);
  arma::vec v = zeros<arma::vec>(n);
  arma::vec u = zeros<arma::vec>(p);
  arma::vec w = zeros<arma::vec>(n);
  
  // Precompute matrices
  arma::mat X0_T = X0.t();
  arma::mat I_p = eye<arma::mat>(p, p);
  arma::mat I_n = eye<arma::mat>(n, n);
  
  // ADMM iterations
  for (int j = 0; j < max_iter; ++j) {
    // Update x and r
    arma::mat A = join_vert(join_horiz(X0_T * X0 + rho * I_p, X0_T),
                            join_horiz(X0, I_n + rho * I_n));
    arma::vec b = join_vert(X0_T * (y0 - X0 * w_hat_A0) + rho * (z - u),
                            (y0 - X0 * w_hat_A0) + rho * (v - w));
    arma::vec xr_update = solve(A, b);
    arma::vec x_next = xr_update.subvec(0, p - 1);
    arma::vec r_next = xr_update.subvec(p, p + n - 1);
    
    // Update z
    arma::vec z_next = soft_threshold(x_next + u, lambda_beta / rho);
    
    // Update v
    arma::vec v_next = soft_threshold(r_next + w, lambda_e / rho);
    
    // Update dual variables u and w
    u += x_next - z_next;
    w += r_next - v_next;
    
    // Check convergence
    double primal_residual = sqrt(sum(square(x_next - z_next)) + sum(square(r_next - v_next)));
    double dual_residual = rho * (sqrt(sum(square(z_next - z)) + sqrt(sum(square(v_next - v)))));
    
    // Update variables
    x = x_next;
    r = r_next;
    z = z_next;
    v = v_next;
    u = u;
    w = w;
    
    // Check stopping criteria
    if (primal_residual < tol && dual_residual < tol) {
      break;
    }
  }
  
  // Compute final beta_hat
  arma::vec beta_hat = w_hat_A0 + x;
  
  // Return results
  return List::create(Named("beta_hat") = beta_hat, Named("r") = r);
}


// [[Rcpp::export]]
List cv_RTL_cpp(int core, arma::uvec folds,const arma::mat& X_test,const arma::vec& y_test , const arma::mat& X0, const arma::vec& y0, double ctilde, int r, int k0, int r2, const arma::vec& w_hat_A0, const arma::vec& lambda_beta_grid, const arma::vec& lambda_e_grid) {
  double best_score = datum::inf;
  List best_params = List::create(Named("lambda_beta") = NA_REAL, Named("lambda_e") = NA_REAL);
  
  // Define opf2 function
  std::function<double(const arma::vec&)> opf2;
  if (r2 > r) {
    // arma::mat X_test = Xs[test_data];
    //  arma::vec y_test = ys[test_data];
    opf2 = [&](const arma::vec& z) {
      double tau_x = z(0);
      double tau_r = z(1);
      List result1 = RTL_admm_solver_cpp(y0, X0, w_hat_A0, tau_x, tau_r, 0.1, 40, 1e-6);
      arma::vec beta_hat = as<arma::vec>(result1["beta_hat"]);
      double ab0 = mean(pow(y_test - X_test * beta_hat, 2));
      double penalty = (sum(abs(beta_hat - w_hat_A0)) < ctilde) ? 10 : 0;
      return ab0 + penalty;
    };
  } else {
    opf2 = [&](const arma::vec& z) {
      double tau_x = z(0);
      double tau_r = z(1);
      //arma::uvec folds = arma::randi<arma::uvec>(X0.n_rows, arma::distr_param(0, k0-1));
      arma::vec scores(k0, arma::fill::zeros);
      for (int i = 1; i < (k0+1); ++i) {
        arma::uvec test_idx = find(folds == i);
        arma::uvec train_idx = find(folds != i);
        arma::mat X_test = X0.rows(test_idx);
        arma::vec y_test = y0.elem(test_idx);
        arma::mat X_train = X0.rows(train_idx);
        arma::vec y_train = y0.elem(train_idx);
        List result = RTL_admm_solver_cpp(y_train, X_train, w_hat_A0, tau_x, tau_r, 0.1, 40, 1e-6);
        arma::vec beta_hat = as<arma::vec>(result["beta_hat"]);
        scores(i-1) = mean(pow(y_test - X_test * beta_hat, 2));
      }
      return mean(scores);
    };
  }
  
  // Grid search over lambda_beta_grid and lambda_e_grid
  omp_set_num_threads(core); // 设置线程数为4
#pragma omp parallel for collapse(2)
  for (double lambda_beta : lambda_beta_grid) {
    for (double lambda_e : lambda_e_grid) {
      double score;
      score = opf2({lambda_beta, lambda_e});
      //Rcout << lambda_beta << " " << lambda_e << " " << score ;
      if (score < best_score) {
        best_score = score;
        best_params["lambda_beta"] = lambda_beta;
        best_params["lambda_e"] = lambda_e;
      }
    }
  }
  
  return best_params;
}
