#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


double soft_threshold(double a, double lambda) {
    if (a > lambda) return a - lambda;
    if (a < -lambda) return a + lambda;
    return 0.0;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List admm_lasso_cpp2(arma::mat X, arma::vec Y, double tau_x, double tau_r, double rho , int max_iter , double tol = 1e-6) {
  int n = X.n_rows; // Number of samples
  int p = X.n_cols; // Number of features
  
  // Initialize variables for ADMM algorithm
  vec q(p, fill::zeros); // Coefficient vector
  vec r(n, fill::zeros); // Residual vector
  vec z_q(p, fill::zeros); // Auxiliary variable for q
  vec z_r(n, fill::zeros); // Auxiliary variable for r
  vec u_q(p, fill::zeros); // Dual variable for q
  vec u_r(n, fill::zeros); // Dual variable for r
  
  // Precompute matrices for efficiency
  mat XT = X.t(); // Transpose of X
  mat XTX = XT * X; // X^T X matrix
  vec XTY = XT * Y; // X^T Y vector
  
  // ADMM iterations
  for (int k = 0; k < max_iter; k++) {
    // Update q and r by solving linear system
    mat A = XTX + rho * eye(p, p); // Regularized coefficient matrix
    vec b_q = XTY + rho * (z_q - u_q); // Right-hand side for q
    q = solve(A, b_q); // Solve using Armadillo's linear system solver
    
    vec b_r = Y - X * q + rho * (z_r - u_r);
    r = b_r / (1 + rho); // Update residual vector
    
    // Apply soft thresholding to auxiliary variables
    vec z_q_prev = z_q;
    vec z_r_prev = z_r;
    for (int i = 0; i < p; ++i) {
      z_q[i] = soft_threshold(q[i] + u_q[i], tau_x / rho);
    }
    for (int i = 0; i < n; ++i) {
      z_r[i] = soft_threshold(r[i] + u_r[i], tau_r / rho);
    }
    
    // Update dual variables
    u_q = u_q + q - z_q;
    u_r = u_r + r - z_r;
    
    // Check convergence criteria
    double primal_residual = norm(q - z_q, 2) + norm(r - z_r, 2);
    double dual_residual = rho * (norm(z_q - z_q_prev, 2) + norm(z_r - z_r_prev, 2));
    
    if (primal_residual < tol && dual_residual < tol) {
      Rcpp::Rcout << "Converged at iteration "<< k << std::endl;
      break;
    }
  }
  
  // Return results as a list
  return List::create(
    Named("q") = q,
    Named("r") = r,
    Named("z_q") = z_q,
    Named("z_r") = z_r
  );
}


// [[Rcpp::export]]

List delta_estimation(List Xs, List ys, double tau_r, double tau_x, arma::mat w_hat_A, arma::mat X0, arma::vec y0, int L, int p, int n, double rho = 0.1, int max_iter = 40) {
  NumericVector hhat1(L);
  std::vector<arma::vec> delta; 
  
  for (int j = 0; j < L; j++) {
    arma::vec w_hat = w_hat_A.col(j);
    arma::mat X1 = join_cols(as<arma::mat>(Xs[j]), X0);
    arma::vec Y1 = join_cols(as<arma::vec>(ys[j]), y0);
    List q_result = admm_lasso_cpp2(X1, Y1, tau_x, tau_r, rho, max_iter);
    arma::vec q = q_result["q"];
    
    
    for (uword i = 0; i < q.n_elem; ++i) {
      if (std::abs(q(i)) > 1) {
        q(i) = arma::sign(q(i)); 
      }
    }
    
    hhat1[j] = sum(abs(q - w_hat));
    delta.push_back(2 * (q - w_hat)); 
  }
  
  return List::create(Named("hhat") = hhat1, Named("delta") = delta); 
}
