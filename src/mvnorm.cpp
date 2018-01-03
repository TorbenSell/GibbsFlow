#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

const double log2pi = std::log(2.0 * arma::datum::pi);

//' @rdname mvnpdf
//' @title Evaluate multivariate Gaussian density
//' @param x evaluation points  
//' @param mu mean vector  
//' @param sigma covariance matrix
//' @return density values
//' @export
// [[Rcpp::export]]
arma::vec mvnpdf(arma::mat x, arma::rowvec mu, arma::mat sigma){
  int n = x.n_rows;
  int dim = x.n_cols;
  arma::vec output(n);
  arma::mat rooti = arma::trans(arma::inv(arma::trimatu(arma::chol(sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(dim)/2.0) * log2pi + rootisum;
  
  for (int i = 0; i < n; i++){
    arma::vec z = rooti * arma::trans(x.row(i) - mu);
    output(i) = constants - 0.5 * arma::sum(z % z);
  }
  
  return(output);
}

//' @rdname mvnpdf_chol
//' @title Evaluate multivariate Gaussian density (with pre-computed Cholesky factor)
//' @param x evaluation points  
//' @param mu mean vector  
//' @param rooti inverse of Cholesky factor
//' @param rootisum normalizing constant term
//' @return density values
//' @export
// [[Rcpp::export]]
arma::vec mvnpdf_chol(arma::mat x, arma::rowvec mu, arma::mat rooti, double rootisum){
  int n = x.n_rows;
  int dim = x.n_cols;
  arma::vec output(n);
  double constants = -(static_cast<double>(dim)/2.0) * log2pi + rootisum;
  
  for (int i = 0; i < n; i++){
    arma::vec z = rooti * arma::trans(x.row(i) - mu);
    output(i) = constants - 0.5 * arma::sum(z % z);
  }
  
  return(output);
}

//' @rdname mvnrnd
//' @title Simulate from multivariate Gaussian distribution
//' @param n number of samples
//' @param mu mean vector  
//' @param sigma covariance matrix
//' @return samples 
//' @export
// [[Rcpp::export]]
arma::mat mvnrnd(int n, arma::rowvec mu, arma::mat sigma){
  int dim = mu.n_elem;
  arma::mat Z = arma::randn(n, dim);
  return(arma::repmat(mu,n,1) + Z * arma::chol(sigma));
}

