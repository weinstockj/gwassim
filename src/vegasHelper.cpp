#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::export]]
Eigen::MatrixXd get_chol(const Eigen::Map<Eigen::MatrixXd> & A){
  // code uses ideas from fastGP package: https://github.com/ggopalan/FastGP/blob/master/src/rcppeigen.cpp
  Eigen::MatrixXd chol = A.llt().matrixL();
  return chol;
}

// [[Rcpp::export]]
NumericMatrix rmvnorm_(const int n,
                       const NumericVector mu,
                       Eigen::Map<Eigen::MatrixXd> sigma) {
  // NumericMatrix chol_sigma(wrap(get_chol(sigma)));
  Eigen::MatrixXd chol_sigma = get_chol(sigma);

  int ncols = sigma.cols();

  NumericMatrix randvals(n, ncols);
  Eigen::MatrixXd chol_and_product_eigen;

  for(int row = 0; row < n; row++) {
    for(int col = 0; col < ncols; col++) {
      randvals(row, col) = R::rnorm(0, 1);
    }
  }

  Eigen::Map<Eigen::MatrixXd> randvals_eigen = as<Eigen::Map<Eigen::MatrixXd> >(randvals);

  chol_and_product_eigen = randvals_eigen * chol_sigma;

  return wrap(chol_and_product_eigen);

}
