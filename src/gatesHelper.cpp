#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

double changeR(double x);

Eigen::VectorXd getEigenValues(NumericMatrix mat) {
  Eigen::MatrixXd M(as<Eigen::MatrixXd>(mat));
  SelfAdjointEigenSolver<Eigen::MatrixXd> es(M);
  return es.eigenvalues();
}

NumericMatrix convertToCor(NumericMatrix mat) {
  int nrow = mat.nrow();
  int ncol = mat.ncol();

  for(int row = 0; row < nrow; row++) {
    for(int col = 0; col < ncol; col++) {
      mat(row, col) = changeR(mat(row, col));
    }
  }

  return mat;
}

// [[Rcpp::export]]
double gatesHelper_(NumericMatrix mat) {
  int n_snps = mat.ncol();
  mat = convertToCor(mat);
  Eigen::VectorXd vals = getEigenValues(mat);

  NumericVector eigenvals(wrap(vals));

  double total = 0;
  double M; // from GATES paper

  for(int i = 0; i < eigenvals.length(); i++) {
    if(eigenvals(i) > 1) {
      total = total + eigenvals(i) - 1;
    }
  }

  M = n_snps - total;

  return M;
}


double changeR(double x) {
  double val = 0.2982 * pow(x, 6) - 0.0127 * pow(x, 5) + 0.0588 * pow(x, 4) +
      0.0099 * pow(x, 3) + 0.628 * pow(x, 2) - 0.0009 * x;
    return val;
}

