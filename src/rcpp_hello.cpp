// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::interfaces(r, cpp)]]
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
std::vector<int> traverse_cor(NumericMatrix x, float maxcor) {
  std::vector<int> elements;
  int ncol = x.ncol();
  for(int i = 0; i < ncol; i++) {
    for(int j = 0; j < ncol; j++) {
      if(i < j) {
        if(x(i, j) > maxcor && x(i, j) < 1){
          elements.push_back(i + 1); // shift index up 1 for R
        }
      }
    }
  }
  return elements;
}
